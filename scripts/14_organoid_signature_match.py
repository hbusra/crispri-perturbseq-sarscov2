#!/usr/bin/env python
"""
14 - Which perturbations move Calu-3 cells along the organoid signatures?

Writes
    results/tables/14_signature_overview.tsv   the organoid signatures used
    results/tables/14_element_signature.tsv    per-element correlation with each
    results/tables/14_top_elements.tsv         the elements that mimic or reverse

Design
------
The organoid study (Nat Commun 2025) gives depth-matched differential-expression
profiles from ALO monolayers. This screen gives 237 perturbations in Calu-3. The
bridge is a profile-to-profile comparison: for each library element, compute its
transcriptional effect in Calu-3 and ask how it lines up with each organoid
signature.

Three decisions make the comparison honest rather than decorative.

1. **Both sides are put on a common footing.** The organoid LFCs used here are
   the downsampled ("ds") columns the organoid reanalysis provides, and the
   Calu-3 cells are likewise brought to a common host depth, so the two profiles
   are built the same way before they are compared.
2. **Perturbation effects are measured within one cell population.** Each
   element's profile is computed in bystander cells only, against non-targeting
   bystanders, so both arms share infection status and the comparison is
   perturbation against control rather than one cell state against another.
3. **The null comes from the controls.** Each non-targeting element is scored
   the same way against the other non-targeting cells, giving the distribution
   of signature correlations expected with no perturbation at all. Element
   p-values are computed against that distribution, not against zero.

Signatures (organoid side, all depth-matched):
    W3  infected vs viral-RNA-free cells, within one run, 24 h, vehicle
            - the clean infection signature
    W4  the same contrast under drug
    P6  drug vs vehicle in infected cells, 24 h
            - the therapeutic direction
    P4  the cross-run infection contrast the paper used
            - carried only to show what batch does to it
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

import _common as C

# The organoid differential-expression tables are not public and are not in this
# repository. Point ALO_LFC at the matrix of per-comparison log fold changes.
ALO = os.environ.get(
    "ALO_LFC", "/mnt/project-files/alo_reanalysis/results/04_all_gene_lfc.tsv.gz")
SIGNATURES = {
    "organoid_infection_within_run": "W3_within_run3_24h_infDMSO_VHC_vs_zero",
    "organoid_infection_under_drug": "W4_within_run3_24h_RMC_VHC_vs_zero",
    "organoid_drug_in_infected": "P6_24h_RMC_vs_DMSO_inf(run3)",
    "organoid_infection_cross_run": "P4_24h_inf(run3)_vs_uninf(run1)_DMSO [paper]",
}
HOST_DEPTH = 5000      # bystander median host UMI is ~26,000
MIN_CELLS = 30
MIN_ABS_LFC = 0.5      # genes the organoid signature actually moves
MIN_GENES = 40


def load_signatures() -> pd.DataFrame:
    if not os.path.exists(ALO):
        raise SystemExit(
            f"organoid LFC table not found at {ALO}. These tables are not public; "
            "set ALO_LFC to your own copy to run this step."
        )
    raw = pd.read_csv(ALO, sep="\t", low_memory=False)
    kind = raw.iloc[0]                      # row 0 labels each column raw / ds
    body = raw.iloc[1:].set_index(raw.columns[0])
    out = {}
    for label, col in SIGNATURES.items():
        # the depth-matched column is the ".1" partner labelled "ds"
        cand = [c for c in body.columns if c == col or c == col + ".1"]
        ds = [c for c in cand if str(kind.get(c)) == "ds"]
        use = ds[0] if ds else cand[0]
        out[label] = pd.to_numeric(body[use], errors="coerce")
        print(f"  {label}: column '{use}' ({kind.get(use)}), "
              f"{int(out[label].notna().sum())} genes")
    return pd.DataFrame(out)


def main() -> None:
    print("organoid signatures (depth-matched where available):")
    sig = load_signatures()
    rows = []
    for label in sig.columns:
        s = sig[label].dropna()
        rows.append(dict(signature=label, source_column=SIGNATURES[label],
                         n_genes=len(s), n_strong=int((s.abs() >= MIN_ABS_LFC).sum()),
                         median_abs_lfc=float(s.abs().median())))
    pd.DataFrame(rows).to_csv(C.tbl("14_signature_overview.tsv"), sep="\t", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))

    # ------------------------------------------------- Calu-3, depth-matched
    a = C.load(assigned_only=True)
    state = C.read_tbl("03_cell_state.tsv.gz").set_index("cell_barcode")
    a.obs["state"] = state.reindex(a.obs_names)["state"].to_numpy()
    h = a[:, (a.var.kind == "host").to_numpy()].copy()
    keep = (np.asarray(h.X.sum(axis=1)).ravel() >= HOST_DEPTH) & \
           (a.obs.state == "uninfected").to_numpy()
    h = h[keep].copy()
    obs = a.obs[keep].copy()
    print(f"\nbystander cells at >= {HOST_DEPTH} host UMI: {h.shape[0]} "
          f"({obs.element.nunique()} elements, {int(obs.is_ntc.sum())} control cells)")
    sc.pp.downsample_counts(h, counts_per_cell=HOST_DEPTH, random_state=0)

    X = h.X.tocsr()
    genes = np.array(h.var_names)
    is_ntc = obs.is_ntc.to_numpy()
    ntc_idx = np.where(is_ntc)[0]
    detected = np.asarray((X > 0).mean(axis=0)).ravel()
    expressed = detected >= 0.10
    print(f"genes detected in >=10% of bystanders: {int(expressed.sum())}")

    def profile(idx, ref_idx):
        """log2 fold change of mean downsampled expression, element vs reference."""
        m = np.asarray(X[idx].mean(axis=0)).ravel()
        r = np.asarray(X[ref_idx].mean(axis=0)).ravel()
        return np.log2((m + 0.01) / (r + 0.01))

    shared = {}
    for label in sig.columns:
        s = sig[label]
        strong = s[(s.abs() >= MIN_ABS_LFC) & s.notna()]
        idx = np.where(np.isin(genes, strong.index) & expressed)[0]
        shared[label] = (idx, strong.reindex(genes[idx]).to_numpy())
        print(f"  {label}: {len(idx)} genes shared, expressed and moved in the organoids")

    rows = []
    for element, sub in obs.groupby("element", observed=True):
        if len(sub) < MIN_CELLS:
            continue
        idx = np.where(obs.element.to_numpy() == element)[0]
        ref = np.setdiff1d(ntc_idx, idx)        # leave-one-out for the controls
        if len(ref) < 100:
            continue
        p = profile(idx, ref)
        row = dict(element=element, target=sub.target.iloc[0],
                   is_ntc=bool(sub.is_ntc.iloc[0]), n_cells=int(len(sub)))
        for label, (gidx, slfc) in shared.items():
            if len(gidx) < MIN_GENES:
                continue
            rho = stats.spearmanr(p[gidx], slfc).statistic
            row[f"rho_{label}"] = float(rho)
        rows.append(row)
    el = pd.DataFrame(rows)

    # ----------------------------------------- control-based null and testing
    for label in sig.columns:
        col = f"rho_{label}"
        if col not in el:
            continue
        null = el.loc[el.is_ntc, col].dropna()
        mu, sd = float(null.mean()), float(null.std(ddof=1))
        el[f"z_{label}"] = (el[col] - mu) / sd
        el[f"p_{label}"] = 2 * stats.norm.sf(el[f"z_{label}"].abs())
        q = np.full(len(el), np.nan)
        m = (~el.is_ntc).to_numpy()
        q[m] = C.bh_fdr(el.loc[m, f"p_{label}"].to_numpy())
        el[f"q_{label}"] = q
        print(f"\n{label}: control correlations {mu:+.3f} +- {sd:.3f} "
              f"(n = {len(null)} non-targeting elements)")
        t = el[~el.is_ntc]
        hits = t[t[f"q_{label}"] < 0.1]
        print(f"  {len(hits)} of {len(t)} elements differ from the control band at q < 0.1")
        if len(hits):
            print("  most like the organoid signature:")
            print(hits.nlargest(6, col)[["target", "n_cells", col, f"q_{label}"]]
                  .round(3).to_string(index=False))
            print("  most opposite:")
            print(hits.nsmallest(6, col)[["target", "n_cells", col, f"q_{label}"]]
                  .round(3).to_string(index=False))
    el.to_csv(C.tbl("14_element_signature.tsv"), sep="\t", index=False)

    # --------------------------------------------------------- the headline
    inf_col = "rho_organoid_infection_within_run"
    drug_col = "rho_organoid_drug_in_infected"
    t = el[~el.is_ntc].copy()
    if inf_col in t and drug_col in t:
        print(f"\nthe two organoid axes are correlated across elements at "
              f"r = {np.corrcoef(t[inf_col], t[drug_col])[0,1]:.3f} "
              "(a perturbation that mimics infection tends to oppose the drug)")
        t["mimics_infection"] = t[f"q_organoid_infection_within_run"] < 0.1
        top = pd.concat([t.nlargest(10, inf_col), t.nsmallest(10, inf_col)])
        top.to_csv(C.tbl("14_top_elements.tsv"), sep="\t", index=False)

    # autophagy-module elements specifically
    MODULE = ["ATG14", "ATG9A", "BECN1", "RB1CC1", "UVRAG", "VMP1", "TMEM41B",
              "RAB7A", "VPS11", "VPS16", "VPS33A", "VPS39", "CCZ1", "CCZ1B",
              "FYCO1", "ATP6AP1", "ATP6V1A", "ATP6V1B2", "M6PR", "SPNS1",
              "TMEM106B", "MFSD12", "SIGMAR1", "VAC14"]
    mod = el[el.target.isin(MODULE)]
    if len(mod) and inf_col in mod:
        print("\nthe autophagy-endolysosomal module, against the organoid infection "
              "signature:")
        print(mod.sort_values(inf_col)[
            ["target", "n_cells", inf_col, f"q_organoid_infection_within_run"]
        ].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
