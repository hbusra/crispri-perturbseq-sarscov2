#!/usr/bin/env python
"""
09 - Extension: does a perturbation do the same thing in infected and in
     bystander cells?

Writes
    results/tables/09_state_specific_effects.tsv  per element, the two profiles
    results/tables/09_program_scores_by_state.tsv program-level effects
    results/tables/09_divergent_examples.tsv      the elements that diverge most

Every per-element phenotype in 05-07 is measured on one population at a time:
the infection phenotype on all cells, the ISG phenotype on bystanders only.
That leaves an obvious question unanswered. A host factor could be required for
the cell's response to infection, in which case its knockdown signature should
appear only in infected cells; or it could set the cell's resting state, in
which case the signature should be the same either side of infection.

Method
------
For each element with enough cells in both states, a pseudobulk profile is
built separately for its infected and its bystander cells: summed counts over
cells, divided by summed host UMIs, times 1e6. The same is done for the pooled
non-targeting cells of that state, and the log2 ratio is the element's effect
in that state. Genes are restricted to those with a non-trivial NTC level in
both states, so the comparison is not dominated by genes that are simply off in
infected cells.

The headline statistic is the correlation between an element's bystander and
infected effect profiles. A high correlation means the perturbation's effect
survives infection; a low one with a large magnitude in a single state means the
effect is conditional on infection state.

Caveats, stated because they matter for how far this can be pushed:
  - infected cells carry a fraction of the host mRNA that bystanders do, so
    their profiles are noisier at the same cell count; the correlation is
    therefore attenuated, and the table reports the cell counts alongside.
  - this is a descriptive comparison of effect profiles, not a formal
    interaction test. An interaction test would need a model with a
    perturbation-by-state term, and with one infection experiment the
    between-experiment variance for it is not estimable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common as C

MIN_PER_STATE = 25
MIN_NTC_CPM = 5.0

PROGRAMS = {
    "isg": C.ISG_GENES,
    "nfkb": C.NFKB_GENES,
    "mhc_class_i": ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1",
                    "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5", "CALR", "CANX"],
    "unfolded_protein": ["HSPA5", "DDIT3", "XBP1", "ATF4", "HERPUD1", "EDEM1",
                         "SEL1L", "DNAJB9", "PDIA4", "PDIA6", "HSP90B1"],
}


def pseudobulk(X, host_umi, rows) -> np.ndarray:
    """Summed counts / summed host UMIs * 1e6 for the cells in `rows`."""
    tot = float(host_umi[rows].sum())
    if tot <= 0:
        return None
    return np.asarray(X[rows].sum(axis=0)).ravel() / tot * 1e6


def main() -> None:
    a = C.load(assigned_only=True)
    state = C.read_tbl("03_cell_state.tsv.gz").set_index("cell_barcode")
    a.obs["state"] = state.reindex(a.obs_names)["state"].to_numpy()

    host = (a.var.kind == "host").to_numpy()
    b = a[:, host]
    X = b.X.tocsr()
    genes = np.array(b.var_names)
    host_umi = a.obs.host_umi.to_numpy(float)

    inf_mask = (a.obs.state == "infected").to_numpy()
    byst_mask = (a.obs.state == "uninfected").to_numpy()
    ntc_mask = a.obs.is_ntc.to_numpy()

    ntc_inf = pseudobulk(X, host_umi, np.where(ntc_mask & inf_mask)[0])
    ntc_byst = pseudobulk(X, host_umi, np.where(ntc_mask & byst_mask)[0])
    print(
        f"NTC reference: {int((ntc_mask & inf_mask).sum())} infected, "
        f"{int((ntc_mask & byst_mask).sum())} bystander cells"
    )

    keep = (ntc_inf >= MIN_NTC_CPM) & (ntc_byst >= MIN_NTC_CPM)
    print(f"genes with NTC CPM >= {MIN_NTC_CPM} in both states: {int(keep.sum())} of {len(genes)}")

    rows, profiles = [], {}
    for element in a.obs.element.unique():
        m = (a.obs.element == element).to_numpy()
        if a.obs.is_ntc[m].iloc[0]:
            continue
        ri = np.where(m & inf_mask)[0]
        rb = np.where(m & byst_mask)[0]
        if len(ri) < MIN_PER_STATE or len(rb) < MIN_PER_STATE:
            continue
        pi = pseudobulk(X, host_umi, ri)
        pb = pseudobulk(X, host_umi, rb)
        if pi is None or pb is None:
            continue
        lfc_inf = np.log2((pi[keep] + 1) / (ntc_inf[keep] + 1))
        lfc_byst = np.log2((pb[keep] + 1) / (ntc_byst[keep] + 1))
        profiles[element] = (lfc_byst, lfc_inf)
        r = float(np.corrcoef(lfc_byst, lfc_inf)[0, 1])
        rows.append(
            dict(
                element=element,
                target=a.obs.target[m].iloc[0],
                n_infected=len(ri),
                n_bystander=len(rb),
                r_profiles=r,
                rmse_bystander=float(np.sqrt(np.mean(lfc_byst**2))),
                rmse_infected=float(np.sqrt(np.mean(lfc_inf**2))),
                n_genes_up_bystander_only=int(((lfc_byst > 1) & (lfc_inf.__abs__() < 0.5)).sum()),
                n_genes_up_infected_only=int(((lfc_inf > 1) & (lfc_byst.__abs__() < 0.5)).sum()),
            )
        )
    tab = pd.DataFrame(rows).sort_values("r_profiles")
    tab["magnitude_ratio_inf_over_byst"] = tab.rmse_infected / tab.rmse_bystander
    tab.to_csv(C.tbl("09_state_specific_effects.tsv"), sep="\t", index=False)

    print(f"\nelements with >={MIN_PER_STATE} cells in both states: {len(tab)}")
    print(
        f"profile correlation between states: median {tab.r_profiles.median():.3f}, "
        f"IQR {tab.r_profiles.quantile(.25):.3f} to {tab.r_profiles.quantile(.75):.3f}"
    )
    # noise floor: how well do two halves of the NTC pool agree with each other?
    rng = np.random.default_rng(0)
    ntc_cells_i = np.where(ntc_mask & inf_mask)[0]
    ntc_cells_b = np.where(ntc_mask & byst_mask)[0]
    rng.shuffle(ntc_cells_i)
    rng.shuffle(ntc_cells_b)
    hi, hb = np.array_split(ntc_cells_i, 2), np.array_split(ntc_cells_b, 2)
    sp_i = np.log2((pseudobulk(X, host_umi, hi[0])[keep] + 1) / (pseudobulk(X, host_umi, hi[1])[keep] + 1))
    sp_b = np.log2((pseudobulk(X, host_umi, hb[0])[keep] + 1) / (pseudobulk(X, host_umi, hb[1])[keep] + 1))
    print(
        "noise floor from splitting the NTC pool in half: "
        f"bystander RMSE {np.sqrt(np.mean(sp_b**2)):.3f}, infected RMSE {np.sqrt(np.mean(sp_i**2)):.3f} "
        "(an element's RMSE has to clear this to mean anything)"
    )

    print("\nmost divergent between states (lowest profile correlation):")
    print(
        tab.head(12)[["element", "target", "n_infected", "n_bystander", "r_profiles",
                      "rmse_bystander", "rmse_infected"]].round(3).to_string(index=False)
    )
    print("\nmost consistent across states:")
    print(
        tab.tail(8)[["element", "target", "n_infected", "n_bystander", "r_profiles",
                     "rmse_bystander", "rmse_infected"]].round(3).to_string(index=False)
    )
    tab.head(15).to_csv(C.tbl("09_divergent_examples.tsv"), sep="\t", index=False)

    # ---------------------------------------------------- program-level effects
    gi = {g: i for i, g in enumerate(genes[keep])}
    prows = []
    for element, (lfc_byst, lfc_inf) in profiles.items():
        rec = dict(element=element,
                   target=tab.loc[tab.element == element, "target"].iloc[0])
        for name, gl in PROGRAMS.items():
            idx = [gi[g] for g in gl if g in gi]
            if not idx:
                continue
            rec[f"{name}_bystander"] = float(np.mean(lfc_byst[idx]))
            rec[f"{name}_infected"] = float(np.mean(lfc_inf[idx]))
            rec[f"{name}_n_genes"] = len(idx)
        prows.append(rec)
    pt = pd.DataFrame(prows)
    pt.to_csv(C.tbl("09_program_scores_by_state.tsv"), sep="\t", index=False)

    print("\nprogram-level effects: correlation between the two states across elements")
    for name in PROGRAMS:
        cb, ci = f"{name}_bystander", f"{name}_infected"
        if cb in pt and pt[[cb, ci]].dropna().shape[0] > 5:
            sub = pt[[cb, ci]].dropna()
            print(
                f"  {name:>17}: r = {np.corrcoef(sub[cb], sub[ci])[0, 1]:+.3f} "
                f"(n = {len(sub)} elements, {int(pt[f'{name}_n_genes'].iloc[0])} genes)"
            )

    print("\nlargest interferon effects, side by side:")
    print(
        pt.reindex(pt.isg_bystander.abs().sort_values(ascending=False).index)
        .head(12)[["target", "isg_bystander", "isg_infected",
                   "mhc_class_i_bystander", "mhc_class_i_infected"]]
        .round(3).to_string(index=False)
    )


if __name__ == "__main__":
    main()
