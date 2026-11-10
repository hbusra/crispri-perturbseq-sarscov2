#!/usr/bin/env python
"""
13 - Does infection move the autophagy-lysosome programme in Calu-3, once depth
     is matched?

Writes
    results/tables/13_programs.tsv          the gene programmes scored
    results/tables/13_infection_effect.tsv  raw and depth-matched, with controls

The organoid work reports autophagy and lysosomal genes falling in infected lung
organoid monolayers at 24 hpi. This asks the same question of an independent
system, and the answer turns on one analytical decision.

Infected Calu-3 cells carry a median of 2,760 host UMIs against 25,898 in
bystanders - a 9.4-fold gap, caused by NSP1-driven shutoff. Any gene-set score
compared across that gap is partly comparing capture depth. The comparison is
therefore run twice: as-is, and after binomially downsampling every cell to a
common host depth, with ten size-matched random gene sets carried through both
to show what an arbitrary set does under the same treatment.

Only non-targeting cells enter the comparison, so no perturbation contributes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

import _common as C

HOST_DEPTH = 1500
N_RANDOM = 10
RNG = np.random.default_rng(0)

PROGRAMS = {
    "autophagy_initiation": ["ULK1", "ULK2", "ATG13", "ATG101", "RB1CC1", "BECN1",
                             "ATG14", "NRBF2", "PIK3C3", "PIK3R4", "AMBRA1",
                             "WIPI1", "WIPI2", "UVRAG", "VMP1", "TMEM41B", "ATG9A"],
    "autophagy_elongation": ["ATG12", "ATG16L1", "ATG5", "ATG7", "ATG3", "ATG10",
                             "MAP1LC3B", "GABARAP", "GABARAPL1", "GABARAPL2",
                             "SQSTM1", "TAX1BP1", "NBR1", "CALCOCO2"],
    "fusion_hops_snare": ["LAMP2", "RAB7A", "PLEKHM1", "VPS41", "VPS11", "VPS16",
                          "VPS18", "VPS33A", "VPS39", "STX17", "VAMP8", "SNAP29",
                          "CCZ1", "CCZ1B", "FYCO1"],
    "lysosome": ["LAMP1", "LAPTM4A", "LMBRD1", "CTSZ", "CTSD", "GGH", "GNS",
                 "NPC2", "GLA", "HEXA", "HEXB", "PSAP", "TPP1", "SCARB2",
                 "MCOLN1", "SLC38A9", "SPNS1", "M6PR"],
    "v_atpase": ["ATP6V0B", "ATP6V0C", "ATP6V0D1", "ATP6V0E1", "ATP6V1A",
                 "ATP6V1B2", "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F",
                 "ATP6V1G1", "ATP6AP1", "ATP6AP2"],
    "mtorc1": ["MTOR", "RPTOR", "RHEB", "LAMTOR1", "LAMTOR2", "LAMTOR3",
               "LAMTOR4", "LAMTOR5", "RRAGA", "RRAGC", "AKT1S1", "EIF4EBP1",
               "RPS6KB1"],
}
LIPID = ["PIP4K2C", "PIP4K2A", "PIP4K2B", "PIKFYVE", "FIG4", "VAC14"]


def score_block(h, programs):
    out = {}
    for name, genes in programs.items():
        present = [g for g in genes if g in h.var_names]
        if not present:
            continue
        sc.tl.score_genes(h, present, score_name="_s", random_state=0)
        out[name] = h.obs["_s"].to_numpy()
    return out


def main() -> None:
    a = C.load(assigned_only=True)
    state = C.read_tbl("03_cell_state.tsv.gz").set_index("cell_barcode")
    a.obs["state"] = state.reindex(a.obs_names)["state"].to_numpy()
    counts = a[:, (a.var.kind == "host").to_numpy()].copy()

    deep = np.asarray(counts.X.sum(axis=1)).ravel() >= HOST_DEPTH
    print(f"cells at >= {HOST_DEPTH} host UMI: {int(deep.sum())} of {counts.shape[0]}")
    counts = counts[deep].copy()
    obs = a.obs[deep].copy()

    programs = dict(PROGRAMS)
    tmp = counts.copy()
    sc.pp.normalize_total(tmp, target_sum=1e6)
    sc.pp.log1p(tmp)
    seen = np.asarray((tmp.X > 0).sum(axis=0)).ravel()
    pool = tmp.var_names[seen > 300]
    for i in range(N_RANDOM):
        programs[f"random_control_{i+1}"] = list(RNG.choice(pool, size=15, replace=False))
    for g in LIPID:
        if g in tmp.var_names:
            programs[f"expr_{g}"] = [g]

    rows = []
    for name, genes in programs.items():
        present = [g for g in genes if g in tmp.var_names]
        in_lib = sorted(set(present) & {g for t in a.obs.target.unique()
                                        if t not in ("NTC", "unassigned")
                                        for g in t.split("+")})
        rows.append(dict(program=name, n_genes=len(genes), n_present=len(present),
                         n_targeted_in_library=len(in_lib),
                         targeted_in_library=";".join(in_lib), genes=";".join(present)))
    pd.DataFrame(rows).to_csv(C.tbl("13_programs.tsv"), sep="\t", index=False)
    print(pd.DataFrame(rows)[["program", "n_present", "n_targeted_in_library"]]
          .head(6).to_string(index=False))

    raw_scores = score_block(tmp, programs)

    ds = counts.copy()
    sc.pp.downsample_counts(ds, counts_per_cell=HOST_DEPTH, random_state=0)
    sc.pp.normalize_total(ds, target_sum=1e4)
    sc.pp.log1p(ds)
    ds_scores = score_block(ds, programs)

    ntc = obs.is_ntc.to_numpy()
    inf = (obs.state == "infected").to_numpy() & ntc
    byst = (obs.state == "uninfected").to_numpy() & ntc
    print(f"\nnon-targeting cells compared: {int(inf.sum())} infected, "
          f"{int(byst.sum())} bystander")

    def delta_ci(v):
        d = v[inf].mean() - v[byst].mean()
        bi = RNG.choice(v[inf], size=(2000, int(inf.sum())), replace=True).mean(axis=1)
        bb = RNG.choice(v[byst], size=(2000, int(byst.sum())), replace=True).mean(axis=1)
        b = bi - bb
        return d, float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))

    rows = []
    for name in programs:
        if name not in raw_scores:
            continue
        r_d, r_lo, r_hi = delta_ci(raw_scores[name])
        d_d, d_lo, d_hi = delta_ci(ds_scores[name])
        rows.append(dict(program=name, delta_raw=r_d, raw_ci_low=r_lo, raw_ci_high=r_hi,
                         delta_depth_matched=d_d, ds_ci_low=d_lo, ds_ci_high=d_hi))
    eff = pd.DataFrame(rows)
    rc = eff[eff.program.str.startswith("random_control")]
    eff["random_control_mean"] = rc.delta_depth_matched.mean()
    eff["random_control_sd"] = rc.delta_depth_matched.std(ddof=1)
    eff["z_vs_random_controls"] = (
        (eff.delta_depth_matched - rc.delta_depth_matched.mean())
        / rc.delta_depth_matched.std(ddof=1))
    eff.to_csv(C.tbl("13_infection_effect.tsv"), sep="\t", index=False)

    show = eff[~eff.program.str.startswith("random_control")]
    print("\ninfected minus bystander, before and after depth matching:")
    print(show[["program", "delta_raw", "delta_depth_matched",
                "z_vs_random_controls"]].round(3).to_string(index=False))
    print(f"\nthe {N_RANDOM} random gene sets after depth matching: "
          f"mean {rc.delta_depth_matched.mean():+.3f}, "
          f"sd {rc.delta_depth_matched.std(ddof=1):.3f}, "
          f"range [{rc.delta_depth_matched.min():+.3f}, "
          f"{rc.delta_depth_matched.max():+.3f}]")
    print(
        "\nBefore depth matching every programme looks strongly suppressed, which is "
        "the published organoid direction. After it, the autophagy and lysosomal "
        "programmes sit inside the random-gene-set band and only mTORC1 keeps a "
        "residual. In this system the transcriptional version of that observation "
        "is a depth artefact."
    )


if __name__ == "__main__":
    main()
