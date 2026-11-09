#!/usr/bin/env python
"""
02 - Cell-level QC of the assigned population, reported per 10x channel.

Writes
    results/tables/02_qc_per_cell_summary.tsv     overall and per-gem-group QC
    results/tables/02_qc_by_gem_group.tsv         lane effects on depth/viral load
    results/tables/02_cells_per_element.tsv       representation per library element

What this script is for
-----------------------
Two things matter downstream and both are checked here rather than assumed:

1. Depth is not constant across lanes. Any per-element phenotype is at risk of
   confounding if elements are unevenly distributed across lanes, so we report
   the lane composition of every element and the depth spread between lanes.

2. The viral fraction of a cell's library is large (up to most of the library).
   That means "total UMI" is not a depth covariate: it is partly the phenotype.
   Every downstream model therefore uses host-only UMI as the offset, and this
   script quantifies how different the two denominators are.

No cells are filtered here. The authors' filtered matrix is already QC'd; this
script documents the state of the analysis set and flags anything that would
invalidate a model, rather than re-filtering.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common as C


def q(s, p):
    return float(np.nanpercentile(s, p))


def main() -> None:
    a = C.load(assigned_only=False)
    obs = a.obs
    asg = obs[obs.cell_set == "assigned"]
    wt = obs[obs.cell_set == "wt_candidate"]

    # ------------------------------------------------------------- overall QC
    rows = []
    for label, d in [("assigned", asg), ("wt_candidate", wt)]:
        rows.append(
            dict(
                population=label,
                n_cells=len(d),
                host_umi_p25=q(d.host_umi, 25),
                host_umi_median=q(d.host_umi, 50),
                host_umi_p75=q(d.host_umi, 75),
                total_umi_median=q(d.total_umi, 50),
                n_genes_median=q(d.n_genes, 50),
                pct_mito_median=q(d.pct_mito, 50),
                viral_umi_median=q(d.viral_umi, 50),
                viral_frac_median=q(d.viral_frac, 50),
                viral_frac_p90=q(d.viral_frac, 90),
                lenti_umi_median=q(d.lenti_umi, 50),
                frac_lenti_zero=float((d.lenti_umi == 0).mean()),
            )
        )
    summ = pd.DataFrame(rows)
    summ.to_csv(C.tbl("02_qc_per_cell_summary.tsv"), sep="\t", index=False)
    print(summ.to_string(index=False))

    # ------------------------------------------------- how different is total vs host
    # If total UMI were a valid depth covariate it would track host UMI; it does
    # not, because viral reads dominate a subset of cells.
    r_all = np.corrcoef(asg.total_umi, asg.host_umi)[0, 1]
    hi = asg[asg.viral_frac > 0.5]
    r_hi = np.corrcoef(hi.total_umi, hi.host_umi)[0, 1] if len(hi) > 10 else np.nan
    print(
        f"\ncorr(total_umi, host_umi): all assigned {r_all:.3f}; "
        f"cells with viral_frac>0.5 (n={len(hi)}) {r_hi:.3f}"
    )

    # ------------------------------------------------------------ lane effects
    g = (
        asg.groupby("gem_group", observed=True)
        .agg(
            n_cells=("host_umi", "size"),
            host_umi_median=("host_umi", "median"),
            n_genes_median=("n_genes", "median"),
            pct_mito_median=("pct_mito", "median"),
            viral_frac_median=("viral_frac", "median"),
            viral_umi_median=("viral_umi", "median"),
            n_elements=("element", "nunique"),
        )
        .reset_index()
    )
    g.to_csv(C.tbl("02_qc_by_gem_group.tsv"), sep="\t", index=False)
    print("\nper-lane:\n", g.to_string(index=False))
    spread = g.host_umi_median.max() / g.host_umi_median.min()
    print(f"\nlane median host-UMI max/min = {spread:.2f}x")

    # --------------------------------------------------- element representation
    ct = (
        asg.groupby("element", observed=True)
        .agg(
            target=("target", "first"),
            is_ntc=("is_ntc", "first"),
            n_target_genes=("n_target_genes", "first"),
            n_cells=("host_umi", "size"),
            n_lanes=("gem_group", "nunique"),
            host_umi_median=("host_umi", "median"),
            viral_frac_median=("viral_frac", "median"),
        )
        .reset_index()
        .sort_values("n_cells")
    )
    # lane imbalance: largest share of an element's cells coming from one lane
    share = (
        asg.groupby(["element", "gem_group"], observed=True).size()
        .groupby(level=0).apply(lambda s: s.max() / s.sum())
        .rename("max_lane_share")
    )
    ct = ct.merge(share, left_on="element", right_index=True, how="left")
    ct.to_csv(C.tbl("02_cells_per_element.tsv"), sep="\t", index=False)

    print(
        f"\nelements: {len(ct)}  (NTC {int(ct.is_ntc.sum())}, "
        f"single-gene {int(((~ct.is_ntc) & (ct.n_target_genes == 1)).sum())}, "
        f"multi-gene {int((ct.n_target_genes > 1).sum())})"
    )
    print(
        f"cells per element: min {ct.n_cells.min()}, median {int(ct.n_cells.median())}, "
        f"max {ct.n_cells.max()};  {(ct.n_cells < 55).sum()} elements under 55 cells"
    )
    print(f"all elements present in all 8 lanes: {bool((ct.n_lanes == 8).all())}")
    print(f"worst single-lane share for any element: {ct.max_lane_share.max():.2f}")


if __name__ == "__main__":
    main()
