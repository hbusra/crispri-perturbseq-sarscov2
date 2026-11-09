#!/usr/bin/env python
"""
08 - Sensitivity to perturbation representation.

Writes
    results/tables/08_downsample_curve.tsv      effect and detection vs cells per element
    results/tables/08_detectable_effect.tsv     smallest detectable effect at each cell count
    results/tables/08_underpowered_elements.tsv elements where the screen cannot answer

This is deliberately *not* a formal power analysis. There is one infection
experiment here, so there is no replicate structure to estimate a between-
experiment variance component from, and any "power" number would be a
statement about cell-level sampling noise dressed up as a statement about the
experiment. What can be said honestly is narrower and still useful:

  - how the estimated effect and its uncertainty behave as an element's cell
    count is reduced, holding the biology fixed (subsampling real cells);
  - what effect size the NTC-based null can separate from zero at a given cell
    count, which is a property of the screen's coverage;
  - which library elements have too few cells for either phenotype, so that
    their absence from the hit list is not read as evidence of no effect.

The last point is the one that matters for interpreting a screen: 47 of 237
elements have fewer than 55 cells, and the paper drops 48 for the same reason.
Those elements are mostly essential genes, i.e. the ones whose knockdown kills
cells, so their under-representation is informative rather than random.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common as C
from _hitcalling_core import mh_log_or, lane_counts

DRAWS = 400
SUBSAMPLE_SIZES = [25, 40, 60, 80, 120, 160]
RNG = np.random.default_rng(0)


def main() -> None:
    d = C.read_tbl("03_cell_state.tsv.gz")
    d["gem_group"] = d.gem_group.astype(str)
    lanes = sorted(d.gem_group.unique())
    ntc = d[d.is_ntc]
    c_l, d_l = lane_counts(ntc, lanes, "state")
    p_l = np.divide(c_l, c_l + d_l, out=np.zeros(len(lanes)), where=(c_l + d_l) > 0)

    # Pick real hits spanning the whole effect-size range, not just the biggest
    # ones: a subsampling curve built only on ACE2-sized effects says nothing
    # about the elements whose detection is actually marginal.
    hits = C.read_tbl("07_hits_infection.tsv")
    sig = hits[(~hits.is_ntc) & (hits.q_perm < 0.1)].copy()
    sig["abs_eff"] = sig.log_or_mh.abs()
    sig = sig.sort_values("abs_eff")
    picks = np.unique(np.linspace(0, len(sig) - 1, 8).round().astype(int))
    strong = sig.iloc[picks]
    print("elements used for the subsampling curve (weakest to strongest):")
    print(strong[["target", "n_cells", "log_or_mh", "q_perm"]].round(3).to_string(index=False))

    # ---------------------------------------- subsample real cells of real hits
    rows = []
    for _, h in strong.iterrows():
        sub = d[d.element == h.element]
        for k in SUBSAMPLE_SIZES + [len(sub)]:
            if k > len(sub):
                continue
            eff, hit = [], 0
            for _ in range(60):
                take = sub.sample(k, random_state=int(RNG.integers(1e9)))
                a_l, b_l = lane_counts(take, lanes, "state")
                lo, _se = mh_log_or(a_l, b_l, c_l, d_l)
                if not np.isfinite(lo):
                    continue
                n_l = a_l + b_l
                a_null = RNG.binomial(np.repeat(n_l[None, :], DRAWS, axis=0), p_l[None, :])
                b_null = n_l[None, :] - a_null
                with np.errstate(divide="ignore", invalid="ignore"):
                    null = np.log(
                        np.sum(a_null * d_l / (n_l + c_l + d_l), axis=1)
                        / np.sum(b_null * c_l / (n_l + c_l + d_l), axis=1)
                    )
                p = (1 + np.sum(np.abs(null[np.isfinite(null)]) >= abs(lo))) / (DRAWS + 1)
                eff.append(lo)
                hit += int(p < 0.05)
            if eff:
                rows.append(
                    dict(element=h.element, target=h.target, cells=k,
                         mean_log_or=float(np.mean(eff)),
                         sd_log_or=float(np.std(eff)),
                         detection_rate=hit / len(eff),
                         full_log_or=h.log_or_mh)
                )
    curve = pd.DataFrame(rows)
    curve.to_csv(C.tbl("08_downsample_curve.tsv"), sep="\t", index=False)
    print("\nsubsampling real hits: detection rate at p<0.05, by cells retained")
    print(
        curve[curve.cells.isin(SUBSAMPLE_SIZES)]
        .pivot_table(index="cells", columns="target", values="detection_rate")
        .round(2).to_string()
    )
    print("\nand the effect estimate itself (mean log odds ratio):")
    print(
        curve[curve.cells.isin(SUBSAMPLE_SIZES)]
        .pivot_table(index="cells", columns="target", values="mean_log_or")
        .round(2).to_string()
    )

    # -------------------------------- smallest detectable effect per cell count
    # The null's spread at a given cell count sets the resolution of the screen.
    rows = []
    for k in [25, 40, 60, 80, 120, 160, 240]:
        n_l = np.full(len(lanes), k / len(lanes))
        n_l = np.floor(n_l).astype(int)
        n_l[: k - n_l.sum()] += 1
        a_null = RNG.binomial(np.repeat(n_l[None, :], 4000, axis=0), p_l[None, :])
        b_null = n_l[None, :] - a_null
        with np.errstate(divide="ignore", invalid="ignore"):
            null = np.log(
                np.sum(a_null * d_l / (n_l + c_l + d_l), axis=1)
                / np.sum(b_null * c_l / (n_l + c_l + d_l), axis=1)
            )
        null = null[np.isfinite(null)]
        rows.append(
            dict(cells=k,
                 null_sd=float(np.std(null)),
                 abs_log_or_at_p05=float(np.percentile(np.abs(null), 95)),
                 fold_change_at_p05=float(np.exp(np.percentile(np.abs(null), 95)))))
    det = pd.DataFrame(rows)
    det.to_csv(C.tbl("08_detectable_effect.tsv"), sep="\t", index=False)
    print("\nresolution of the screen, from the NTC null alone:")
    print(det.round(3).to_string(index=False))

    # ------------------------------------------------- underpowered elements
    rep = C.read_tbl("04_element_representation.tsv")
    rep["testable_infection"] = rep.n_cells >= 55
    rep["testable_burden"] = rep.n_infected >= 20
    rep["testable_isg"] = rep.n_bystander >= 30
    under = rep[~(rep.testable_infection & rep.testable_burden & rep.testable_isg)]
    under.to_csv(C.tbl("08_underpowered_elements.tsv"), sep="\t", index=False)
    print(
        f"\nelements of {len(rep)} not testable for at least one phenotype: {len(under)}"
        f"\n  infection (>=55 cells)      : {int((~rep.testable_infection).sum())} fail"
        f"\n  burden   (>=20 infected)    : {int((~rep.testable_burden).sum())} fail"
        f"\n  ISG      (>=30 bystanders)  : {int((~rep.testable_isg).sum())} fail"
    )
    print("\nthe 15 worst-represented elements (absence from a hit list means nothing here):")
    print(
        under.nsmallest(15, "n_cells")[
            ["element", "target", "n_cells", "n_infected", "n_bystander"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
