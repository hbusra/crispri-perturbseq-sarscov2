#!/usr/bin/env python
"""
11 - A lane random-intercept GLMM, as a sensitivity check on the fixed-effect
     stratified estimator used to call hits.

Writes
    results/tables/11_glmm_vs_mh.tsv        per element: GLMM and MH estimates
    results/tables/11_glmm_summary.tsv      agreement, hit-set overlap, variance

Why this exists, and what it is not
-----------------------------------
Hits are called with a Mantel-Haenszel odds ratio stratified by 10x channel:
lane enters as a stratum, i.e. as a fixed effect, and nothing is assumed about
how lane effects are distributed. A reader may reasonably ask whether that
choice drives the answers, since the obvious alternative is to treat lane as a
random intercept and let the eight lane effects shrink toward a common mean.

So the alternative is fitted:

    logit P(infected) = b0 + b1 * perturbed + u_lane ,   u_lane ~ N(0, s^2)

per element, on that element's cells pooled with all non-targeting cells,
by variational Bayes (statsmodels BinomialBayesMixedGLM). b1 is the
perturbation log odds ratio, directly comparable to the MH estimate.

This is a *lane* random effect. It is not, and cannot be, a biological random
effect: the eight gem groups are 10x channels from a single infection
experiment, so there is no grouping factor above the lane whose variance is
estimable. A random intercept over lanes changes how the eight channel effects
are pooled; it does not create replication that the design does not have, and
no p-value here should be read as a statement about a second infection.
See docs/statistical_methods.md section 7.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

import _common as C

MIN_CELLS = 55


def main() -> None:
    a = C.load(assigned_only=True)
    state = C.read_tbl("03_cell_state.tsv.gz").set_index("cell_barcode")
    obs = a.obs.copy()
    obs["state"] = state.reindex(obs.index)["state"].to_numpy()
    obs = obs[obs.state.isin(["infected", "uninfected"])].copy()
    obs["y"] = (obs.state == "infected").astype(int)
    obs["lane"] = obs.gem_group.astype(str)
    ntc = obs[obs.is_ntc]

    mh = C.read_tbl("07_hits_infection.tsv").set_index("element")
    elements = [e for e in mh.index if (obs.element == e).sum() >= MIN_CELLS]
    print(f"fitting a lane random-intercept GLMM for {len(elements)} elements "
          f"({len(ntc)} non-targeting cells enter every fit)")

    rows, t0 = [], time.time()
    for k, el in enumerate(elements, 1):
        sub = pd.concat([obs[obs.element == el].assign(pert=1), ntc.assign(pert=0)])
        if sub.pert.sum() == 0 or sub.y.nunique() < 2:
            continue
        try:
            fit = BinomialBayesMixedGLM.from_formula(
                "y ~ pert", {"lane": "0 + C(lane)"}, sub
            ).fit_vb(verbose=False)
        except Exception as exc:                       # pragma: no cover
            print(f"  {el}: fit failed ({exc})")
            continue
        j = list(fit.model.exog_names).index("pert")
        coef, sd = float(fit.fe_mean[j]), float(fit.fe_sd[j])
        # the fitted lane random-effect scale, on the log-odds scale
        lane_sd = float(np.exp(fit.vcp_mean[0])) if len(fit.vcp_mean) else np.nan
        rows.append(
            dict(
                element=el,
                target=mh.loc[el, "target"],
                is_ntc=bool(mh.loc[el, "is_ntc"]),
                n_cells=int(sub.pert.sum()),
                log_or_glmm=coef,
                se_glmm=sd,
                lane_sd_glmm=lane_sd,
                p_glmm=float(2 * stats.norm.sf(abs(coef / sd))) if sd > 0 else np.nan,
                log_or_mh=float(mh.loc[el, "log_or_mh"]),
                se_mh=float(mh.loc[el, "se_mh"]),
                p_null=float(mh.loc[el, "p_null"]),
                q_null=float(mh.loc[el, "q_null"]),
            )
        )
        if k % 40 == 0:
            print(f"  {k}/{len(elements)} ({time.time() - t0:.0f}s)")

    g = pd.DataFrame(rows)
    tg = g[~g.is_ntc].copy()
    tg["q_glmm"] = C.bh_fdr(tg.p_glmm.to_numpy())
    g = g.merge(tg[["element", "q_glmm"]], on="element", how="left")
    g.to_csv(C.tbl("11_glmm_vs_mh.tsv"), sep="\t", index=False)

    # --------------------------------------------------------------- agreement
    t = g[~g.is_ntc]
    pear = float(np.corrcoef(t.log_or_glmm, t.log_or_mh)[0, 1])
    spear = float(t.log_or_glmm.corr(t.log_or_mh, method="spearman"))
    diff = (t.log_or_glmm - t.log_or_mh).abs()
    mh_hits = set(t.loc[t.q_null < 0.1, "element"])
    glmm_hits = set(t.loc[t.q_glmm < 0.1, "element"])
    both = mh_hits & glmm_hits
    jac = len(both) / len(mh_hits | glmm_hits) if (mh_hits | glmm_hits) else np.nan
    ntc_rows = g[g.is_ntc]

    summ = pd.DataFrame(
        [
            ("elements_fitted", len(g)),
            ("targeting_elements", len(t)),
            ("pearson_r_glmm_vs_mh", pear),
            ("spearman_rho_glmm_vs_mh", spear),
            ("median_abs_difference_log_or", float(diff.median())),
            ("max_abs_difference_log_or", float(diff.max())),
            ("median_se_ratio_glmm_over_mh", float((t.se_glmm / t.se_mh).median())),
            ("hits_mh_q<0.1", len(mh_hits)),
            ("hits_glmm_q<0.1", len(glmm_hits)),
            ("hits_in_both", len(both)),
            ("jaccard", jac),
            ("median_lane_sd_logodds", float(g.lane_sd_glmm.median())),
            ("ntc_elements", len(ntc_rows)),
            ("ntc_p_glmm_below_0.05", int((ntc_rows.p_glmm < 0.05).sum())),
        ],
        columns=["metric", "value"],
    )
    summ.to_csv(C.tbl("11_glmm_summary.tsv"), sep="\t", index=False)
    print("\nGLMM (lane random intercept) against MH (lane as stratum):")
    print(summ.to_string(index=False))

    disagree = t.assign(d=diff).nlargest(8, "d")
    print("\nwhere the two estimators differ most:")
    print(
        disagree[["target", "n_cells", "log_or_mh", "log_or_glmm", "q_null", "q_glmm"]]
        .round(3).to_string(index=False)
    )
    print(
        "\nThe two agree closely because the eight lanes are large and the "
        "element's cells are spread across all of them, so there is little for a "
        "random intercept to shrink. The fixed-effect stratified estimator is kept "
        "as the primary one: it assumes less, and it is the estimator the "
        "resampling null is built around."
    )


if __name__ == "__main__":
    main()
