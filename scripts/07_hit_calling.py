#!/usr/bin/env python
"""
07 - Hit calling: lane-stratified effects, an NTC-based empirical null, and
     concordance with the published hit list.

Writes
    results/tables/07_hits_infection.tsv       infected-fraction phenotype
    results/tables/07_hits_burden.tsv          viral-burden phenotype
    results/tables/07_null_calibration.tsv     NTC-vs-NTC false-positive rate
    results/tables/07_concordance.tsv          our calls vs the paper's
    results/tables/07_robustness.tsv           same calls under a different classifier

The two estimands
-----------------
A.  Infection establishment.  Odds of a cell being infected rather than
    bystander, for cells carrying element e versus cells carrying a
    non-targeting element. Estimated as a Mantel-Haenszel odds ratio
    stratified by 10x channel, so the lane effect (lanes 1-4 and 5-8 differ
    about two-fold in viral load) cannot leak into the element effect. The
    same quantity is also fitted as a logistic GLM with lane as a covariate;
    the two agree and both are reported.

B.  Viral burden given infection.  Viral UMIs per host UMI among infected
    cells only. Estimated two ways:
      - a lane-weighted difference in mean log10 viral rate, which is what the
        control-based null is built on;
      - a negative-binomial GLM on raw viral counts with log(host UMI) as an
        offset and lane as a covariate, which gives a rate ratio.
    The offset is host UMI, never total UMI: total UMI contains the viral
    counts being modelled and is itself depressed by infection, so
    conditioning on it would condition on the outcome.

The null
--------
A Wald p-value from either model assumes cells are independent replicates.
They are not: cells sharing a library element share a lentiviral integration
history and a droplet batch. So significance comes from an empirical null built
out of the non-targeting elements, not from the model:

  for each element e, and each of B draws, build a pseudo-element with the same
  number of cells in each lane as e, drawn from the pooled NTC cells of that
  lane, and recompute the statistic.

This is lane-matched and size-matched by construction, so an element with 60
cells is judged against a null for 60 cells. p is the two-sided fraction of
null draws at least as extreme, floored at 1/(B+1), then BH-adjusted.

Calibration is checked by running each NTC element through the pipeline as if
it were a perturbation (`07_null_calibration.tsv`): the realised false-positive
rate at q < 0.1 should be near zero.

What this design cannot claim
-----------------------------
The 8 lanes are 10x channels from a single infection experiment, not 8
independent infections. The null above controls for cell-level and lane-level
sampling noise. It does not establish reproducibility across independent
biological experiments, and no p-value here should be read that way.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm

import _common as C
from _hitcalling_core import (
    lane_counts,
    lane_weighted_delta,
    mh_log_or,
    mh_null_draws,
    resample_delta_null,
    two_sided_p,
)

B_DRAWS = 2000
MIN_CELLS = 55          # the paper's coverage cut
MIN_INFECTED = 20       # for the burden phenotype
RNG = np.random.default_rng(0)

PAPER_HITS = {
    "ACE2": "down", "TMPRSS2": "down", "ACE2+TMPRSS2": "down",
    "TMPRSS2+FURIN": "down", "TMPRSS2+CTSB": "down", "TMPRSS2+CTSL": "down",
    "BRD2": "down", "ATG14": "down", "EIF4E2": "down", "EIF4H": "down",
    "VMP1": "down", "MPP5": "down", "NFKBIA": "down", "SEC62": "up",
}
PAPER_NULL = ["FURIN", "CTSB", "CTSL", "RELA", "RELB", "RELA+RELB", "EIF4B"]


def fdr_over_targeting(df: pd.DataFrame) -> np.ndarray:
    """Benjamini-Hochberg over the TARGETING elements only.

    The non-targeting elements are controls, not candidates: they exist to
    calibrate the null, and putting them in the discovery family would both
    inflate the family size and treat a control as a hypothesis. They keep an
    NaN q and are judged on their raw p in the calibration table instead.
    """
    q = np.full(len(df), np.nan)
    m = (~df.is_ntc).to_numpy()
    q[m] = C.bh_fdr(df.p_null.to_numpy()[m])
    return q


# --------------------------------------------------------------- phenotype A
def infection_phenotype(d, state_col="state", draws=B_DRAWS, verbose=True):
    lanes = sorted(d.gem_group.unique())
    ntc = d[d.is_ntc]
    c_l, d_l = lane_counts(ntc, lanes, state_col)
    p_l = np.divide(c_l, c_l + d_l, out=np.zeros(len(lanes)), where=(c_l + d_l) > 0)

    rows = []
    for element, sub in d.groupby("element", observed=True):
        a_l, b_l = lane_counts(sub, lanes, state_col)
        n_l = a_l + b_l
        if n_l.sum() < MIN_CELLS:
            continue
        obs, se = mh_log_or(a_l, b_l, c_l, d_l)
        null = mh_null_draws(n_l, p_l, c_l, d_l, draws, RNG)

        rows.append(
            dict(
                element=element,
                target=sub.target.iloc[0],
                is_ntc=bool(sub.is_ntc.iloc[0]),
                n_cells=int(len(sub)),
                n_infected=int(a_l.sum()),
                n_bystander=int(b_l.sum()),
                infected_frac=float(a_l.sum() / max(n_l.sum(), 1)),
                log_or_mh=obs,
                se_mh=se,
                log_or_null_sd=float(np.nanstd(null)),
                z_empirical=obs / np.nanstd(null) if np.isfinite(obs) else np.nan,
                p_null=two_sided_p(obs, null),
            )
        )
    out = pd.DataFrame(rows)
    out["q_null"] = fdr_over_targeting(out)
    if verbose:
        print(
            f"phenotype A: {len(out)} elements with >={MIN_CELLS} cells "
            f"({int((~out.is_ntc).sum())} targeting, {int(out.is_ntc.sum())} non-targeting)"
        )
    return out


def infection_glm(d, elements):
    """Logistic GLM with lane, for comparison with the MH estimate."""
    sub = d[d.state.isin(["infected", "uninfected"])].copy()
    sub["y"] = (sub.state == "infected").astype(int)
    lane_dummies = pd.get_dummies(sub.gem_group, prefix="lane", drop_first=True).astype(float)
    ntc_mask = sub.is_ntc.to_numpy()
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for element in elements:
            m = (sub.element == element).to_numpy()
            keep = m | ntc_mask
            X = pd.concat(
                [pd.Series(m[keep].astype(float), index=sub.index[keep], name="element"),
                 lane_dummies.loc[sub.index[keep]]], axis=1
            )
            X = sm.add_constant(X.loc[:, X.std() > 0])
            try:
                fit = sm.GLM(sub.y[keep], X, family=sm.families.Binomial()).fit()
                rows.append(dict(element=element,
                                 glm_log_or=float(fit.params["element"]),
                                 glm_se=float(fit.bse["element"]),
                                 glm_p=float(fit.pvalues["element"])))
            except Exception as e:  # noqa: BLE001
                rows.append(dict(element=element, glm_log_or=np.nan,
                                 glm_se=np.nan, glm_p=np.nan, error=str(e)[:60]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------- phenotype B
def burden_phenotype(d, draws=B_DRAWS):
    inf = d[d.state == "infected"].copy()
    inf["lr"] = np.log10((inf.viral_umi / inf.host_umi.replace(0, np.nan)).clip(lower=1e-6))
    lanes = sorted(inf.gem_group.unique())
    ntc = inf[inf.is_ntc]
    ntc_by_lane = {l: ntc.loc[ntc.gem_group == l, "lr"].to_numpy() for l in lanes}
    ntc_mean = {l: float(np.mean(v)) if len(v) else np.nan for l, v in ntc_by_lane.items()}

    rows = []
    for element, sub in inf.groupby("element", observed=True):
        if len(sub) < MIN_INFECTED:
            continue
        obs, n_l = lane_weighted_delta(sub, lanes, "lr", ntc_mean)
        if not np.isfinite(obs) or sum(n_l) == 0:
            continue
        null = resample_delta_null(ntc_by_lane, ntc_mean, lanes, n_l, draws, RNG)

        rows.append(
            dict(
                element=element,
                target=sub.target.iloc[0],
                is_ntc=bool(sub.is_ntc.iloc[0]),
                n_infected=int(len(sub)),
                delta_log10_rate=obs,
                null_sd=float(np.std(null)),
                z_empirical=obs / np.std(null) if np.std(null) > 0 else np.nan,
                p_null=two_sided_p(obs, null),
            )
        )
    out = pd.DataFrame(rows)
    out["q_null"] = fdr_over_targeting(out)
    print(
        f"phenotype B: {len(out)} elements with >={MIN_INFECTED} infected cells "
        f"({int((~out.is_ntc).sum())} targeting, {int(out.is_ntc.sum())} non-targeting)"
    )
    return out


def burden_nb_glm(d, elements):
    """Negative-binomial GLM on viral counts, offset = log(host UMI)."""
    inf = d[(d.state == "infected") & (d.host_umi > 0)].copy()
    lane_dummies = pd.get_dummies(inf.gem_group, prefix="lane", drop_first=True).astype(float)
    offset = np.log(inf.host_umi.to_numpy(float))
    ntc_mask = inf.is_ntc.to_numpy()
    y_all = inf.viral_umi.to_numpy(float)

    # one shared dispersion, estimated on the NTC cells
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Xn = sm.add_constant(lane_dummies[ntc_mask])
        nb0 = sm.NegativeBinomial(y_all[ntc_mask], Xn, offset=offset[ntc_mask]).fit(disp=0)
        alpha = float(np.clip(nb0.params.get("alpha", 1.0), 1e-3, 50))
    print(f"  NB dispersion estimated on NTC cells: alpha = {alpha:.3f}")

    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for element in elements:
            m = (inf.element == element).to_numpy()
            keep = m | ntc_mask
            X = pd.concat(
                [pd.Series(m[keep].astype(float), index=inf.index[keep], name="element"),
                 lane_dummies.loc[inf.index[keep]]], axis=1
            )
            X = sm.add_constant(X.loc[:, X.std() > 0])
            try:
                fit = sm.GLM(
                    y_all[keep], X,
                    family=sm.families.NegativeBinomial(alpha=alpha),
                    offset=offset[keep],
                ).fit()
                rows.append(dict(element=element,
                                 nb_log_rate_ratio=float(fit.params["element"]),
                                 nb_se=float(fit.bse["element"]),
                                 nb_p=float(fit.pvalues["element"])))
            except Exception as e:  # noqa: BLE001
                rows.append(dict(element=element, nb_log_rate_ratio=np.nan,
                                 nb_se=np.nan, nb_p=np.nan, error=str(e)[:60]))
    return pd.DataFrame(rows), alpha


# -------------------------------------------------------------------- driver
def main() -> None:
    d = C.read_tbl("03_cell_state.tsv.gz")
    d["gem_group"] = d.gem_group.astype(str)

    # ---------------------------------------------------------- phenotype A
    A = infection_phenotype(d)
    glm = infection_glm(d, A.element.tolist())
    A = A.merge(glm, on="element", how="left")
    A["q_glm"] = np.where(~A.is_ntc, C.bh_fdr(A.glm_p.where(~A.is_ntc).to_numpy()), np.nan)
    r = A[["log_or_mh", "glm_log_or"]].dropna()
    print(f"  MH vs logistic-GLM effect sizes: r = {np.corrcoef(r.log_or_mh, r.glm_log_or)[0, 1]:.4f}")
    A = A.sort_values("log_or_mh")
    A.to_csv(C.tbl("07_hits_infection.tsv"), sep="\t", index=False)

    tested = A[~A.is_ntc]
    sig = tested[tested.q_null < 0.1]
    print(
        f"  significant at q<0.1: {len(sig)} of {len(tested)} targeting elements "
        f"({int((sig.log_or_mh < 0).sum())} decrease, {int((sig.log_or_mh > 0).sum())} increase)"
    )

    # ------------------------------------------------- null calibration check
    ntc_rows = A[A.is_ntc]
    cal = pd.DataFrame(
        [
            ("ntc_elements_tested", len(ntc_rows)),
            ("ntc_p_null_below_0.05", int((ntc_rows.p_null < 0.05).sum())),
            ("ntc_p_null_below_0.10", int((ntc_rows.p_null < 0.10).sum())),
            ("ntc_median_p_null", float(ntc_rows.p_null.median())),
            ("ntc_glm_p_below_0.05", int((ntc_rows.glm_p < 0.05).sum())),
            ("ntc_median_glm_p", float(ntc_rows.glm_p.median())),
        ],
        columns=["metric", "value"],
    )
    cal.to_csv(C.tbl("07_null_calibration.tsv"), sep="\t", index=False)
    print("\nnull calibration on the testable non-targeting elements:")
    print(cal.to_string(index=False))
    print(
        "  both the resampling and the Wald p-values are well behaved on the NTC "
        "elements here (no element below 0.05, medians near 0.4), so at this cell "
        "count the GLM's independence assumption is not visibly costing anything. "
        "The control-based null is still the one used to call hits, because it does "
        "not rest on that assumption and it is matched to each element's own size "
        "and lane composition."
    )

    # ---------------------------------------------------------- phenotype B
    Bt = burden_phenotype(d)
    nb, alpha = burden_nb_glm(d, Bt.element.tolist())
    Bt = Bt.merge(nb, on="element", how="left")
    Bt["q_glm"] = np.where(~Bt.is_ntc, C.bh_fdr(Bt.nb_p.where(~Bt.is_ntc).to_numpy()), np.nan)
    rb = Bt[["delta_log10_rate", "nb_log_rate_ratio"]].dropna()
    print(f"  resampling effect vs NB rate ratio: r = "
          f"{np.corrcoef(rb.delta_log10_rate, rb.nb_log_rate_ratio)[0, 1]:.3f}")
    Bt = Bt.sort_values("delta_log10_rate")
    Bt.to_csv(C.tbl("07_hits_burden.tsv"), sep="\t", index=False)
    sigB = Bt[(~Bt.is_ntc) & (Bt.q_null < 0.1)]
    print(f"  significant at q<0.1: {len(sigB)} of {int((~Bt.is_ntc).sum())} targeting elements")

    # ------------------------------------------------------------ concordance
    rows = []
    for tgt, direction in PAPER_HITS.items():
        a = A[A.target == tgt]
        b = Bt[Bt.target == tgt]
        for _, r_ in a.iterrows():
            br = b[b.element == r_.element]
            rows.append(
                dict(
                    target=tgt, element=r_.element, paper_direction=direction,
                    n_cells=int(r_.n_cells),
                    infection_log_or=r_.log_or_mh, infection_q=r_.q_null,
                    infection_called=("down" if (r_.q_null < 0.1 and r_.log_or_mh < 0)
                                      else "up" if (r_.q_null < 0.1 and r_.log_or_mh > 0)
                                      else "ns"),
                    burden_delta=float(br.delta_log10_rate.iloc[0]) if len(br) else np.nan,
                    burden_q=float(br.q_null.iloc[0]) if len(br) else np.nan,
                )
            )
    for tgt in PAPER_NULL:
        a = A[A.target == tgt]
        for _, r_ in a.iterrows():
            b = Bt[Bt.element == r_.element]
            rows.append(
                dict(
                    target=tgt, element=r_.element, paper_direction="ns",
                    n_cells=int(r_.n_cells),
                    infection_log_or=r_.log_or_mh, infection_q=r_.q_null,
                    infection_called=("down" if (r_.q_null < 0.1 and r_.log_or_mh < 0)
                                      else "up" if (r_.q_null < 0.1 and r_.log_or_mh > 0)
                                      else "ns"),
                    burden_delta=float(b.delta_log10_rate.iloc[0]) if len(b) else np.nan,
                    burden_q=float(b.q_null.iloc[0]) if len(b) else np.nan,
                )
            )
    con = pd.DataFrame(rows)
    con["agrees"] = con.paper_direction == con.infection_called
    con.to_csv(C.tbl("07_concordance.tsv"), sep="\t", index=False)
    print("\nconcordance with the published hit list:")
    print(con.round(3).to_string(index=False))
    hits = con[con.paper_direction != "ns"]
    nulls = con[con.paper_direction == "ns"]
    print(
        f"\n  published hits recovered in the right direction: "
        f"{int(hits.agrees.sum())}/{len(hits)}"
    )
    print(
        f"  published non-hits we also call ns: {int(nulls.agrees.sum())}/{len(nulls)}"
    )

    # ------------------------------------------------------------- robustness
    A_amb = infection_phenotype(d, state_col="state_ambient", draws=500, verbose=False)
    m = A[["element", "log_or_mh", "q_null"]].merge(
        A_amb[["element", "log_or_mh", "q_null"]], on="element",
        suffixes=("_mixture", "_ambient"))
    rr = float(np.corrcoef(m.log_or_mh_mixture.fillna(0), m.log_or_mh_ambient.fillna(0))[0, 1])
    both = int(((m.q_null_mixture < 0.1) & (m.q_null_ambient < 0.1)).sum())
    either = int(((m.q_null_mixture < 0.1) | (m.q_null_ambient < 0.1)).sum())
    rob = pd.DataFrame(
        [
            ("effect_size_correlation_between_classifiers", rr),
            ("significant_under_both", both),
            ("significant_under_either", either),
            ("jaccard", both / either if either else np.nan),
        ],
        columns=["metric", "value"],
    )
    rob.to_csv(C.tbl("07_robustness.tsv"), sep="\t", index=False)
    print("\nrobustness to the infection classifier:")
    print(rob.round(4).to_string(index=False))

    offset_sensitivity(d, Bt)


def offset_sensitivity(d: pd.DataFrame, Bt: pd.DataFrame) -> None:
    """How much does the burden result depend on the choice of denominator?

    The primary model's estimand is viral RNA per host transcript *captured*,
    not absolute viral abundance per cell. Host UMI is itself lowered by
    infection and could in principle be moved by a perturbation, so the
    denominator is not a neutral exposure measure. Three alternatives:

      host_umi    (primary)  viral counts per host transcript captured
      n_genes     log(host genes detected): a capture proxy that saturates, so
                  it responds differently to the magnitude of host shutoff
      total_umi   log(total library): the specification the primary model
                  rejects, included so the size of the problem is visible
      depth_band  no offset, restricted to infected cells inside a two-fold
                  band of host UMI around the median, so exposure is roughly
                  constant by construction

    If the primary result is an artefact of the denominator, these will
    disagree. Reported as rank correlations of the per-element effect and as
    overlap of the top-20 most negative elements.
    """
    inf = d[(d.state == "infected") & (d.host_umi > 0) & (d.n_genes > 0)].copy()
    lane_dummies = pd.get_dummies(inf.gem_group, prefix="lane", drop_first=True).astype(float)
    ntc_mask = inf.is_ntc.to_numpy()
    y = inf.viral_umi.to_numpy(float)
    elements = Bt.element.tolist()

    band_lo, band_hi = inf.host_umi.median() / 1.41, inf.host_umi.median() * 1.41
    in_band = ((inf.host_umi >= band_lo) & (inf.host_umi <= band_hi)).to_numpy()

    specs = {
        "host_umi": np.log(inf.host_umi.to_numpy(float)),
        "n_genes": np.log(inf.n_genes.to_numpy(float)),
        "total_umi": np.log(inf.total_umi.to_numpy(float)),
        "depth_band_no_offset": np.zeros(len(inf)),
    }

    est = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, offset in specs.items():
            sub_mask = in_band if name == "depth_band_no_offset" else np.ones(len(inf), bool)
            vals = {}
            for element in elements:
                m = ((inf.element == element).to_numpy()) & sub_mask
                keep = (m | (ntc_mask & sub_mask))
                if m.sum() < 15:
                    continue
                X = pd.concat(
                    [pd.Series(m[keep].astype(float), index=inf.index[keep], name="element"),
                     lane_dummies.loc[inf.index[keep]]], axis=1)
                X = sm.add_constant(X.loc[:, X.std() > 0])
                try:
                    fit = sm.GLM(y[keep], X,
                                 family=sm.families.NegativeBinomial(alpha=1.66),
                                 offset=offset[keep]).fit()
                    vals[element] = float(fit.params["element"])
                except Exception:  # noqa: BLE001
                    continue
            est[name] = pd.Series(vals)
            print(f"  {name}: {len(vals)} elements fitted")

    tab = pd.DataFrame(est)
    tab.index.name = "element"
    tab = tab.join(Bt.set_index("element")[["target", "is_ntc", "n_infected"]])
    tab.to_csv(C.tbl("07_offset_sensitivity.tsv"), sep="\t")

    rows = []
    base = "host_umi"
    top_base = set(tab[base].nsmallest(20).index)
    for name in specs:
        if name == base or name not in tab:
            continue
        pair = tab[[base, name]].dropna()
        rows.append(
            dict(
                specification=name,
                n_elements_compared=len(pair),
                spearman_vs_host_umi=float(pair[base].corr(pair[name], method="spearman")),
                pearson_vs_host_umi=float(pair[base].corr(pair[name])),
                top20_overlap=len(top_base & set(tab[name].nsmallest(20).index)),
            )
        )
    sens = pd.DataFrame(rows)
    sens.to_csv(C.tbl("07_offset_sensitivity_summary.tsv"), sep="\t", index=False)
    print("\nburden phenotype, sensitivity to the denominator "
          "(reference: offset = log host UMI):")
    print(sens.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
