#!/usr/bin/env python
"""
Correctness tests for the estimators the hit calls rest on.

Run with:  python -m pytest tests -q      (or: python tests/test_estimators.py)

Each test checks one of the claims the README makes, against an independent
implementation or against a simulation with a known answer:

  * the Mantel-Haenszel log odds ratio and its standard error agree with
    statsmodels' StratifiedTable on the same data;
  * the MH estimator recovers a known odds ratio from simulated data, and is
    not fooled by a stratum effect that an unstratified estimator falls for;
  * the negative-binomial GLM with a host-UMI offset recovers a known rate
    ratio;
  * the Benjamini-Hochberg implementation matches statsmodels' multipletests;
  * the resampling null is calibrated: feeding it null data gives uniform
    p-values and the nominal false-positive rate;
  * the two-sided empirical p carries the plus-one correction, so it can never
    be zero.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "notebook"))

from methods import bh_fdr  # noqa: E402
from methods import (  # noqa: E402
    lane_weighted_delta,
    mh_log_or,
    mh_null_draws,
    resample_delta_null,
    two_sided_p,
)


# ------------------------------------------------------- Mantel-Haenszel
def test_mh_matches_statsmodels():
    """Our MH point estimate and SE must equal statsmodels' StratifiedTable."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        k = rng.integers(3, 9)
        a = rng.integers(5, 80, k)
        b = rng.integers(5, 80, k)
        c = rng.integers(5, 80, k)
        d = rng.integers(5, 80, k)

        ours, our_se = mh_log_or(a, b, c, d)

        # statsmodels wants 2x2xK with rows = group, cols = outcome
        tables = np.array([[[a[i], b[i]], [c[i], d[i]]] for i in range(k)])
        st = StratifiedTable(np.moveaxis(tables, 0, -1))

        assert ours == pytest.approx(st.logodds_pooled, rel=1e-10)
        assert our_se == pytest.approx(st.logodds_pooled_se, rel=1e-8)


def test_mh_recovers_known_odds_ratio():
    """Simulate a known OR with strong stratum effects and recover it."""
    rng = np.random.default_rng(1)
    true_log_or = np.log(0.4)
    # eight lanes whose baseline infection rate varies a lot, as ours does
    base = np.array([0.18, 0.20, 0.19, 0.21, 0.34, 0.33, 0.35, 0.34])
    n_el, n_ntc = 400, 2000

    a, b, c, d = [], [], [], []
    for p0 in base:
        odds0 = p0 / (1 - p0)
        p1 = (odds0 * np.exp(true_log_or)) / (1 + odds0 * np.exp(true_log_or))
        ai = rng.binomial(n_el, p1)
        ci = rng.binomial(n_ntc, p0)
        a.append(ai); b.append(n_el - ai)
        c.append(ci); d.append(n_ntc - ci)

    est, se = mh_log_or(a, b, c, d)
    assert abs(est - true_log_or) < 3 * se, f"{est:.3f} vs {true_log_or:.3f} (se {se:.3f})"


def test_stratification_actually_matters():
    """A confounded design: no true effect, but the element sits in the
    high-infection lanes. Stratified must say null; pooled must not."""
    rng = np.random.default_rng(2)
    base = np.array([0.18, 0.18, 0.18, 0.18, 0.40, 0.40, 0.40, 0.40])
    # element cells only in lanes 5-8, NTC cells spread evenly
    n_el = np.array([0, 0, 0, 0, 300, 300, 300, 300])
    n_ntc = np.full(8, 1000)

    a = rng.binomial(n_el, base)
    c = rng.binomial(n_ntc, base)
    b, d = n_el - a, n_ntc - c

    strat, se = mh_log_or(a, b, c, d)
    pooled = np.log((a.sum() / b.sum()) / (c.sum() / d.sum()))

    assert abs(strat) < 3 * se, "stratified estimate should be null"
    # the pooled estimate sees a ~1.6x spurious odds ratio purely from where
    # the element's cells sit
    assert pooled > 0.4, f"pooled estimate should be confounded, got {pooled:.3f}"
    assert pooled > 4 * abs(strat), "pooling should be much worse than stratifying"


# -------------------------------------------------- negative binomial GLM
def test_nb_glm_recovers_rate_ratio():
    """Viral counts ~ NB with log(host_umi) offset: recover a known ratio."""
    rng = np.random.default_rng(3)
    n = 4000
    true_log_rr = np.log(0.5)
    alpha = 1.5

    host_umi = rng.lognormal(mean=np.log(8000), sigma=0.6, size=n)
    is_element = rng.binomial(1, 0.25, n)
    lane = rng.integers(0, 8, n)
    lane_effect = np.linspace(-0.3, 0.3, 8)[lane]

    log_rate = np.log(0.02) + true_log_rr * is_element + lane_effect
    mu = np.exp(log_rate + np.log(host_umi))
    # NB2 as a gamma-Poisson mixture with dispersion alpha
    y = rng.poisson(rng.gamma(shape=1 / alpha, scale=mu * alpha))

    X = pd.DataFrame({"element": is_element.astype(float)})
    for l in range(1, 8):
        X[f"lane_{l}"] = (lane == l).astype(float)
    X = sm.add_constant(X)

    fit = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha),
                 offset=np.log(host_umi)).fit()
    est = float(fit.params["element"])
    se = float(fit.bse["element"])
    assert abs(est - true_log_rr) < 3 * se, f"{est:.3f} vs {true_log_rr:.3f}"


def test_offset_choice_changes_the_answer():
    """The host-UMI offset is not cosmetic: using a denominator that contains
    the response biases the estimate toward null. This is the behaviour the
    repository's choice of offset is meant to avoid."""
    rng = np.random.default_rng(4)
    n = 4000
    true_log_rr = np.log(0.4)
    host_umi = rng.lognormal(np.log(8000), 0.5, n)
    is_element = rng.binomial(1, 0.3, n)
    mu = np.exp(np.log(0.05) + true_log_rr * is_element + np.log(host_umi))
    y = rng.poisson(mu)
    total_umi = host_umi + y            # what a library actually contains

    X = sm.add_constant(pd.DataFrame({"element": is_element.astype(float)}))
    good = sm.GLM(y, X, family=sm.families.Poisson(),
                  offset=np.log(host_umi)).fit().params["element"]
    bad = sm.GLM(y, X, family=sm.families.Poisson(),
                 offset=np.log(total_umi)).fit().params["element"]

    assert abs(good - true_log_rr) < 0.05
    assert abs(bad - true_log_rr) > abs(good - true_log_rr), (
        "the total-UMI offset should be the worse estimator"
    )


# ------------------------------------------------------------------ FDR
def test_bh_matches_statsmodels():
    rng = np.random.default_rng(5)
    for _ in range(10):
        p = np.concatenate([rng.uniform(0, 1, 150), rng.uniform(0, 1e-3, 15)])
        ours = bh_fdr(p)
        theirs = multipletests(p, method="fdr_bh")[1]
        assert np.allclose(ours, theirs)


def test_bh_handles_nan():
    p = np.array([0.01, np.nan, 0.5, 0.001])
    q = bh_fdr(p)
    assert np.isnan(q[1])
    assert np.allclose(q[[0, 2, 3]], multipletests(p[[0, 2, 3]], method="fdr_bh")[1])


# ------------------------------------------------------- the empirical null
def test_two_sided_p_has_plus_one_and_is_never_zero():
    null = np.linspace(-1, 1, 1000)
    assert two_sided_p(99.0, null) == pytest.approx(1 / 1001)
    assert two_sided_p(0.0, null) == pytest.approx(1.0)
    assert two_sided_p(5.0, null) > 0


def test_mh_null_is_calibrated():
    """Under the null, the resampling p-values must be ~uniform, so the
    false-positive rate at 0.05 is ~5%."""
    rng = np.random.default_rng(6)
    lanes = list(range(8))
    p_l = np.array([0.18, 0.20, 0.19, 0.21, 0.34, 0.33, 0.35, 0.34])
    n_ntc = np.full(8, 1200)
    c_l = rng.binomial(n_ntc, p_l)
    d_l = n_ntc - c_l

    ps = []
    for _ in range(300):
        n_l = rng.integers(8, 30, 8)           # a ~150-cell pseudo element
        a_l = rng.binomial(n_l, p_l)           # drawn under the null
        b_l = n_l - a_l
        obs, _se = mh_log_or(a_l, b_l, c_l, d_l)
        null = mh_null_draws(n_l, p_l, c_l, d_l, 400, rng)
        ps.append(two_sided_p(obs, null))
    ps = np.array([p for p in ps if np.isfinite(p)])

    fpr = float((ps < 0.05).mean())
    assert 0.01 < fpr < 0.11, f"false positive rate {fpr:.3f} is not near 0.05"
    assert 0.40 < ps.mean() < 0.60, f"mean p {ps.mean():.3f} is not near 0.5"


def test_mh_null_detects_a_real_effect():
    rng = np.random.default_rng(7)
    p_l = np.full(8, 0.30)
    c_l = rng.binomial(np.full(8, 1200), p_l)
    d_l = 1200 - c_l
    n_l = np.full(8, 25)                        # 200 cells
    p_el = np.full(8, 0.10)                     # a real, large effect
    a_l = rng.binomial(n_l, p_el)
    b_l = n_l - a_l
    obs, _ = mh_log_or(a_l, b_l, c_l, d_l)
    p = two_sided_p(obs, mh_null_draws(n_l, p_l, c_l, d_l, 2000, rng))
    assert obs < 0 and p < 0.01


def test_resample_delta_null_is_calibrated():
    """The bootstrap null used for the burden and ISG phenotypes."""
    rng = np.random.default_rng(8)
    lanes = [str(i) for i in range(1, 9)]
    pool = {l: rng.normal(float(i) * 0.1, 1.0, 400) for i, l in enumerate(lanes)}
    base = {l: float(np.mean(v)) for l, v in pool.items()}

    ps = []
    for _ in range(300):
        # a pseudo element drawn from the same pools: null is true
        frames = []
        for l in lanes:
            k = rng.integers(3, 12)
            frames.append(pd.DataFrame({
                "gem_group": l,
                "x": pool[l][rng.integers(0, len(pool[l]), k)],
            }))
        sub = pd.concat(frames)
        obs, n_l = lane_weighted_delta(sub, lanes, "x", base)
        null = resample_delta_null(pool, base, lanes, n_l, 400, rng)
        ps.append(two_sided_p(obs, null))

    ps = np.array(ps)
    fpr = float((ps < 0.05).mean())
    assert 0.01 < fpr < 0.11, f"false positive rate {fpr:.3f}"


def test_lane_weighted_delta_removes_a_lane_effect():
    """If an element sits only in lanes with a high baseline, the lane-weighted
    statistic must still be ~0 when there is no element effect."""
    rng = np.random.default_rng(9)
    lanes = [str(i) for i in range(1, 9)]
    pool = {l: rng.normal(0 if int(l) <= 4 else 3.0, 0.5, 500) for l in lanes}
    base = {l: float(np.mean(v)) for l, v in pool.items()}

    sub = pd.concat([
        pd.DataFrame({"gem_group": l, "x": pool[l][rng.integers(0, 500, 40)]})
        for l in lanes[4:]                      # only the high-baseline lanes
    ])
    obs, _ = lane_weighted_delta(sub, lanes, "x", base)
    naive = sub.x.mean() - np.mean(np.concatenate(list(pool.values())))
    assert abs(obs) < 0.15, f"lane-weighted statistic should be ~0, got {obs:.3f}"
    assert naive > 1.0, "the naive difference should be badly confounded"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
