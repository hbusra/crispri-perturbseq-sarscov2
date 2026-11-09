"""Estimators and null-construction helpers shared by 07 and 08.

Kept in one place so the hit-calling statistic used to call hits and the one
used to study the screen's resolution are provably the same function.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def mh_log_or(a, b, c, d) -> tuple[float, float]:
    """Mantel-Haenszel log odds ratio across strata, with the
    Robins-Breslow-Greenland standard error.

    a, b : per-stratum infected / bystander counts for the element
    c, d : per-stratum infected / bystander counts for the pooled NTC cells

    Stratifying on the 10x channel keeps the lane effect out of the element
    effect without spending a degree of freedom per lane on every element.
    """
    a, b, c, d = (np.asarray(x, float) for x in (a, b, c, d))
    n = a + b + c + d
    ok = n > 0
    a, b, c, d, n = a[ok], b[ok], c[ok], d[ok], n[ok]
    R = np.sum(a * d / n)
    S = np.sum(b * c / n)
    if R <= 0 or S <= 0:
        return np.nan, np.nan
    log_or = float(np.log(R / S))
    P = (a + d) / n
    Q = (b + c) / n
    var = (
        np.sum(P * a * d / n) / (2 * R**2)
        + np.sum(P * b * c / n + Q * a * d / n) / (2 * R * S)
        + np.sum(Q * b * c / n) / (2 * S**2)
    )
    return log_or, float(np.sqrt(var))


def lane_counts(df: pd.DataFrame, lanes, state_col: str = "state"):
    """Infected / bystander counts per lane, in the order of `lanes`."""
    inf = df.loc[df[state_col] == "infected", "gem_group"].value_counts()
    byst = df.loc[df[state_col] == "uninfected", "gem_group"].value_counts()
    return (
        inf.reindex(lanes, fill_value=0).to_numpy(),
        byst.reindex(lanes, fill_value=0).to_numpy(),
    )


def mh_null_draws(n_l, p_l, c_l, d_l, draws: int, rng) -> np.ndarray:
    """Lane- and size-matched null for the MH log odds ratio.

    For each draw, the element's cells in lane l are replaced by a binomial
    draw at the pooled NTC infection rate of that lane. This is the exchangeable
    null: same number of cells, same lanes, no perturbation effect.
    """
    n_l = np.asarray(n_l)
    a_null = rng.binomial(np.repeat(n_l[None, :], draws, axis=0), np.asarray(p_l)[None, :])
    b_null = n_l[None, :] - a_null
    denom = n_l + c_l + d_l
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(
            np.sum(a_null * d_l / denom, axis=1) / np.sum(b_null * c_l / denom, axis=1)
        )


def two_sided_p(obs: float, null: np.ndarray) -> float:
    """Empirical two-sided p, floored at 1/(B+1)."""
    null = np.asarray(null)
    null = null[np.isfinite(null)]
    if len(null) == 0 or not np.isfinite(obs):
        return np.nan
    return (1 + int(np.sum(np.abs(null) >= abs(obs)))) / (len(null) + 1)


def lane_weighted_delta(sub: pd.DataFrame, lanes, col: str, base: dict):
    """Lane-weighted mean difference of `col` against a per-lane NTC baseline.

    Returns (statistic, per-lane cell counts). Lanes where the element has no
    cells, or the NTC baseline is undefined, contribute nothing.
    """
    n_l, num, w = [], 0.0, 0.0
    for l in lanes:
        s = sub.loc[sub.gem_group == l, col].to_numpy()
        if len(s) == 0 or not np.isfinite(base.get(l, np.nan)):
            n_l.append(0)
            continue
        n_l.append(len(s))
        num += len(s) * (s.mean() - base[l])
        w += len(s)
    return (num / w if w else np.nan), n_l


def resample_delta_null(pool: dict, base: dict, lanes, n_l, draws: int, rng) -> np.ndarray:
    """Null for `lane_weighted_delta`: resample the NTC pool within each lane,
    matching the element's per-lane cell counts."""
    w = sum(n_l)
    out = np.empty(draws)
    for i in range(draws):
        num = 0.0
        for l, k in zip(lanes, n_l):
            if k == 0:
                continue
            p = pool[l]
            num += k * (p[rng.integers(0, len(p), k)].mean() - base[l])
        out[i] = num / w
    return out
