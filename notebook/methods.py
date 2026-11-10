"""Shared methods for the notebooks: paths, loading, estimators, gene sets, style.

Only things that are reused by more than one notebook live here, so that the
notebooks themselves stay readable. The estimators are the part worth keeping
in a module: tests/test_estimators.py checks them against statsmodels and
against simulated data where the truth is known.
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"
FIGURES = REPO / "results" / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)


def table(name: str) -> pd.DataFrame:
    """Read one of the tables written by an earlier notebook."""
    return pd.read_csv(TABLES / name, sep="\t")

# ------------------------------------------------------------------ loading
def h5ad_path() -> Path:
    """Where 01_build_anndata.py wrote the object. Override with PERTURBSEQ_H5AD."""
    env = os.environ.get("PERTURBSEQ_H5AD")
    return Path(env) if env else REPO / "data" / "perturbseq.h5ad"


def load_cells(assigned_only: bool = True) -> ad.AnnData:
    a = ad.read_h5ad(h5ad_path())
    if assigned_only:
        a = a[a.obs.cell_set == "assigned"].copy()
    return a


def host_cpm(a: ad.AnnData):
    """Counts per million with a host-only denominator, keeping all columns.

    The denominator excludes viral and lentiviral reads, so that viral takeover
    of a library does not rescale host expression.
    """
    host = (a.var.kind == "host").to_numpy()
    denom = np.asarray(a[:, host].X.sum(axis=1)).ravel()
    denom[denom == 0] = 1
    return a.X.multiply(1e6 / denom[:, None]).tocsr()

# --------------------------------------------------------------- estimators
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


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q

# ----------------------------------------------------------------- gene sets
# Interferon-stimulated genes. Core type-I ISGs that are unambiguous markers of
# an autocrine/paracrine IFN response in epithelial cells. Compiled from the
# MSigDB Hallmark IFN-alpha response leading edge, restricted to genes detected
# in this dataset at 01. Script 06 reports the score with and without the
# genes that are themselves CRISPRi targets in this library, since knocking a
# gene down mechanically lowers any score that contains it.
ISG_GENES = [
    "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "MX1", "MX2",
    "OAS1", "OAS2", "OAS3", "OASL", "IFI6", "IFI27", "IFI35", "IFI44",
    "IFI44L", "IFITM1", "IFITM2", "IFITM3", "STAT1", "STAT2", "IRF7", "IRF9",
    "BST2", "RSAD2", "CMPK2", "XAF1", "HERC5", "HERC6", "EIF2AK2", "DDX58",
    "IFIH1", "PARP9", "DTX3L", "SAMD9", "SAMD9L", "EPSTI1", "LY6E", "LAP3",
    "PLSCR1", "TRIM22", "SP100", "SP110", "B2M", "PSMB8", "PSMB9", "TAP1",
    "NMI", "CXCL10", "CXCL11", "IFI16", "ZBP1", "APOL6", "ADAR",
]

# Type I/III interferons and the chemokines that mark the rare "IFN-producing"
# subset (bystander cluster L / infected cluster M in the paper).
IFN_PRODUCING = [
    "IFNB1", "IFNL1", "IFNL2", "IFNL3", "IFNA1", "IL6", "CXCL8", "CCL5",
    "TNF", "IL1B",
]

# NF-kB target genes reported as induced in infected cells (paper Fig 2E).
NFKB_GENES = [
    "NFKBIA", "NFKBIE", "NFKBIZ", "NFKB1", "NFKB2", "REL", "RELA", "RELB",
    "EGR1", "TNFAIP3", "BIRC3", "CXCL1", "CXCL2", "CXCL3", "SOD2", "ICAM1",
]

VIRAL_GENES = [
    "SCV_ORF1ab", "SCV_S", "SCV_ORF3a", "SCV_E", "SCV_M",
    "SCV_ORF6", "SCV_ORF7ab", "SCV_ORF8", "SCV_N",
]

# ---------------------------------------------------------- palette and style
import matplotlib.pyplot as plt

# ----------------------------------------------------------------- palette
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED = "#e34948"
NEUTRAL = "#b4b3ad"
GRID = "#e6e5e1"
INK = "#0b0b0b"
INK2 = "#52514e"

STATE_COLOR = {"uninfected": BLUE, "infected": ORANGE, "borderline": NEUTRAL}
STATE_LABEL = {"uninfected": "bystander", "infected": "infected", "borderline": "borderline"}

plt.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.5,
    "axes.titleweight": "semibold",
    "axes.edgecolor": "#c9c8c3",
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK2,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.frameon": False,
    "legend.fontsize": 7.5,
    "lines.linewidth": 2.0,
})


def style(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
