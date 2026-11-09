"""Shared helpers: paths, loading, gene sets, and the pseudobulk/offset conventions."""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# --------------------------------------------------------------------- paths
REPO = Path(__file__).resolve().parent.parent
TABLES = REPO / "results" / "tables"
FIGURES = REPO / "results" / "figures"


def h5ad_path() -> Path:
    """Where 01 wrote the object. Override with PERTURBSEQ_H5AD."""
    env = os.environ.get("PERTURBSEQ_H5AD")
    if env:
        return Path(env)
    return REPO / "data" / "perturbseq.h5ad"


def load(assigned_only: bool = True) -> ad.AnnData:
    a = ad.read_h5ad(h5ad_path())
    if assigned_only:
        a = a[a.obs.cell_set == "assigned"].copy()
    return a


def tbl(name: str) -> Path:
    TABLES.mkdir(parents=True, exist_ok=True)
    return TABLES / name


def read_tbl(name: str) -> pd.DataFrame:
    p = TABLES / name
    return pd.read_csv(p, sep="\t")


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


# ----------------------------------------------------------------- utilities
def host_cpm(a: ad.AnnData) -> np.ndarray:
    """Host-only counts per million: the denominator excludes viral and lenti
    reads, so that viral takeover does not rescale host expression."""
    host = (a.var.kind == "host").to_numpy()
    X = a[:, host].X
    denom = np.asarray(X.sum(axis=1)).ravel()
    denom[denom == 0] = 1
    return X.multiply(1e6 / denom[:, None]).tocsr()


def score_gene_set(a: ad.AnnData, genes: list[str], layer_cpm=None) -> np.ndarray:
    """Mean log1p(host CPM) over the gene set minus the mean over a size- and
    expression-matched background, i.e. the scanpy score_genes construction,
    computed on host-only CPM."""
    import scanpy as sc

    b = a.copy()
    b.X = layer_cpm if layer_cpm is not None else host_cpm_full(a)
    sc.pp.log1p(b)
    present = [g for g in genes if g in b.var_names]
    sc.tl.score_genes(b, present, score_name="s", random_state=0)
    return b.obs["s"].to_numpy(), present


def host_cpm_full(a: ad.AnnData):
    """Host-normalised CPM but keeping all columns, so gene indexing is stable."""
    host = (a.var.kind == "host").to_numpy()
    denom = np.asarray(a[:, host].X.sum(axis=1)).ravel()
    denom[denom == 0] = 1
    return a.X.multiply(1e6 / denom[:, None]).tocsr()


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
