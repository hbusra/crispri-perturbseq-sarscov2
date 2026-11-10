"""Loading the processed object, and host-only normalisation."""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import numpy as np

REPO = Path(__file__).resolve().parents[2]


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
