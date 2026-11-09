#!/usr/bin/env python
"""
01 - Build the analysis AnnData from the GSE208240 CellRanger matrix.

Inputs  (from scripts/00_download.sh)
    <RAW>/matrix_dir/{barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz,cell_identities.csv}

Output
    <OUT>/perturbseq.h5ad
    <OUT>/X_kept.npz                 cached subset matrix, so reruns are fast
    results/tables/01_build_summary.tsv

Design notes
------------
* The published matrix is 90,380 barcodes x 33,549 features: 33,538 human genes
  + 9 SARS-CoV-2 genes (SCV_*) + 2 lentiviral constructs (dCas9-BFP-KRAB,
  Puro-t2a-GFP).
* Two populations of cells are kept:
    1. `assigned`     - match_type == exact_match and good_coverage == True.
                        This is the 27,882-cell analysis set of the paper.
    2. `wt_candidate` - barcodes with no guide detected. The authors spiked in
                        ~1% wild-type (no dCas9, no sgRNA) cells and used them
                        to measure ambient viral RNA for infection calling.
                        The lentivirus-negative and depth filters are applied
                        in script 03, so the reference population is derived
                        in the open rather than hard-coded here.
  Multiplets, scrambled pairs and likely-but-not-exact matches are dropped,
  which is what the authors' analysis set does.
* Guides are never collapsed to genes at build time. `element` (e.g.
  CRISPRi_022_BRD2) is the experimental unit carried forward. Elements 216-239
  are 12 genes with two *independent* single-guide vectors each, so keeping
  `element` distinct lets script 07 use them as an internal reproducibility
  check.
* `gem_group` (1-8) is the 10x channel: a real blocking variable, carried into
  every model and every control-based null downstream.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

CHUNK = 20_000_000  # mtx lines per read_csv chunk

EXPECTED_CELLS_ASSIGNED = 27_882
EXPECTED_ELEMENTS = 237
EXPECTED_FEATURES = 33_549


def parse_element(name: str) -> tuple[str, bool, int]:
    """CRISPRi_022_BRD2 -> ('BRD2', False, 1)
    CRISPRi_182_IRF9_IRF3 -> ('IRF9+IRF3', False, 2)
    CRISPRi_212_non_targeting_pair_19 -> ('NTC', True, 0)"""
    if name in ("n/a", "", None) or pd.isna(name):
        return ("unassigned", False, 0)
    parts = name.split("_", 2)
    rest = parts[2] if len(parts) == 3 else name
    if rest.startswith("non_targeting"):
        return ("NTC", True, 0)
    genes = rest.split("_")
    return ("+".join(genes), False, len(genes))


def read_matrix_subset(mtx_gz: Path, remap: np.ndarray, n_cells: int, n_feat: int,
                       n_bc: int) -> sp.csr_matrix:
    """Stream the .mtx and keep only the columns mapped to >= 0 in `remap`."""
    print("streaming matrix.mtx.gz ...", flush=True)
    with gzip.open(mtx_gz, "rt") as fh:
        assert fh.readline().startswith("%%MatrixMarket")
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_feat_mtx, n_bc_mtx, n_nz = (int(x) for x in line.split())
        assert (n_feat_mtx, n_bc_mtx) == (n_feat, n_bc), "matrix header does not match the tsv files"
        print(f"  header: {n_feat_mtx} x {n_bc_mtx}, {n_nz} non-zeros")

        rows_l, cols_l, vals_l = [], [], []
        seen = 0
        reader = pd.read_csv(
            fh, sep=" ", header=None, names=["g", "c", "v"],
            dtype={"g": np.int32, "c": np.int32, "v": np.int32},
            chunksize=CHUNK, engine="c",
        )
        for i, chunk in enumerate(reader):
            seen += len(chunk)
            newc = remap[chunk.c.to_numpy()]
            m = newc >= 0
            if m.any():
                rows_l.append(newc[m].astype(np.int32))
                cols_l.append(chunk.g.to_numpy()[m] - 1)
                vals_l.append(chunk.v.to_numpy()[m])
            print(f"  chunk {i}: {seen}/{n_nz} lines, kept {int(m.sum())}", flush=True)
        assert seen == n_nz, f"read {seen} entries, header said {n_nz}"

    rows = np.concatenate(rows_l)
    cols = np.concatenate(cols_l)
    vals = np.concatenate(vals_l)
    print(f"  kept {len(vals)} non-zeros ({100 * len(vals) / n_nz:.1f}% of the file)")
    X = sp.coo_matrix(
        (vals.astype(np.float32), (rows, cols)), shape=(n_cells, n_feat)
    ).tocsr()
    return X


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True, help="directory containing matrix_dir/")
    ap.add_argument("--out", required=True, help="output directory for the .h5ad")
    ap.add_argument("--tables", default="results/tables")
    args = ap.parse_args()

    raw = Path(args.raw)
    mdir = raw / "matrix_dir"
    if not mdir.exists():
        cand = list(raw.rglob("CRISPRi_perturbseq_sarscov2_filtered"))
        if not cand:
            sys.exit(f"could not find the matrix directory under {raw}")
        mdir = cand[0]

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tables = Path(args.tables)
    tables.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- features
    feat = pd.read_csv(
        mdir / "features.tsv.gz", sep="\t", header=None,
        names=["gene_id", "gene_name", "feature_type"],
    )
    assert len(feat) == EXPECTED_FEATURES, f"expected {EXPECTED_FEATURES} features, got {len(feat)}"
    feat["kind"] = np.where(
        feat.gene_id.str.startswith("SCV_"), "viral",
        np.where(feat.gene_id.str.startswith("lenti_"), "lenti", "host"),
    )
    print("features by kind:", feat.kind.value_counts().to_dict())

    # ---------------------------------------------------------------- barcodes
    barcodes = pd.read_csv(mdir / "barcodes.tsv.gz", header=None)[0].to_numpy()
    n_bc = len(barcodes)
    print(f"barcodes: {n_bc}")

    # ------------------------------------------------------- guide assignments
    ci = pd.read_csv(mdir / "cell_identities.csv").set_index("cell_barcode").reindex(barcodes)

    assigned = (ci.match_type == "exact_match") & (ci.good_coverage == True)  # noqa: E712
    wt_cand = (ci.match_type == "no_match") & (ci.number_of_guides.fillna(0) == 0)
    keep = (assigned | wt_cand).to_numpy()
    n_cells = int(keep.sum())

    print(f"assigned cells      : {int(assigned.sum())}")
    print(f"wt candidate cells  : {int(wt_cand.sum())}")
    print(f"kept for the h5ad   : {n_cells}")
    assert int(assigned.sum()) == EXPECTED_CELLS_ASSIGNED, (
        f"expected {EXPECTED_CELLS_ASSIGNED} exact-match cells, got {int(assigned.sum())}"
    )

    # old column index (1-based in the mtx) -> new row index, or -1 to drop
    remap = np.full(n_bc + 1, -1, dtype=np.int64)
    remap[np.where(keep)[0] + 1] = np.arange(n_cells)

    # -------------------------------------------------------------- the matrix
    cache = out / "X_kept.npz"
    if cache.exists():
        print(f"loading cached matrix from {cache}")
        X = sp.load_npz(cache).tocsr()
        assert X.shape == (n_cells, len(feat)), "cached matrix shape mismatch; delete the cache"
    else:
        X = read_matrix_subset(mdir / "matrix.mtx.gz", remap, n_cells, len(feat), n_bc)
        sp.save_npz(cache, X)
        print(f"cached subset matrix -> {cache}")

    # ----------------------------------------------------------------- obs/var
    obs = ci.loc[keep].copy()
    obs.index = pd.Index(barcodes[keep], name="cell_barcode")
    obs["cell_set"] = np.where(assigned.to_numpy()[keep], "assigned", "wt_candidate")
    obs["element"] = obs["matched_library_element"].where(obs.cell_set == "assigned", "unassigned")

    parsed = [parse_element(e) for e in obs["element"]]
    obs["target"] = [p[0] for p in parsed]
    obs["is_ntc"] = [p[1] for p in parsed]
    obs["n_target_genes"] = [p[2] for p in parsed]
    obs["gem_group"] = obs["gem_group"].astype(int).astype(str)
    obs["n_guides_detected"] = obs["number_of_guides"].fillna(0).astype(int)
    # the raw guide columns are semicolon-joined per guide; keep as text
    for c in ("read_count", "UMI_count", "coverage", "guide_identity"):
        obs[c] = obs[c].astype(str)

    var = feat.copy()
    var.index = pd.Index(var.gene_name, name="gene_name")
    var = var.drop(columns=["gene_name"])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()

    # ------------------------------------------------- per-cell count summaries
    host = (adata.var.kind == "host").to_numpy()
    viral = (adata.var.kind == "viral").to_numpy()
    lenti = (adata.var.kind == "lenti").to_numpy()
    mito = host & np.char.startswith(adata.var_names.to_numpy().astype(str), "MT-")
    print(f"mitochondrial genes: {int(mito.sum())}")

    def rowsum(mask):
        return np.asarray(adata.X[:, mask].sum(axis=1)).ravel()

    adata.obs["total_umi"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["host_umi"] = rowsum(host)
    adata.obs["viral_umi"] = rowsum(viral)
    adata.obs["lenti_umi"] = rowsum(lenti)
    adata.obs["mito_umi"] = rowsum(mito)
    adata.obs["n_genes"] = np.asarray((adata.X[:, host] > 0).sum(axis=1)).ravel()
    adata.obs["viral_frac"] = adata.obs.viral_umi / adata.obs.total_umi.replace(0, np.nan)
    adata.obs["pct_mito"] = 100 * adata.obs.mito_umi / adata.obs.host_umi.replace(0, np.nan)

    # per-viral-gene counts, used by the infection classifier in 03
    for j in np.where(viral)[0]:
        adata.obs[f"vc_{adata.var_names[j]}"] = np.asarray(adata.X[:, j].todense()).ravel()

    adata.uns["source"] = {
        "accession": "GSE208240",
        "paper": "Sunshine et al., Nat Commun 14:6245 (2023)",
        "matrix_dir": str(mdir),
    }
    adata.uns["feature_kinds"] = {
        k: int(v) for k, v in adata.var.kind.value_counts().items()
    }

    # -------------------------------------------------------------- sanity set
    a = adata[adata.obs.cell_set == "assigned"]
    n_elements = a.obs.element.nunique()
    print(f"library elements with >=1 assigned cell: {n_elements}")
    assert n_elements == EXPECTED_ELEMENTS, f"expected {EXPECTED_ELEMENTS} elements, got {n_elements}"
    n_ntc_elements = a.obs.loc[a.obs.is_ntc, "element"].nunique()
    print(f"non-targeting elements: {n_ntc_elements}")

    summary = pd.DataFrame(
        [
            ("barcodes_in_matrix", n_bc),
            ("features_in_matrix", len(feat)),
            ("host_features", int(host.sum())),
            ("viral_features", int(viral.sum())),
            ("lenti_features", int(lenti.sum())),
            ("cells_assigned", int((adata.obs.cell_set == "assigned").sum())),
            ("cells_wt_candidate", int((adata.obs.cell_set == "wt_candidate").sum())),
            ("library_elements", n_elements),
            ("ntc_elements", n_ntc_elements),
            ("gem_groups", a.obs.gem_group.nunique()),
            ("median_cells_per_element", int(a.obs.element.value_counts().median())),
            ("elements_under_55_cells", int((a.obs.element.value_counts() < 55).sum())),
            ("median_host_umi_assigned", float(a.obs.host_umi.median())),
            ("median_genes_assigned", float(a.obs.n_genes.median())),
            ("median_viral_frac_assigned", float(a.obs.viral_frac.median())),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(tables / "01_build_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))

    dest = out / "perturbseq.h5ad"
    adata.write_h5ad(dest, compression="gzip")
    print(f"\nwrote {dest} ({os.path.getsize(dest) / 1e9:.2f} GB)")


if __name__ == "__main__":
    main()
