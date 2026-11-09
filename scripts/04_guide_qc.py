#!/usr/bin/env python
"""
04 - Guide/element QC: representation, and on-target knockdown.

Writes
    results/tables/04_knockdown.tsv              knockdown per element x target gene
    results/tables/04_knockdown_summary.tsv      the library-level headline numbers
    results/tables/04_element_representation.tsv cells and lanes per element

On-target knockdown is the screen's internal positive control: a CRISPRi
element that does not lower its own target's mRNA cannot be interpreted, and
including it in a hit list inflates the false-positive rate.

Two decisions matter and both are reported rather than assumed.

1. Which cells.  Knockdown is measured in bystander cells only (state
   `uninfected` from 03). Infected cells lose most of their host mRNA to
   NSP1-driven shutoff, so including them compresses every ratio toward zero
   and makes knockdown look better than it is. The table reports the bystander
   estimate and the all-cells estimate side by side.

2. Which estimator.  The published figure quotes a single percentage per
   target, and the natural reading, "normalised target counts relative to
   NTCs", can mean either
       mean-of-ratios : mean over cells of (target CPM), divided by the NTC mean
       pseudobulk     : summed target counts / summed host counts, over the
                        element's cells, divided by the same ratio in NTC cells
   These are not the same statistic when depth varies between cells, and the
   pseudobulk version is the one that is robust to a few very shallow cells.
   Both are reported. The repository quotes the pseudobulk version.

For the 22 two-gene elements, knockdown is reported for each gene separately.
Elements 216-239 are 12 genes carried by two independent single-guide vectors;
their two estimates are an internal reproducibility check and are flagged with
`has_second_element`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common as C


def main() -> None:
    a = C.load(assigned_only=True)
    state = C.read_tbl("03_cell_state.tsv.gz").set_index("cell_barcode")
    a.obs["state"] = state.reindex(a.obs_names)["state"].to_numpy()

    host = (a.var.kind == "host").to_numpy()
    host_umi = a.obs.host_umi.to_numpy().astype(float)
    host_umi[host_umi == 0] = np.nan

    # ------------------------------------------------- target genes per element
    el = (
        a.obs.groupby("element", observed=True)
        .agg(target=("target", "first"), is_ntc=("is_ntc", "first"),
             n_cells=("host_umi", "size"), n_lanes=("gem_group", "nunique"))
        .reset_index()
    )
    pairs = []
    for _, r in el.iterrows():
        if r.is_ntc:
            continue
        for g in r.target.split("+"):
            pairs.append((r.element, g))
    genes = sorted({g for _, g in pairs})
    present = [g for g in genes if g in a.var_names]
    missing = [g for g in genes if g not in a.var_names]
    print(f"target genes: {len(genes)}; present in the matrix: {len(present)}")
    if missing:
        print(f"  not in the reference (reported as NA): {missing}")

    # pull the target-gene columns once
    idx = [a.var_names.get_loc(g) for g in present]
    M = np.asarray(a.X[:, idx].todense())           # cells x target genes
    cpm = M / host_umi[:, None] * 1e6
    col = {g: j for j, g in enumerate(present)}

    is_ntc = a.obs.is_ntc.to_numpy()
    byst = (a.obs.state == "uninfected").to_numpy()

    def stats(mask, j):
        """mean-of-CPM and pseudobulk rate for gene j over the cells in mask."""
        if mask.sum() == 0:
            return np.nan, np.nan, 0
        mean_cpm = float(np.nanmean(cpm[mask, j]))
        pb = float(M[mask, j].sum() / np.nansum(host_umi[mask]) * 1e6)
        return mean_cpm, pb, int(mask.sum())

    rows = []
    for element, g in pairs:
        if g not in col:
            rows.append(dict(element=element, target_gene=g, note="gene absent from reference"))
            continue
        j = col[g]
        m_el = (a.obs.element == element).to_numpy()
        for pop, cells_el, cells_ntc in [
            ("bystander", m_el & byst, is_ntc & byst),
            ("all_cells", m_el, is_ntc),
        ]:
            e_mean, e_pb, n_e = stats(cells_el, j)
            c_mean, c_pb, n_c = stats(cells_ntc, j)
            rows.append(
                dict(
                    element=element,
                    target_gene=g,
                    population=pop,
                    n_cells=n_e,
                    n_ntc_cells=n_c,
                    element_mean_cpm=e_mean,
                    ntc_mean_cpm=c_mean,
                    element_pseudobulk_cpm=e_pb,
                    ntc_pseudobulk_cpm=c_pb,
                    kd_mean_of_cpm=1 - e_mean / c_mean if c_mean else np.nan,
                    kd_pseudobulk=1 - e_pb / c_pb if c_pb else np.nan,
                )
            )
    kd = pd.DataFrame(rows)

    # flag the 12 genes with two independent elements
    tgt_counts = (
        kd[kd.population == "bystander"].groupby("target_gene").element.nunique()
    )
    kd["has_second_element"] = kd.target_gene.map(tgt_counts.gt(1))
    kd = kd.merge(el[["element", "n_cells", "n_lanes"]].rename(
        columns={"n_cells": "element_n_cells"}), on="element", how="left")
    kd.to_csv(C.tbl("04_knockdown.tsv"), sep="\t", index=False)

    # ----------------------------------------------------------- the headline
    b = kd[(kd.population == "bystander") & kd.kd_pseudobulk.notna()]
    # only elements with enough bystander cells to estimate anything
    bb = b[b.n_cells >= 20]
    summ = pd.DataFrame(
        [
            ("element_gene_pairs_tested", len(b)),
            ("pairs_with_20plus_bystander_cells", len(bb)),
            ("median_knockdown_pseudobulk", float(bb.kd_pseudobulk.median())),
            ("median_knockdown_mean_of_cpm", float(bb.kd_mean_of_cpm.median())),
            ("frac_above_50pct_kd", float((bb.kd_pseudobulk > 0.5).mean())),
            ("frac_above_75pct_kd", float((bb.kd_pseudobulk > 0.75).mean())),
            ("frac_above_90pct_kd", float((bb.kd_pseudobulk > 0.90).mean())),
            ("n_negative_kd", int((bb.kd_pseudobulk < 0).sum())),
        ],
        columns=["metric", "value"],
    )
    summ.to_csv(C.tbl("04_knockdown_summary.tsv"), sep="\t", index=False)
    print("\nknockdown (bystander cells, pseudobulk estimator):")
    print(summ.to_string(index=False))

    allc = kd[(kd.population == "all_cells") & kd.kd_pseudobulk.notna() & (kd.n_cells >= 20)]
    print(
        f"\nsame statistic computed on all cells instead of bystanders: median "
        f"{allc.kd_pseudobulk.median():.3f} (bystander {bb.kd_pseudobulk.median():.3f}), "
        f"fraction above 75%: {(allc.kd_pseudobulk > 0.75).mean():.3f} "
        f"(bystander {(bb.kd_pseudobulk > 0.75).mean():.3f})"
    )

    print("\nweakest knockdown (bystander, >=20 cells):")
    print(
        bb.nsmallest(12, "kd_pseudobulk")[
            ["element", "target_gene", "n_cells", "ntc_pseudobulk_cpm",
             "element_pseudobulk_cpm", "kd_pseudobulk"]
        ].round(3).to_string(index=False)
    )

    # ------------------------------------------- construct concordance
    # Two elements hitting the same gene are NOT automatically two independent
    # constructs, and conflating the two categories inflates the concordance.
    # The library gives two genuinely different situations:
    #
    #   independent_vectors   elements 216-239: 24 separately cloned
    #                         single-guide vectors covering 12 genes, two each.
    #                         These are independent constructs: different guide
    #                         sequences, separate clones. This is the honest
    #                         reagent-reproducibility panel.
    #   single_vs_combination a single-gene element and a two-gene combination
    #                         element that contains the same gene (e.g.
    #                         CRISPRi_001_ACE2 vs CRISPRi_172_ACE2_TMPRSS2).
    #                         These typically share guide design for the shared
    #                         gene, so agreement between them is a weaker claim
    #                         and partly measures the same reagent twice.
    #
    # Both are reported, separately and labelled. The repository quotes the
    # 12-gene independent-vector number as the reproducibility figure.
    bb = bb.copy()
    bb["element_number"] = bb.element.str.split("_").str[1].astype(int)
    bb["element_is_combination"] = bb.element.map(
        lambda e: "+" in el.set_index("element").target.get(e, "")
    )

    pair_rows = []
    for gene, sub in bb.groupby("target_gene"):
        if sub.target_gene.count() < 2:
            continue
        sub = sub.sort_values("element_number")
        if len(sub) != 2:
            kind = "more_than_two_elements"
        elif sub.element_number.between(216, 239).all() and not sub.element_is_combination.any():
            kind = "independent_vectors"
        elif sub.element_is_combination.any():
            kind = "single_vs_combination"
        else:
            kind = "two_single_gene_elements"
        pair_rows.append(
            dict(
                target_gene=gene,
                pair_type=kind,
                element_A=sub.element.iloc[0],
                element_B=sub.element.iloc[1],
                kd_A=sub.kd_pseudobulk.iloc[0],
                kd_B=sub.kd_pseudobulk.iloc[1],
                n_cells_A=int(sub.n_cells.iloc[0]),
                n_cells_B=int(sub.n_cells.iloc[1]),
                abs_difference=abs(sub.kd_pseudobulk.iloc[0] - sub.kd_pseudobulk.iloc[1]),
            )
        )
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(C.tbl("04_construct_concordance.tsv"), sep="\t", index=False)

    print("\nconstruct concordance, by what the two elements actually are:")
    summary_rows = []
    for kind, sub in pairs.groupby("pair_type"):
        if len(sub) < 3:
            print(f"  {kind}: n={len(sub)} (too few to correlate)")
            continue
        r = float(np.corrcoef(sub.kd_A, sub.kd_B)[0, 1])
        rho = float(sub.kd_A.corr(sub.kd_B, method="spearman"))
        print(
            f"  {kind}: n={len(sub)}  Pearson r={r:.3f}  Spearman rho={rho:.3f}  "
            f"median |difference| = {sub.abs_difference.median():.3f}"
        )
        summary_rows.append(dict(pair_type=kind, n_genes=len(sub), pearson_r=r,
                                 spearman_rho=rho,
                                 median_abs_difference=float(sub.abs_difference.median())))
    pd.DataFrame(summary_rows).to_csv(
        C.tbl("04_construct_concordance_summary.tsv"), sep="\t", index=False)

    indep = pairs[pairs.pair_type == "independent_vectors"]
    print(
        f"\nthe reproducibility panel: {len(indep)} genes carried by two separately "
        "cloned single-guide vectors (library elements 216-239)"
    )
    print(indep[["target_gene", "element_A", "element_B", "kd_A", "kd_B",
                 "abs_difference"]].round(3).to_string(index=False))
    print(
        "\nPooling all pair types into one correlation inflates it: the genes "
        "appearing in three or more elements agree almost perfectly with each other "
        "and drag the pooled number up. Only the independent-vector row is a "
        "reagent-to-reagent comparison, and it is the one quoted."
    )

    # ----------------------------------------------------- representation
    rep_tbl = (
        a.obs.groupby("element", observed=True)
        .agg(
            target=("target", "first"),
            is_ntc=("is_ntc", "first"),
            n_cells=("host_umi", "size"),
            n_bystander=("state", lambda s: int((s == "uninfected").sum())),
            n_infected=("state", lambda s: int((s == "infected").sum())),
            n_lanes=("gem_group", "nunique"),
            host_umi_median=("host_umi", "median"),
        )
        .reset_index()
        .sort_values("n_cells")
    )
    rep_tbl.to_csv(C.tbl("04_element_representation.tsv"), sep="\t", index=False)
    print(
        f"\nelements with >=20 bystander and >=20 infected cells: "
        f"{int(((rep_tbl.n_bystander >= 20) & (rep_tbl.n_infected >= 20)).sum())} of {len(rep_tbl)}"
    )


if __name__ == "__main__":
    main()
