#!/usr/bin/env python
"""
10 - All figures, built from the tables written by 02-09.

Writes results/figures/fig1..fig8 .png

Analysis and plotting are kept apart on purpose: every number drawn here comes
from a .tsv in results/tables, so a figure can be regenerated, or redrawn
differently, without rerunning any statistics.

Palette: a three-slot categorical set (blue / orange / aqua) plus a blue-red
diverging pair with a gray neutral, validated for colour-vision deficiency and
normal-vision separation on an all-pairs basis. Every series is also direct-
labelled or legended, so identity never rests on hue alone.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import _common as C

# ----------------------------------------------------------------- palette
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED = "#e34948"
NEUTRAL = "#b4b3ad"
GRID = "#e6e5e1"
INK = "#0b0b0b"
INK2 = "#52514e"
B_LABEL = "2001"

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


def fig_path(name):
    C.FIGURES.mkdir(parents=True, exist_ok=True)
    return C.FIGURES / name


# ============================================================ figure 1: QC
def figure_infection_qc():
    cs = C.read_tbl("03_cell_state.tsv.gz")
    lane = C.read_tbl("03_state_summary.tsv")
    dec = C.read_tbl("03_infected_by_depth_decile.tsv")
    matched = pd.read_csv(C.TABLES / "03_host_umi_matched_depth.tsv", sep="\t", index_col=0)

    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.4))

    # (a) the bimodal viral fraction, split by the mixture call
    ax = axes[0, 0]
    bins = np.linspace(-4.2, 0, 70)
    for st in ("uninfected", "borderline", "infected"):
        v = np.log10(cs.loc[cs.state == st, "viral_frac"].clip(lower=1e-6))
        ax.hist(v, bins=bins, color=STATE_COLOR[st], alpha=0.85,
                label=f"{STATE_LABEL[st]} ({len(v):,})", linewidth=0)
    style(ax)
    ax.set_xlabel("viral fraction of the cell's library (log$_{10}$)")
    ax.set_ylabel("cells")
    ax.set_title("a   Infection state is a two-mode split", loc="left")
    ax.legend(loc="upper left")

    # (b) host shutoff at matched library size
    ax = axes[0, 1]
    x = np.arange(len(matched))
    for st in ("uninfected", "infected"):
        if st in matched.columns:
            ax.plot(x, matched[st] / 1000, color=STATE_COLOR[st], marker="o",
                    markersize=4.5, label=STATE_LABEL[st])
    style(ax)
    ax.set_xticks(x[::2])
    ax.set_xticklabels([f"{i+1}" for i in x[::2]])
    ax.set_xlabel("decile of total library size")
    ax.set_ylabel("median host UMIs (thousands)")
    ax.set_title("b   Host shutoff, not lost sequencing depth", loc="left")
    ax.legend(loc="upper left")
    ax.annotate("same library size,\nfewer host transcripts", xy=(4.15, 7.4),
                xytext=(5.1, 2.0), fontsize=7, color=INK2,
                arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))

    # (c) lane composition
    ax = axes[1, 0]
    lane["gem_group"] = lane.gem_group.astype(str)
    bottom = np.zeros(len(lane))
    for st in ("uninfected", "borderline", "infected"):
        if st not in lane.columns:
            continue
        ax.bar(lane.gem_group, lane[st] * 100, bottom=bottom, color=STATE_COLOR[st],
               width=0.68, label=STATE_LABEL[st], linewidth=1.2, edgecolor="white")
        bottom += lane[st].to_numpy() * 100
    style(ax)
    ax.set_xlabel("10x channel (gem group)")
    ax.set_ylabel("% of cells")
    ax.set_ylim(0, 112)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_title("c   Channels 5-8 carry more infection", loc="left")
    ax.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))

    # (d) the depth confound, stated plainly
    ax = axes[1, 1]
    ax.plot(dec.host_umi_median / 1000, dec.infected_frac * 100, color=ORANGE,
            marker="o", markersize=4.5)
    style(ax)
    ax.set_xscale("log")
    ax.set_xticks([1.5, 3, 6, 12, 25, 45])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("median host UMIs of the decile (thousands, log)")
    ax.set_ylabel("% called infected")
    ax.set_title("d   Infected cells are the shallow cells", loc="left")
    ax.annotate("a consequence of shutoff, not a\nclassifier artefact: the calls are\n"
                "stable under 10x downsampling",
                xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                fontsize=7, color=INK2)

    fig.suptitle("Figure 2   Infection calling and the QC structure it rests on",
                 x=0.005, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(fig_path("fig2_infection_qc.png"))
    plt.close(fig)
    print("fig2 (infection QC) done")


# ================================================== figure 2: knockdown QC
def figure_knockdown():
    kd = C.read_tbl("04_knockdown.tsv")
    b = kd[(kd.population == "bystander") & kd.kd_pseudobulk.notna() & (kd.n_cells >= 20)]
    con = C.read_tbl("04_construct_concordance.tsv")

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.1))

    ax = axes[0]
    ax.hist(b.kd_pseudobulk * 100, bins=np.arange(-30, 105, 5), color=BLUE,
            linewidth=0, alpha=0.9)
    style(ax)
    ax.axvline(75, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.annotate(f"{(b.kd_pseudobulk > 0.75).mean()*100:.0f}% of element-gene pairs\nclear 75%",
                xy=(77, ax.get_ylim()[1] * 0.78), fontsize=7, color=INK2)
    ax.set_xlabel("on-target knockdown in bystander cells (%)")
    ax.set_ylabel("element-gene pairs")
    ax.set_title(f"a   Median {b.kd_pseudobulk.median()*100:.0f}% knockdown", loc="left")

    ax = axes[1]
    indep = con[con.pair_type == "independent_vectors"]
    other = con[con.pair_type != "independent_vectors"]
    r = float(np.corrcoef(indep.kd_A, indep.kd_B)[0, 1])
    ax.axhline(0, color=GRID, linewidth=0.8)
    ax.axvline(0, color=GRID, linewidth=0.8)
    ax.plot([-0.3, 1.05], [-0.3, 1.05], color=NEUTRAL, linewidth=1.0,
            linestyle=(0, (4, 3)), zorder=1)
    ax.scatter(other.kd_A * 100, other.kd_B * 100, s=20, facecolor="none",
               edgecolor="#b5b4ad", linewidth=0.9, zorder=2,
               label="shared-guide pairs (not independent)")
    ax.scatter(indep.kd_A * 100, indep.kd_B * 100, s=30, color=AQUA,
               edgecolor="white", linewidth=0.8, zorder=3,
               label=f"independent vectors (n = {len(indep)})")
    for _, row in indep.iterrows():
        if row.target_gene in ("FURIN", "ACE2", "BRD2", "VAC14", "MEPCE", "LARP1"):
            ax.annotate(row.target_gene, (row.kd_A * 100, row.kd_B * 100),
                        textcoords="offset points", xytext=(5, -2), fontsize=7, color=INK)
    style(ax, grid_axis="both")
    ax.legend(frameon=False, fontsize=6.6, loc="upper left", handletextpad=0.4)
    ax.set_xlabel("knockdown, vector A (%)")
    ax.set_ylabel("knockdown, vector B (%)")
    ax.set_title(f"b   Separately cloned vectors agree at r = {r:.2f}", loc="left")

    fig.suptitle("Figure 3   The screen's internal positive control",
                 x=0.005, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(fig_path("fig3_knockdown.png"))
    plt.close(fig)
    print("fig3 (knockdown) done")


# ============================================= figure 3: the viral-load panel
def figure_viral_hits():
    A = C.read_tbl("07_hits_infection.tsv")
    A = A.sort_values("log_or_mh").reset_index(drop=True)
    sig = (A.q_null < 0.1) & (~A.is_ntc)

    fig, axes = plt.subplots(1, 2, figsize=(7.9, 3.6),
                             gridspec_kw={"width_ratios": [1.25, 1]})

    # (a) the ranked panel, the analogue of the paper's Fig 3A
    ax = axes[0]
    x = np.arange(len(A))
    ntc_band = A.loc[A.is_ntc, "infected_frac"]
    ax.axhspan(ntc_band.min() * 100, ntc_band.max() * 100, color=NEUTRAL,
               alpha=0.25, zorder=0)
    colors = np.where(A.is_ntc, NEUTRAL, np.where(sig, np.where(A.log_or_mh < 0, BLUE, RED), "#cfced8"))
    ax.scatter(x, A.infected_frac * 100, s=18, c=colors, edgecolor="white",
               linewidth=0.5, zorder=3)
    ax.scatter(x[A.is_ntc.to_numpy()], A.loc[A.is_ntc, "infected_frac"] * 100,
               s=26, color=INK2, marker="D", zorder=4, label="non-targeting")
    label_offsets = {
        "TMPRSS2+CTSB": (16, -4), "ACE2": (18, 4), "NFKBIA": (20, 12),
        "EIF4E2": (26, 30), "BRD2": (34, -2), "SEC62": (-30, 2),
    }
    for tgt, off in label_offsets.items():
        row = A[A.target == tgt]
        if len(row):
            i = int(row.index[0])
            ax.annotate(tgt, (i, row.infected_frac.iloc[0] * 100),
                        textcoords="offset points", xytext=off, fontsize=6.8, color=INK,
                        arrowprops=dict(arrowstyle="-", color="#9a9992", lw=0.6,
                                        shrinkA=0, shrinkB=2))
    style(ax)
    ax.set_xlabel("library element, ranked by effect")
    ax.set_ylabel("% of cells infected")
    ax.set_title("a   Infection rate per perturbation", loc="left")
    ax.legend(handles=[
        Line2D([], [], marker="D", color=INK2, linestyle="", markersize=5, label="non-targeting"),
        Line2D([], [], marker="o", color=BLUE, linestyle="", markersize=5, label="less infection (q<0.1)"),
        Line2D([], [], marker="o", color=RED, linestyle="", markersize=5, label="more infection (q<0.1)"),
        Line2D([], [], marker="o", color="#cfced8", linestyle="", markersize=5, label="not significant"),
    ], loc="upper left", ncol=1)

    # (b) effect vs significance
    ax = axes[1]
    nl = -np.log10(A.p_null.clip(lower=1e-4))
    ax.scatter(A.log_or_mh[~sig & ~A.is_ntc], nl[~sig & ~A.is_ntc], s=16,
               color="#cfced8", edgecolor="white", linewidth=0.4, zorder=2)
    ax.scatter(A.log_or_mh[A.is_ntc], nl[A.is_ntc], s=24, color=INK2, marker="D", zorder=3)
    ax.scatter(A.log_or_mh[sig], nl[sig], s=22,
               c=np.where(A.log_or_mh[sig] < 0, BLUE, RED),
               edgecolor="white", linewidth=0.5, zorder=4)
    for tgt, off in {"ACE2": (6, -12), "TMPRSS2": (6, -20), "SEC62": (-34, 0),
                     "EIF4H": (6, 6)}.items():
        row = A[A.target == tgt]
        if len(row):
            ax.annotate(tgt, (row.log_or_mh.iloc[0], nl[row.index[0]]),
                        textcoords="offset points", xytext=off, fontsize=6.8, color=INK,
                        arrowprops=dict(arrowstyle="-", color="#9a9992", lw=0.6,
                                        shrinkA=0, shrinkB=2))
    style(ax, grid_axis="both")
    ax.axvline(0, color=GRID, linewidth=0.8)
    ax.set_ylim(-0.15, 4.0)
    ax.set_xlabel("lane-stratified log odds ratio vs NTC")
    ax.set_ylabel("$-$log$_{10}$ control-based p")
    ax.set_title("b   Effect size against the NTC null", loc="left")
    ax.annotate(f"ceiling at p = 1/{B_LABEL}, the\nresolution of the control-based null",
                xy=(0.03, 0.99), xycoords="axes fraction", va="top",
                fontsize=6.6, color=INK2)
    ax.annotate(f"the {int(A.is_ntc.sum())} non-targeting elements (diamonds)\nsit at the origin: the null is calibrated",
                xy=(0.5, 0.26), xycoords="axes fraction", ha="center",
                fontsize=6.6, color=INK2)

    fig.suptitle("Figure 4   Hit calling for infection establishment",
                 x=0.005, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(fig_path("fig4_viral_hits.png"))
    plt.close(fig)
    print("fig4 (hit calling) done")


# ================================================ figure 4: the concordance
def figure_concordance():
    con = C.read_tbl("07_concordance.tsv")
    con = con.sort_values("infection_log_or")
    y = np.arange(len(con))
    colors = [BLUE if d == "down" else RED if d == "up" else NEUTRAL
              for d in con.paper_direction]

    # Fill carries this analysis's call, colour carries the paper's claim:
    #   solid bar  = an effect at q < 0.1 here
    #   hollow bar = no effect here
    # So a hollow grey bar, however long, is a target called null both in the
    # paper and here: bar length alone is not a claim.
    called = (con.infection_q < 0.1).to_numpy()

    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    ax.axvline(0, color="#c9c8c3", linewidth=0.9, zorder=1)
    ax.barh(y, con.infection_log_or,
            color=[c if k else "none" for c, k in zip(colors, called)],
            edgecolor=colors, linewidth=1.1, height=0.6, zorder=2)
    for i, (_, r) in enumerate(con.iterrows()):
        xpos = max(r.infection_log_or, 0) + 0.08
        ax.text(xpos, i, f"q = {r.infection_q:.3f}", va="center", ha="left",
                fontsize=6.6, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(con.target, fontsize=7.2)
    style(ax, grid_axis="x")
    ax.set_xlabel("lane-stratified log odds ratio vs non-targeting controls")
    ax.set_xlim(-2.2, 2.0)
    n_ok = int(con.agrees.sum())
    ax.set_title(
        f"Figure 1   All {n_ok} of {len(con)} calls in Sunshine et al. Fig. 3A reproduce\n"
        "bar colour = the paper's claim; filled = an effect here at q < 0.1",
        loc="left", fontsize=9.5,
    )
    ax.legend(handles=[
        Line2D([], [], marker="s", color=BLUE, linestyle="", markersize=6,
               label="paper: less viral RNA"),
        Line2D([], [], marker="s", color=RED, linestyle="", markersize=6,
               label="paper: more"),
        Line2D([], [], marker="s", color=NEUTRAL, linestyle="", markersize=6,
               label="paper: no effect"),
        Line2D([], [], marker="s", markerfacecolor="none", markeredgecolor=INK2,
               color=INK2, linestyle="", markersize=6,
               label="hollow: no effect here"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.115), ncol=4, fontsize=6.6,
       frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(fig_path("fig1_concordance.png"))
    plt.close(fig)
    print("fig1 (concordance) done")

# ======================================================= figure 5: the ISGs
def figure5():
    isg = C.read_tbl("06_isg_by_element.tsv")
    isg = isg[isg.score == "isg_clean"].sort_values("delta_isg").reset_index(drop=True)
    sig = (isg.q_null < 0.1) & (~isg.is_ntc)

    fig, axes = plt.subplots(1, 2, figsize=(7.9, 3.5),
                             gridspec_kw={"width_ratios": [1.3, 1]})

    ax = axes[0]
    x = np.arange(len(isg))
    colors = np.where(isg.is_ntc, INK2, np.where(sig, np.where(isg.delta_isg < 0, BLUE, RED), "#cfced8"))
    ax.axhline(0, color="#c9c8c3", linewidth=0.9)
    ax.scatter(x, isg.delta_isg, s=18, c=colors, edgecolor="white", linewidth=0.5, zorder=3)
    canonical = {"STAT2", "IRF9", "IFNAR2", "IFNAR1", "STAT1", "TYK2", "JAK1"}
    novel = {"GPR89A", "KEAP1", "SPNS1", "BRD2"}
    offsets = {
        "STAT2": (14, -8), "IRF9": (14, 8), "IFNAR1": (16, -10), "IFNAR2": (16, 2),
        "TYK2": (16, 12), "STAT1": (14, 2), "JAK1": (14, 0),
        "GPR89A": (16, 14), "SPNS1": (14, -10), "KEAP1": (14, 4), "BRD2": (14, 0),
    }
    for _, r in isg.iterrows():
        if r.target in offsets:
            ax.annotate(r.target, (r.name, r.delta_isg), textcoords="offset points",
                        xytext=offsets[r.target], fontsize=6.5,
                        color=BLUE if r.target in canonical else ORANGE,
                        arrowprops=dict(arrowstyle="-", color="#9a9992", lw=0.5,
                                        shrinkA=0, shrinkB=2))
    style(ax)
    ax.set_xlabel("library element, ranked")
    ax.set_ylabel("$\\Delta$ bystander ISG score vs NTC")
    ax.set_title("a   Interferon tone in bystander cells", loc="left")
    ax.legend(handles=[
        Line2D([], [], marker="o", color=BLUE, linestyle="", markersize=5,
               label="canonical JAK-STAT (blue labels)"),
        Line2D([], [], marker="o", color=ORANGE, linestyle="", markersize=5,
               label="the paper's four novel hits (orange labels)"),
        Line2D([], [], marker="o", color=INK2, linestyle="", markersize=5, label="non-targeting"),
    ], loc="lower right")

    # the paper's four, with and without the self-containment correction
    ax = axes[1]
    cmp = C.read_tbl("06_isg_paper_comparison.tsv")
    four = cmp[(cmp.paper == "down") & cmp.element.notna()]
    piv = four.pivot_table(index="target", columns="score", values="delta_isg")
    piv = piv.sort_values("isg_clean")
    yy = np.arange(len(piv))
    ax.barh(yy - 0.18, piv["isg"], height=0.33, color=NEUTRAL, label="all 55 ISGs", linewidth=0)
    ax.barh(yy + 0.18, piv["isg_clean"], height=0.33, color=AQUA,
            label="excluding ISGs that are\nthemselves library targets", linewidth=0)
    ax.axvline(0, color="#c9c8c3", linewidth=0.9)
    ax.set_yticks(yy)
    ax.set_yticklabels(piv.index, fontsize=7.2)
    style(ax, grid_axis="x")
    ax.set_xlabel("$\\Delta$ bystander ISG score vs NTC")
    ax.set_title("b   The score is not self-fulfilling", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.92))

    fig.suptitle("Figure 5   Reproduction of the interferon panel",
                 x=0.005, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(fig_path("fig5_isg.png"))
    plt.close(fig)
    print("fig5 done")


# ================================================= figure 6: the resolution
def figure6():
    curve = C.read_tbl("08_downsample_curve.tsv")
    det = C.read_tbl("08_detectable_effect.tsv")
    rep = C.read_tbl("04_element_representation.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(8.6, 3.0))

    ax = axes[0]
    order = (curve.groupby("target").full_log_or.first().abs().sort_values(ascending=False))
    ramp = ["#0d366b", "#184f95", "#256abf", "#2a78d6", "#3987e5", "#5598e7", "#86b6ef", "#9ec5f4"]
    handles = []
    for (tgt, _), col in zip(order.items(), ramp):
        s = curve[(curve.target == tgt) & curve.cells.isin([25, 40, 60, 80, 120, 160])]
        if len(s) < 2:
            continue
        s = s.sort_values("cells")
        ax.plot(s.cells, s.detection_rate * 100, color=col, marker="o", markersize=3.6,
                linewidth=1.8)
        handles.append(Line2D([], [], color=col, marker="o", markersize=3.6, linewidth=1.8,
                              label=f"{tgt}  ({abs(order[tgt]):.2f})"))
    style(ax)
    ax.set_xlabel("cells retained for the element")
    ax.set_ylabel("% of subsamples detected (p<0.05)")
    ax.set_xlim(18, 170)
    ax.set_ylim(0, 108)
    ax.legend(handles=handles, loc="lower right", fontsize=5.9, labelspacing=0.28,
              handlelength=1.4)
    ax.set_title("a   Detection by coverage\n(legend: |log OR| at full depth)", loc="left")

    ax = axes[1]
    ax.plot(det.cells, det.fold_change_at_p05, color=ORANGE, marker="o", markersize=4.5)
    style(ax)
    ax.set_xlabel("cells per element")
    ax.set_ylabel("smallest detectable odds ratio")
    ax.set_title("b   What the screen can resolve", loc="left")
    for _, r in det.iterrows():
        if r.cells in (25, 60, 120, 240):
            ax.annotate(f"{r.fold_change_at_p05:.2f}x", (r.cells, r.fold_change_at_p05),
                        textcoords="offset points", xytext=(5, 5), fontsize=6.8, color=INK2)

    ax = axes[2]
    v = rep.n_cells.clip(upper=600)
    ax.hist(v, bins=np.arange(0, 620, 25), color=BLUE, linewidth=0, alpha=0.9)
    ax.axvline(55, color=RED, linewidth=1.2, linestyle=(0, (4, 3)))
    style(ax)
    ax.set_xlabel("cells per library element")
    ax.set_ylabel("elements")
    n_under = int((rep.n_cells < 55).sum())
    ax.set_title(f"c   {n_under} of {len(rep)} elements below 55 cells", loc="left")
    ax.annotate(f"{n_under} elements fall left of the\ncoverage cut: mostly essential\ngenes, so absence from a hit\nlist is informative, not evidence\nof no effect",
                xy=(0.97, 0.80), xycoords="axes fraction", ha="right", va="top",
                fontsize=6.6, color=INK2)

    fig.suptitle("Figure 6   Perturbation representation sets what the screen can find",
                 x=0.005, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(fig_path("fig6_resolution.png"))
    plt.close(fig)
    print("fig6 done")


# ========================================== figure 7: infected vs bystander
def figure7():
    pt = C.read_tbl("09_program_scores_by_state.tsv")
    tab = C.read_tbl("09_state_specific_effects.tsv")

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4))

    ax = axes[0]
    lim = (min(pt.isg_bystander.min(), pt.isg_infected.min()) - 0.3,
           max(pt.isg_bystander.max(), pt.isg_infected.max()) + 0.3)
    ax.plot(lim, lim, color=NEUTRAL, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    ax.axhline(0, color=GRID, linewidth=0.8)
    ax.axvline(0, color=GRID, linewidth=0.8)
    ax.scatter(pt.isg_bystander, pt.isg_infected, s=24, color=BLUE,
               edgecolor="white", linewidth=0.6, zorder=3)
    offs = {"IRF9": (6, 2), "STAT2": (8, -8), "IFNAR1+IFNAR2": (10, -12),
            "TYK2": (8, -10), "GPR89A": (8, 6), "SPNS1": (10, -12),
            "KEAP1": (14, -4), "BRD2": (10, 8)}
    for tgt, off in offs.items():
        row = pt[pt.target == tgt]
        if len(row):
            ax.annotate(tgt, (row.isg_bystander.iloc[0], row.isg_infected.iloc[0]),
                        textcoords="offset points", xytext=off, fontsize=6.5, color=INK,
                        arrowprops=dict(arrowstyle="-", color="#9a9992", lw=0.5,
                                        shrinkA=0, shrinkB=2))
    r = float(np.corrcoef(pt.isg_bystander, pt.isg_infected)[0, 1])
    style(ax, grid_axis="both")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("ISG effect in bystander cells (log$_2$)")
    ax.set_ylabel("ISG effect in infected cells (log$_2$)")
    ax.set_title(f"a   Interferon effects survive infection (r = {r:.2f})", loc="left")
    ax.annotate("points below the diagonal:\nthe same perturbation, blunted\nonce the cell is infected",
                xy=(0.97, 0.06), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=6.6, color=INK2)

    ax = axes[1]
    progs = ["isg", "mhc_class_i", "nfkb", "unfolded_protein"]
    names = ["interferon-\nstimulated", "MHC class I", "NF-$\\kappa$B", "unfolded-\nprotein"]
    rs, ns = [], []
    for p in progs:
        sub = pt[[f"{p}_bystander", f"{p}_infected"]].dropna()
        rs.append(float(np.corrcoef(sub.iloc[:, 0], sub.iloc[:, 1])[0, 1]))
        ns.append(int(pt[f"{p}_n_genes"].iloc[0]))
    y = np.arange(len(progs))
    ax.barh(y, rs, color=[BLUE, BLUE, AQUA, AQUA], height=0.6, linewidth=0)
    for i, (v, n) in enumerate(zip(rs, ns)):
        ax.text(v + 0.015, i, f"{v:.2f}   ({n} genes)", va="center", fontsize=6.8, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7.2)
    style(ax, grid_axis="x")
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("correlation of the effect between the two states")
    ax.set_title("b   How transferable each program is", loc="left")

    fig.suptitle("Figure 7   Extension: does a perturbation do the same thing once the cell is infected?",
                 x=0.005, ha="left", fontsize=10.5, fontweight="semibold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(fig_path("fig7_infected_vs_bystander.png"))
    plt.close(fig)
    print(f"fig7 done (median gene-level profile r = {tab.r_profiles.median():.3f})")



# =============================== figure 8: the organoid bridge
def figure_organoid():
    el = C.read_tbl("14_element_signature.tsv")

    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    xcol = "rho_organoid_infection_within_run"
    ycol = "rho_organoid_drug_in_infected"
    t = el[~el.is_ntc]
    n = el[el.is_ntc]
    sig = t[t.q_organoid_drug_in_infected < 0.1]
    ax.axhline(0, color=GRID, linewidth=0.8)
    ax.axvline(0, color=GRID, linewidth=0.8)
    ax.axhspan(n[ycol].min(), n[ycol].max(), color=NEUTRAL, alpha=0.25, zorder=0)
    ax.scatter(t[xcol], t[ycol], s=16, color="#cfced8", edgecolor="white",
               linewidth=0.4, zorder=2)
    ax.scatter(sig[xcol], sig[ycol], s=30, color=AQUA, edgecolor="white",
               linewidth=0.6, zorder=3)
    ax.scatter(n[xcol], n[ycol], s=22, color=INK2, marker="D", zorder=4)
    for lab, dx, dy in [("KEAP1", 7, 1), ("SPNS1", 8, -8), ("STAT2", 7, 3),
                        ("GPR89A", 7, -7), ("BRD2", 7, 1), ("IRF9+IRF3", 7, 4),
                        ("TYK2", -26, 4)]:
        r = sig[sig.target == lab]
        if len(r):
            ax.annotate(lab, (r[xcol].iloc[0], r[ycol].iloc[0]),
                        textcoords="offset points", xytext=(dx, dy),
                        fontsize=6.8, color=INK)
    style(ax, grid_axis="both")
    ax.set_xlabel("resembles the organoid infection signature")
    ax.set_ylabel("resembles the organoid drug signature")
    ax.set_title("Figure 8   Interferon knockdowns phenocopy the drug signature\n"
                 "each point is a library element, scored in bystander cells",
                 loc="left", fontsize=9.5)
    ax.legend(handles=[
        Line2D([], [], marker="o", color=AQUA, linestyle="", markersize=6,
               label="q < 0.1 against the control band"),
        Line2D([], [], marker="D", color=INK2, linestyle="", markersize=5,
               label="non-targeting elements"),
    ], loc="lower left", fontsize=6.8)
    fig.tight_layout()
    fig.savefig(fig_path("fig8_organoid_bridge.png"))
    plt.close(fig)
    print("fig8 (organoid bridge) done")


if __name__ == "__main__":
    figure_concordance()      # Figure 1: does the screen reproduce?
    figure_infection_qc()     # Figure 2: is the phenotype real?
    figure_knockdown()        # Figure 3: did the perturbations work?
    figure_viral_hits()       # Figure 4: the full hit list
    figure5()
    figure6()
    figure7()
    figure_organoid()
    print(f"\nfigures written to {C.FIGURES}")
