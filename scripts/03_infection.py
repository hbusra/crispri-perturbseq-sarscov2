#!/usr/bin/env python
"""
03 - Infection-state classification, three ways, and how much of it is depth.

Writes
    results/tables/03_cell_state.tsv.gz          per-cell state, merged by later scripts
    results/tables/03_ambient_background.tsv     ambient viral rate per viral gene
    results/tables/03_classifier_agreement.tsv   concordance of the three rules
    results/tables/03_state_summary.tsv          state composition per lane
    results/tables/03_host_umi_matched_depth.tsv host shutoff at matched library size
    results/tables/03_depth_sensitivity.tsv      calls after downsampling
    results/tables/03_infected_by_depth_decile.tsv

Why this is step 03 and not step 04
-----------------------------------
On-target knockdown has to be measured in bystander cells, because infection
collapses host transcription. The infection call is therefore upstream of the
guide QC in 04, not downstream of it.

The ambient-RNA problem
-----------------------
Every droplet picks up cell-free viral RNA, so a viral count above zero is not
evidence of infection. The paper spiked in ~1% wild-type Calu-3 cells (no
dCas9, no sgRNA) to measure that background. The reference is rebuilt as

    barcodes with no guide detected, zero lentiviral reads, total UMI >= 10,000

which gives 705 cells here. That set is *not* clean: about an eighth of it sits
in the high-viral mode, i.e. real library cells whose guide was not captured.
A mean-and-SD background computed on it is therefore inflated by more than an
order of magnitude (SCV_N: mean 966, SD 2487), and a mean + 2 SD rule calls
almost nothing infected. Median-based statistics are used instead, which tolerate
that contamination, and the resulting rule is checked against the reference set.

Three classifiers
-----------------
1. `mixture`  (primary) - two-component Gaussian mixture on log10 viral
   fraction. The distribution is clearly bimodal (ambient mode near 3e-3,
   infected mode near 0.5), so the mixture finds the split without a
   hand-chosen cutoff. Posterior > 0.9 is infected, < 0.1 uninfected, the rest
   borderline.
2. `ambient`  - per viral gene, infected if the count exceeds the ambient
   expectation for that cell's library size: count_g >= max(6, k * r_g * total_umi)
   where r_g is the median ambient rate in the reference set. This is the
   depth-aware version of the paper's per-gene rule.
3. `fixed`    - viral fraction above the antimode of the bimodal distribution.
   The crude rule, included so the reader can see it barely matters.

The three are reported against each other. Downstream scripts use `mixture`,
and `07_hit_calling.py` repeats its headline test under `ambient` to show the
hit list does not depend on the choice.

The depth check
---------------
A rule with an absolute count floor calls deeply sequenced cells infected more
often, so a perturbation that changes library depth can move the apparent
infection rate with no change in infection biology. Two diagnostics: the
infected fraction across deciles of host depth, and the stability of the calls
when viral counts are binomially downsampled.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

import _common as C

MIN_COUNT = 6              # absolute floor on viral counts for one gene
AMBIENT_K = 5.0            # fold over the depth-matched ambient expectation
REF_MIN_TOTAL_UMI = 10_000
P_HI, P_LO = 0.9, 0.1      # mixture posteriors for infected / uninfected


def fit_mixture(viral_frac: pd.Series) -> tuple[pd.Series, dict]:
    """Two-component Gaussian mixture on log10 viral fraction."""
    x = np.log10(viral_frac.clip(lower=1e-6)).to_numpy().reshape(-1, 1)
    gm = GaussianMixture(n_components=2, random_state=0, n_init=5).fit(x)
    hi = int(np.argmax(gm.means_.ravel()))
    post = gm.predict_proba(x)[:, hi]
    info = dict(
        mode_low=float(10 ** gm.means_.ravel()[1 - hi]),
        mode_high=float(10 ** gm.means_.ravel()[hi]),
        weight_high=float(gm.weights_[hi]),
    )
    # the crossover point, for the `fixed` rule
    grid = np.linspace(x.min(), x.max(), 4000).reshape(-1, 1)
    gp = gm.predict_proba(grid)[:, hi]
    info["antimode_viral_frac"] = float(10 ** grid[np.argmin(np.abs(gp - 0.5))][0])
    return pd.Series(post, index=viral_frac.index, name="p_infected"), info


def main() -> None:
    a = C.load(assigned_only=False)
    obs = a.obs
    viral = [g for g in C.VIRAL_GENES if f"vc_{g}" in obs.columns]
    vcols = [f"vc_{g}" for g in viral]

    # ---------------------------------------------------- ambient reference
    ref = obs[
        (obs.cell_set == "wt_candidate")
        & (obs.lenti_umi == 0)
        & (obs.total_umi >= REF_MIN_TOTAL_UMI)
    ]
    print(f"ambient reference cells: {len(ref)}")
    if len(ref) < 100:
        raise SystemExit("too few reference cells to estimate the ambient background")

    bg = pd.DataFrame(index=pd.Index(viral, name="viral_gene"))
    bg["mean_ref"] = [ref[c].mean() for c in vcols]
    bg["sd_ref"] = [ref[c].std(ddof=1) for c in vcols]
    bg["median_ref"] = [ref[c].median() for c in vcols]
    # median rate per UMI: robust to the infected cells contaminating the set
    bg["ambient_rate"] = [float(np.median(ref[c] / ref.total_umi)) for c in vcols]
    bg.to_csv(C.tbl("03_ambient_background.tsv"), sep="\t")
    print("\nambient background per viral gene:\n", bg.round(6).to_string())
    print(
        "\nnote the mean/SD inflation: a mean + 2 SD rule would put the SCV_N "
        f"threshold at {bg.loc['SCV_N', 'mean_ref'] + 2 * bg.loc['SCV_N', 'sd_ref']:.0f} counts, "
        f"versus a median ambient rate of {bg.loc['SCV_N', 'ambient_rate']:.5f} per UMI."
    )

    # ------------------------------------------------------- the three rules
    # 1. mixture on log10 viral fraction, fitted on the assigned cells
    asg_mask = obs.cell_set == "assigned"
    post_all, info = fit_mixture(obs.viral_frac.fillna(1e-6))
    print(
        f"\nmixture: ambient mode at viral_frac {info['mode_low']:.4f}, "
        f"infected mode at {info['mode_high']:.3f}, "
        f"crossover at {info['antimode_viral_frac']:.4f}, "
        f"high-component weight {info['weight_high']:.3f}"
    )
    state_mix = np.where(post_all >= P_HI, "infected",
                         np.where(post_all <= P_LO, "uninfected", "borderline"))

    # 2. depth-aware per-gene ambient rule
    vc = obs[vcols].to_numpy()
    lam = np.outer(obs.total_umi.to_numpy(), bg.ambient_rate.to_numpy())
    thr = np.maximum(MIN_COUNT, AMBIENT_K * lam)
    state_amb = np.where((vc >= thr).any(axis=1), "infected", "uninfected")

    # 3. fixed cutoff at the mixture crossover
    cut = info["antimode_viral_frac"]
    state_fix = np.where(obs.viral_frac.fillna(0) > cut, "infected", "uninfected")

    # 4. the published rule, implemented literally: per viral gene, infected if
    #    the count clears both an absolute floor of 6 and mean + 2 SD of the
    #    reference. Included to document why it is not used: the reference set
    #    that can be rebuilt from the deposited data is contaminated by real
    #    library cells whose guide was not captured, so its SD is enormous
    #    (SCV_N: mean 966, SD 2487) and the threshold lands above almost every
    #    cell. This is a limitation of the reconstructable reference, not a
    #    claim that the published rule is wrong -- the authors had the true
    #    spike-in identity, which is not recoverable from the deposited matrix.
    thr_literal = np.maximum(MIN_COUNT, bg.mean_ref + 2 * bg.sd_ref).to_numpy()
    state_lit = np.where((obs[vcols].to_numpy() >= thr_literal[None, :]).any(axis=1),
                         "infected", "uninfected")

    obs = obs.assign(state=state_mix, p_infected=post_all.to_numpy(),
                     state_ambient=state_amb, state_fixed=state_fix,
                     state_paper_literal=state_lit)
    asg = obs[asg_mask]

    # ------------------------------------------------------------ concordance
    agree = pd.DataFrame(
        [
            dict(rule="mixture (primary)",
                 infected_frac=float((asg.state == "infected").mean()),
                 borderline_frac=float((asg.state == "borderline").mean()),
                 ref_called_infected=float((obs.loc[ref.index, "state"] == "infected").mean())),
            dict(rule="ambient per-gene",
                 infected_frac=float((asg.state_ambient == "infected").mean()),
                 borderline_frac=0.0,
                 ref_called_infected=float((obs.loc[ref.index, "state_ambient"] == "infected").mean())),
            dict(rule=f"fixed viral_frac > {cut:.3f}",
                 infected_frac=float((asg.state_fixed == "infected").mean()),
                 borderline_frac=0.0,
                 ref_called_infected=float((obs.loc[ref.index, "state_fixed"] == "infected").mean())),
            dict(rule="published rule, literal (mean + 2 SD of reference)",
                 infected_frac=float((asg.state_paper_literal == "infected").mean()),
                 borderline_frac=0.0,
                 ref_called_infected=float(
                     (obs.loc[ref.index, "state_paper_literal"] == "infected").mean())),
        ]
    )
    # pairwise agreement on the non-borderline cells
    nb = asg[asg.state != "borderline"]
    agree["agree_with_mixture"] = [
        1.0,
        float((nb.state_ambient == nb.state).mean()),
        float((nb.state_fixed == nb.state).mean()),
        float((nb.state_paper_literal == nb.state).mean()),
    ]
    agree.to_csv(C.tbl("03_classifier_agreement.tsv"), sep="\t", index=False)
    print("\nclassifier comparison (assigned cells):\n", agree.round(4).to_string(index=False))

    # ---------------------------------------------------------- state exports
    out = asg[["element", "target", "is_ntc", "gem_group", "state", "p_infected",
               "state_ambient", "state_fixed", "state_paper_literal", "host_umi",
               "total_umi", "viral_umi", "viral_frac", "n_genes", "pct_mito"]].copy()
    out.index.name = "cell_barcode"
    out.to_csv(C.tbl("03_cell_state.tsv.gz"), sep="\t", compression="gzip")

    comp = asg.state.value_counts().rename("n_cells").to_frame()
    comp["fraction"] = asg.state.value_counts(normalize=True)
    print("\nassigned-cell state composition:\n", comp.round(4).to_string())

    per_lane = pd.crosstab(asg.gem_group, asg.state, normalize="index").reset_index()
    per_lane["n_cells"] = asg.gem_group.value_counts().reindex(per_lane.gem_group).to_numpy()
    per_lane.to_csv(C.tbl("03_state_summary.tsv"), sep="\t", index=False)
    print("\nstate fractions per lane:\n", per_lane.round(3).to_string(index=False))
    print(
        "\nLanes 1-4 and 5-8 differ systematically in viral load (see 02). The paper "
        "combined a pulsed and a non-pulsed infection condition after finding no "
        "significant difference; that split is what the lane effect reflects, which is "
        "why gem_group stays in every model in 07."
    )

    # ------------------------------------------- host shutoff vs lost depth
    inf = asg[asg.state == "infected"]
    byst = asg[asg.state == "uninfected"]
    print(
        f"\nhost UMI median: infected {inf.host_umi.median():.0f} vs bystander "
        f"{byst.host_umi.median():.0f} (ratio {inf.host_umi.median() / byst.host_umi.median():.2f})"
    )
    print(
        f"total UMI median: infected {inf.total_umi.median():.0f} vs bystander "
        f"{byst.total_umi.median():.0f} (ratio {inf.total_umi.median() / byst.total_umi.median():.2f})"
    )
    matched = (
        asg.assign(bin=pd.qcut(asg.total_umi, 10, duplicates="drop"))
        .pivot_table(index="bin", columns="state", values="host_umi",
                     aggfunc="median", observed=True)
    )
    matched.to_csv(C.tbl("03_host_umi_matched_depth.tsv"), sep="\t")
    print(
        "\nmedian host UMI within deciles of total library size "
        "(a gap here is composition, i.e. host shutoff, not sequencing depth):\n",
        matched.round(0).to_string(),
    )

    # ---------------------------------------------------- depth sensitivity
    rng = np.random.default_rng(0)
    rows = []
    base = asg.state.to_numpy()
    for frac in [1.0, 0.75, 0.5, 0.25, 0.1]:
        if frac == 1.0:
            vf = asg.viral_frac.fillna(1e-6)
        else:
            vsub = rng.binomial(asg[vcols].to_numpy().astype(int), frac).sum(axis=1)
            hsub = rng.binomial(asg.host_umi.to_numpy().astype(int), frac)
            vf = pd.Series(vsub / np.maximum(vsub + hsub, 1), index=asg.index).clip(lower=1e-6)
        st = np.where(vf > cut, "infected", "uninfected")
        rows.append(
            dict(
                downsample_frac=frac,
                infected_frac=float((st == "infected").mean()),
                agree_with_full_depth=float(
                    (st[base != "borderline"] == base[base != "borderline"]).mean()
                ),
            )
        )
    ds = pd.DataFrame(rows)
    ds.to_csv(C.tbl("03_depth_sensitivity.tsv"), sep="\t", index=False)
    print("\nbinomial downsampling of the whole library:\n", ds.round(4).to_string(index=False))

    dec = (
        asg.assign(decile=pd.qcut(asg.host_umi, 10, labels=False, duplicates="drop"))
        .groupby("decile", observed=True)
        .agg(
            host_umi_median=("host_umi", "median"),
            infected_frac=("state", lambda s: float((s == "infected").mean())),
            borderline_frac=("state", lambda s: float((s == "borderline").mean())),
            n=("state", "size"),
        )
        .reset_index()
    )
    dec.to_csv(C.tbl("03_infected_by_depth_decile.tsv"), sep="\t", index=False)
    print("\nby host-depth decile:\n", dec.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
