#!/usr/bin/env python
"""
12 - What are the cells between the two modes?

Writes
    results/tables/12_borderline_by_lane.tsv      where the band sits, per lane
    results/tables/12_borderline_per_element.tsv  per-element odds of landing in it
    results/tables/12_borderline_composition.tsv  viral transcript composition by group
    results/tables/12_borderline_host_response.tsv depth-matched host response
    results/tables/12_borderline_isg_sensitivity.tsv what it does to the ISG result

The question
------------
Viral load in this screen is bimodal, but the two modes do not account for every
cell: 2,696 cells (9.7%) sit between them, at viral fractions of 0.6-1.4% where
the bystander mode tops out at 0.5% and the infected mode starts above 2%. A
binary classifier has to put them somewhere, and Sunshine et al. note that
several proviral knockdowns deplete two *bystander* clusters, which they read as
possibly "very early infection, defective particles, or suppressed viral
transcription" without resolving which.

The band is worth resolving because the bystander compartment is where the
interferon response lives. If the band is ambient RNA, it belongs with the
bystanders and nothing changes. If it is early or abortive infection, then
assigning it to bystanders mixes infected cells into the population used to
measure interferon, which is the measurement that matters most in any model
where only a minority of cells get infected - a differentiated airway culture or
a lung organoid more so than a Calu-3 monolayer.

Four tests, each one a prediction only one explanation can satisfy.

1. **Genetics.** Ambient RNA is cell-extrinsic: how much soup a droplet picks up
   has nothing to do with the guide inside it. Infection is cell-intrinsic, so
   knocking down an entry factor must lower it. Estimator and null are the ones
   used for the main hit calls: lane-stratified Mantel-Haenszel against the
   pooled non-targeting cells, with a lane-matched control-resampling null.
2. **Composition.** Ambient soup is the average lysate of infected cells and
   should carry their mix of viral transcripts. Early infection should not:
   genomic ORF1ab arrives first and the subgenomic mRNAs accumulate later.
3. **Host response.** A cell at the start of infection should have begun to
   respond - host shutoff and the NF-kB/immediate-early genes the paper reports
   in infected cells - compared with depth-matched bystanders.
4. **Consequence.** The published analysis has no borderline category, so these
   cells are scored as bystanders there. Recomputing the bystander interferon
   result with and without them says whether that choice matters.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

import _common as C
from _hitcalling_core import lane_counts, mh_log_or, mh_null_draws, two_sided_p

B_DRAWS = 2000
MIN_CELLS = 55
RNG = np.random.default_rng(0)
VIRAL_COLS = ["vc_SCV_ORF1ab", "vc_SCV_S", "vc_SCV_ORF3a", "vc_SCV_E", "vc_SCV_M",
              "vc_SCV_ORF6", "vc_SCV_ORF7ab", "vc_SCV_ORF8", "vc_SCV_N"]
# the infected-cell signature the paper reports: NF-kB and immediate-early stress
STRESS = ["NFKBIA", "ATF3", "FOS", "JUN", "JUNB", "EGR1", "IER2", "KLF6"]
ENTRY = ["ACE2", "TMPRSS2", "ACE2+TMPRSS2", "TMPRSS2+CTSB", "TMPRSS2+CTSL",
         "TMPRSS2+FURIN", "ACE2+FURIN"]


def main() -> None:
    a = C.load(assigned_only=False)
    state = C.read_tbl("03_cell_state.tsv.gz").set_index("cell_barcode")

    obs = a.obs.copy()
    obs["gem_group"] = obs.gem_group.astype(str)
    obs["state"] = state.reindex(obs.index)["state"].to_numpy()
    assigned = obs[obs.element != "unassigned"].copy()
    lanes = sorted(assigned.gem_group.unique())

    tab = (
        assigned.groupby(["gem_group", "state"], observed=True).size()
        .unstack(fill_value=0)
    )
    tab["frac_borderline"] = tab.borderline / tab.sum(axis=1)
    tab["frac_infected"] = tab.infected / tab.sum(axis=1)
    for s in ("uninfected", "borderline", "infected"):
        q = assigned.loc[assigned.state == s].groupby("gem_group").viral_frac
        tab[f"{s}_median_viral_frac"] = q.median()
    tab.reset_index().to_csv(C.tbl("12_borderline_by_lane.tsv"), sep="\t", index=False)
    print("the band, by lane:")
    print(tab.round(4).to_string())

    # ------------------------------------------- 1. does genotype change it?
    pool = assigned[assigned.state.isin(["uninfected", "borderline"])].copy()
    pool["cstate"] = np.where(pool.state == "borderline", "infected", "uninfected")
    ntc = pool[pool.is_ntc]
    c_l, d_l = lane_counts(ntc, lanes, state_col="cstate")
    p_l = c_l / np.maximum(c_l + d_l, 1)

    rows = []
    for element, sub in pool.groupby("element", observed=True):
        if len(sub) < MIN_CELLS:
            continue
        a_l, b_l = lane_counts(sub, lanes, state_col="cstate")
        log_or, se = mh_log_or(a_l, b_l, c_l, d_l)
        if not np.isfinite(log_or):
            continue
        null = mh_null_draws(a_l + b_l, p_l, c_l, d_l, B_DRAWS, RNG)
        rows.append(dict(
            element=element, target=sub.target.iloc[0], is_ntc=bool(sub.is_ntc.iloc[0]),
            n_cells=int(len(sub)), n_borderline=int(a_l.sum()),
            log_or=log_or, se=se,
            p_null=two_sided_p(log_or, null - np.nanmean(null)),
        ))
    pe = pd.DataFrame(rows)
    tgt_mask = (~pe.is_ntc).to_numpy()
    pe["q_null"] = np.nan
    pe.loc[tgt_mask, "q_null"] = C.bh_fdr(pe.loc[tgt_mask, "p_null"].to_numpy())
    pe = pe.sort_values("log_or")

    A = C.read_tbl("07_hits_infection.tsv")[["element", "log_or_mh", "q_null"]].rename(
        columns={"log_or_mh": "log_or_infection", "q_null": "q_infection"})
    pe = pe.merge(A, on="element", how="left")
    pe.to_csv(C.tbl("12_borderline_per_element.tsv"), sep="\t", index=False)

    t = pe[~pe.is_ntc]
    sig = t[t.q_null < 0.1]
    print(f"\n1. GENETICS: {len(t)} targeting elements tested, {len(sig)} change the odds "
          f"of a cell landing in the band at q < 0.1")
    print(t.nsmallest(10, "log_or")[
        ["target", "n_cells", "n_borderline", "log_or", "q_null", "log_or_infection"]
    ].round(3).to_string(index=False))
    n_ntc = pe[pe.is_ntc]
    print(f"   controls at p < 0.1: {int((n_ntc.p_null < 0.1).sum())} of {len(n_ntc)}")
    m = t.dropna(subset=["log_or_infection"])
    print(f"   agreement with the established-infection effect: "
          f"r = {np.corrcoef(m.log_or, m.log_or_infection)[0,1]:.3f}, "
          f"Spearman = {m.log_or.corr(m.log_or_infection, method='spearman'):.3f} "
          f"over {len(m)} elements")
    ent = t[t.target.isin(ENTRY)]
    print("   entry factors specifically:")
    print(ent[["target", "n_cells", "log_or", "q_null", "log_or_infection"]]
          .round(3).to_string(index=False))

    # ------------------------------------------------------- 2. composition
    groups = {
        "infected": assigned[assigned.state == "infected"],
        "borderline": assigned[assigned.state == "borderline"],
        "bystander": assigned[assigned.state == "uninfected"],
    }
    comp_rows = []
    for name, g in groups.items():
        tot = float(g[VIRAL_COLS].to_numpy().sum())
        share = g[VIRAL_COLS].to_numpy().sum(axis=0) / tot
        comp_rows.append(dict(group=name, n_cells=len(g), viral_umi=tot,
                              **{c.replace("vc_SCV_", ""): float(s)
                                 for c, s in zip(VIRAL_COLS, share)}))
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(C.tbl("12_borderline_composition.tsv"), sep="\t", index=False)
    print("\n2. COMPOSITION: share of viral UMIs carried by each transcript")
    print(comp.round(4).to_string(index=False))
    ratio = (comp.set_index("group").ORF1ab / comp.set_index("group").N)
    print("   ORF1ab : N ratio  " + ", ".join(f"{k} {v:.3f}" for k, v in ratio.items()))

    # ------------------------------------------------------ 3. host response
    host = (a.var.kind == "host").to_numpy()
    h = a[:, host].copy()
    sc.pp.normalize_total(h, target_sum=1e6)
    sc.pp.log1p(h)
    sc.tl.score_genes(h, [g for g in C.ISG_GENES if g in h.var_names],
                      score_name="isg", random_state=0)
    sc.tl.score_genes(h, [g for g in STRESS if g in h.var_names],
                      score_name="stress", random_state=0)
    sco = pd.DataFrame({"isg": h.obs.isg.to_numpy(), "stress": h.obs.stress.to_numpy()},
                       index=h.obs_names)
    assigned = assigned.join(sco)

    assigned["depth_bin"] = pd.qcut(assigned.total_umi, 10, labels=False, duplicates="drop")
    hr = []
    for bin_id, sub in assigned.groupby("depth_bin"):
        g = {s: sub[sub.state == s] for s in ("uninfected", "borderline", "infected")}
        if min(len(v) for v in g.values()) < 15:
            continue
        hr.append(dict(
            depth_bin=int(bin_id), median_total_umi=float(sub.total_umi.median()),
            n_bystander=len(g["uninfected"]), n_borderline=len(g["borderline"]),
            n_infected=len(g["infected"]),
            host_umi_bystander=float(g["uninfected"].host_umi.median()),
            host_umi_borderline=float(g["borderline"].host_umi.median()),
            host_umi_infected=float(g["infected"].host_umi.median()),
            stress_bystander=float(g["uninfected"].stress.mean()),
            stress_borderline=float(g["borderline"].stress.mean()),
            stress_infected=float(g["infected"].stress.mean()),
            isg_bystander=float(g["uninfected"].isg.mean()),
            isg_borderline=float(g["borderline"].isg.mean()),
            isg_infected=float(g["infected"].isg.mean()),
            # doublet proxy: at matched library size, a doublet detects more genes
            n_genes_bystander=float(g["uninfected"].n_genes.median()),
            n_genes_borderline=float(g["borderline"].n_genes.median()),
            n_genes_infected=float(g["infected"].n_genes.median()),
        ))
    hrt = pd.DataFrame(hr)
    hrt.to_csv(C.tbl("12_borderline_host_response.tsv"), sep="\t", index=False)
    print("\n3. HOST RESPONSE, within deciles of total library size")
    print(hrt[["depth_bin", "median_total_umi", "n_borderline",
               "host_umi_bystander", "host_umi_borderline", "host_umi_infected",
               "stress_bystander", "stress_borderline", "stress_infected"]]
          .round(3).to_string(index=False))
    print(hrt[["depth_bin", "isg_bystander", "isg_borderline", "isg_infected"]]
          .round(3).to_string(index=False))
    print("   doublet check - genes detected at matched library size:")
    print(hrt[["depth_bin", "n_genes_bystander", "n_genes_borderline", "n_genes_infected"]]
          .round(1).to_string(index=False))

    # -------------------------------------- 4. consequence for the ISG result
    isg_rows = []
    pops = {
        "bystander_only": assigned[assigned.state == "uninfected"],
        "bystander_plus_borderline": assigned[assigned.state.isin(["uninfected", "borderline"])],
    }
    for label, pop in pops.items():
        ntc_p = pop[pop.is_ntc]
        poolv = {l: ntc_p.loc[ntc_p.gem_group == l, "isg"].to_numpy() for l in lanes}
        base = {l: float(np.mean(v)) if len(v) else np.nan for l, v in poolv.items()}
        for element, sub in pop.groupby("element", observed=True):
            if len(sub) < 30:
                continue
            n_l, num, w = [], 0.0, 0.0
            for l in lanes:
                s = sub.loc[sub.gem_group == l, "isg"].to_numpy()
                if len(s) == 0 or not np.isfinite(base[l]):
                    n_l.append(0)
                    continue
                n_l.append(len(s)); num += len(s) * (s.mean() - base[l]); w += len(s)
            if w == 0:
                continue
            obs_stat = num / w
            null = np.empty(B_DRAWS)
            for i in range(B_DRAWS):
                nn = 0.0
                for l, k in zip(lanes, n_l):
                    if k == 0:
                        continue
                    p = poolv[l]
                    nn += k * (p[RNG.integers(0, len(p), k)].mean() - base[l])
                null[i] = nn / w
            isg_rows.append(dict(
                population=label, element=element, target=sub.target.iloc[0],
                is_ntc=bool(sub.is_ntc.iloc[0]), n_cells=int(len(sub)),
                delta_isg=obs_stat,
                p_null=(1 + int(np.sum(np.abs(null) >= abs(obs_stat)))) / (B_DRAWS + 1)))
    si = pd.DataFrame(isg_rows)
    si["q_null"] = np.nan
    for pop_name, idx in si[~si.is_ntc].groupby("population").groups.items():
        si.loc[idx, "q_null"] = C.bh_fdr(si.loc[idx, "p_null"].to_numpy())
    si.to_csv(C.tbl("12_borderline_isg_sensitivity.tsv"), sep="\t", index=False)

    w = si[~si.is_ntc].pivot(index="element", columns="population",
                             values=["delta_isg", "q_null"]).dropna()
    r = float(np.corrcoef(w[("delta_isg", "bystander_only")],
                          w[("delta_isg", "bystander_plus_borderline")])[0, 1])
    hits_clean = set(si[(si.population == "bystander_only") & (~si.is_ntc)
                        & (si.q_null < 0.1)].element)
    hits_mixed = set(si[(si.population == "bystander_plus_borderline") & (~si.is_ntc)
                        & (si.q_null < 0.1)].element)
    print("\n4. CONSEQUENCE for the bystander interferon result")
    print(f"   effect sizes, clean bystanders vs bystanders+band: r = {r:.4f}")
    print(f"   hits at q < 0.1: {len(hits_clean)} clean, {len(hits_mixed)} with the band "
          f"folded in, {len(hits_clean & hits_mixed)} shared")
    print(f"   only without the band: {sorted(hits_clean - hits_mixed)}")
    print(f"   only with the band:    {sorted(hits_mixed - hits_clean)}")
    for tgt in ["GPR89A", "GPR89A+GPR89B", "KEAP1", "SPNS1", "BRD2", "STAT2", "IRF9"]:
        sub = si[si.target == tgt]
        if len(sub):
            vals = "; ".join(f"{r_.population}: delta {r_.delta_isg:+.3f} q {r_.q_null:.3f}"
                             for _, r_ in sub.iterrows())
            print(f"   {tgt}  {vals}")


if __name__ == "__main__":
    main()
