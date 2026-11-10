#!/usr/bin/env python
"""
05 - Reproduction of the published viral-load panel (paper Fig. 3A).

Writes
    results/tables/05_viral_load_by_element.tsv
    results/tables/05_paper_hits_descriptive.tsv

This is the descriptive reproduction: viral load per library element against the
non-targeting controls, with no model. It exists so the figure can be compared
to the published one by eye, and so the statistical hit calling in 07 has a
plain-language counterpart.

Two phenotypes are kept apart throughout, because they answer different
questions and can move in opposite directions:

  infected_frac   the fraction of an element's cells classified as infected
                  -> how often does infection get established?
  viral_rate      viral UMIs per host UMI, among infected cells only
                  -> once infection is established, how much virus accumulates?

Reporting a single "viral load" per perturbation mixes the two. A knockdown
that blocks entry lowers the first; a knockdown that limits replication lowers
the second.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import _common as C

# Host factors the paper reports as changing viral RNA, for the concordance
# table. Direction is the published one: "down" = less viral RNA on knockdown.
PAPER_HITS = {
    "ACE2": "down", "TMPRSS2": "down", "ACE2+TMPRSS2": "down",
    "TMPRSS2+FURIN": "down", "TMPRSS2+CTSB": "down", "TMPRSS2+CTSL": "down",
    "BRD2": "down", "ATG14": "down", "EIF4E2": "down", "EIF4H": "down",
    "VMP1": "down", "MPP5": "down", "NFKBIA": "down",
    "SEC62": "up",
}
PAPER_NULL = {"FURIN": "ns", "CTSB": "ns", "CTSL": "ns", "RELA": "ns",
              "RELB": "ns", "RELA+RELB": "ns", "EIF4B": "ns"}


def main() -> None:
    d = C.read_tbl("03_cell_state.tsv.gz")
    d["gem_group"] = d.gem_group.astype(str)
    d["viral_rate"] = d.viral_umi / d.host_umi.replace(0, np.nan)

    ntc = d[d.is_ntc]
    print(
        f"cells: {len(d)}; NTC cells: {len(ntc)} across "
        f"{d.loc[d.is_ntc, 'element'].nunique()} non-targeting elements"
    )

    inf = d[d.state == "infected"]

    g = d.groupby("element", observed=True)
    tab = pd.DataFrame(
        {
            "target": g.target.first(),
            "is_ntc": g.is_ntc.first(),
            "n_cells": g.size(),
            "n_lanes": g.gem_group.nunique(),
            "n_infected": g.state.apply(lambda s: int((s == "infected").sum())),
            "n_bystander": g.state.apply(lambda s: int((s == "uninfected").sum())),
            "host_umi_median": g.host_umi.median(),
        }
    )
    tab["infected_frac"] = tab.n_infected / (tab.n_infected + tab.n_bystander)

    gi = inf.groupby("element", observed=True)
    tab["viral_rate_median_infected"] = gi.viral_rate.median()
    tab["log10_viral_rate_mean_infected"] = gi.viral_rate.apply(
        lambda s: float(np.mean(np.log10(s.clip(lower=1e-6))))
    )
    tab["viral_frac_median_infected"] = gi.viral_frac.median()

    # NTC baseline, pooled over the 22 non-targeting elements
    ntc_inf = inf[inf.is_ntc]
    base_if = len(ntc_inf) / len(ntc[ntc.state.isin(["infected", "uninfected"])])
    base_rate = float(np.mean(np.log10(ntc_inf.viral_rate.clip(lower=1e-6))))
    print(f"\nNTC baseline: infected fraction {base_if:.3f}, "
          f"mean log10 viral rate among infected {base_rate:.3f}")

    tab["infected_frac_vs_ntc"] = tab.infected_frac - base_if
    tab["log2_odds_vs_ntc"] = np.log2(
        (tab.infected_frac / (1 - tab.infected_frac)) / (base_if / (1 - base_if))
    )
    tab["delta_log10_viral_rate"] = tab.log10_viral_rate_mean_infected - base_rate

    tab = tab.reset_index().sort_values("log2_odds_vs_ntc")
    tab.to_csv(C.tbl("05_viral_load_by_element.tsv"), sep="\t", index=False)

    usable = tab[(tab.n_cells >= 55) & ~tab.is_ntc]
    print(f"\nelements with >=55 cells (the paper's coverage cut): {len(usable)} targeting + "
          f"{int(((tab.n_cells >= 55) & tab.is_ntc).sum())} NTC")

    print("\nlargest decrease in infected fraction vs NTC:")
    print(
        usable.nsmallest(15, "log2_odds_vs_ntc")[
            ["element", "target", "n_cells", "infected_frac", "log2_odds_vs_ntc",
             "delta_log10_viral_rate"]
        ].round(3).to_string(index=False)
    )
    print("\nlargest increase in infected fraction vs NTC:")
    print(
        usable.nlargest(10, "log2_odds_vs_ntc")[
            ["element", "target", "n_cells", "infected_frac", "log2_odds_vs_ntc",
             "delta_log10_viral_rate"]
        ].round(3).to_string(index=False)
    )

    # ------------------------------------------- where the paper's hits land
    rows = []
    for tgt, direction in {**PAPER_HITS, **PAPER_NULL}.items():
        sub = tab[tab.target == tgt]
        if sub.empty:
            rows.append(dict(target=tgt, paper=direction, note="target not in library"))
            continue
        for _, r in sub.iterrows():
            rows.append(
                dict(
                    target=tgt,
                    paper=direction,
                    element=r.element,
                    n_cells=int(r.n_cells),
                    infected_frac=r.infected_frac,
                    log2_odds_vs_ntc=r.log2_odds_vs_ntc,
                    delta_log10_viral_rate=r.delta_log10_viral_rate,
                    rank_of_237=int((tab.log2_odds_vs_ntc < r.log2_odds_vs_ntc).sum() + 1),
                )
            )
    ph = pd.DataFrame(rows)
    ph.to_csv(C.tbl("05_paper_hits_descriptive.tsv"), sep="\t", index=False)
    print("\nthe paper's reported hits, ranked by the infected-fraction effect here:")
    print(ph.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
