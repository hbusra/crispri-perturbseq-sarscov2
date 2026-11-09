#!/usr/bin/env python
"""
06 - Reproduction of the interferon panel (paper Fig. 5): bystander ISG scores.

Writes
    results/tables/06_isg_by_element.tsv
    results/tables/06_isg_hits.tsv
    results/tables/06_isg_paper_comparison.tsv

The published analysis restricts the ISG score to bystander cells, because
infected cells actively suppress interferon signalling, so a score computed on
a mixed population mostly measures each element's infected fraction rather than
its effect on interferon. We keep that restriction and show what happens
without it.

Two traps, both handled explicitly:

1. Score self-containment. Several ISGs are themselves CRISPRi targets in this
   library (STAT1, STAT2, IRF9, BST2, among others). Knocking a gene down
   mechanically lowers any score that contains it, which looks like a signalling
   phenotype but is just the on-target effect. The score is therefore computed
   both with and without the element's own target genes, and the table carries
   both.

2. Depth. The score is a difference of means of log1p host CPM against an
   expression-matched background, so it is depth-normalised by construction,
   but bystander depth still varies by element. The per-element test keeps lane
   as a stratum and the null is lane- and size-matched, as in 07.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

import _common as C

B_DRAWS = 2000
MIN_BYSTANDER = 30
RNG = np.random.default_rng(0)

# Elements the paper reports as lowering the bystander ISG score.
PAPER_ISG_DOWN = ["GPR89A", "GPR89B", "GPR89A+GPR89B", "KEAP1", "SPNS1", "BRD2"]
# The interferon-pathway group the paper highlights as a functional cluster.
PAPER_IFN_PATHWAY = ["IFNAR1", "IFNAR2", "STAT1", "STAT2", "IRF3", "IRF9"]


def main() -> None:
    a = C.load(assigned_only=True)
    state = C.read_tbl("03_cell_state.tsv.gz").set_index("cell_barcode")
    a.obs["state"] = state.reindex(a.obs_names)["state"].to_numpy()

    # ------------------------------------------------- host-only normalisation
    host = (a.var.kind == "host").to_numpy()
    b = a[:, host].copy()
    sc.pp.normalize_total(b, target_sum=1e6)   # host-only denominator
    sc.pp.log1p(b)

    present = [g for g in C.ISG_GENES if g in b.var_names]
    print(f"ISG set: {len(present)}/{len(C.ISG_GENES)} genes present")
    sc.tl.score_genes(b, present, score_name="isg", random_state=0)
    a.obs["isg"] = b.obs.isg.to_numpy()

    # ISG score excluding any gene that is a CRISPRi target in this library
    targets_in_library = set()
    for t in a.obs.target.unique():
        if t not in ("NTC", "unassigned"):
            targets_in_library.update(t.split("+"))
    clean = [g for g in present if g not in targets_in_library]
    dropped = sorted(set(present) - set(clean))
    print(f"ISGs that are also CRISPRi targets, dropped from the clean score: {dropped}")
    sc.tl.score_genes(b, clean, score_name="isg_clean", random_state=0)
    a.obs["isg_clean"] = b.obs.isg_clean.to_numpy()

    # the rare interferon-producing subset the paper describes
    ifn_present = [g for g in C.IFN_PRODUCING if g in b.var_names]
    a.obs["ifn_producing_umi"] = np.asarray(
        a[:, ifn_present].X.sum(axis=1)
    ).ravel() if ifn_present else 0.0

    d = a.obs.copy()
    d["gem_group"] = d.gem_group.astype(str)
    byst = d[d.state == "uninfected"]
    inf = d[d.state == "infected"]
    print(
        f"\nISG score: bystander mean {byst.isg.mean():.4f}, infected mean {inf.isg.mean():.4f} "
        f"(the paper reports infected cells lower, i.e. active suppression)"
    )
    print(
        f"interferon-producing cells (>=5 UMI of {len(ifn_present)} IFN/chemokine genes): "
        f"bystander {(byst.ifn_producing_umi >= 5).mean() * 100:.2f}%, "
        f"infected {(inf.ifn_producing_umi >= 5).mean() * 100:.2f}%"
    )

    # ------------------------------------------- per-element test in bystanders
    lanes = sorted(byst.gem_group.unique())
    rows = []
    for score_col in ("isg", "isg_clean"):
        ntc = byst[byst.is_ntc]
        pool = {l: ntc.loc[ntc.gem_group == l, score_col].to_numpy() for l in lanes}
        base = {l: float(np.mean(v)) if len(v) else np.nan for l, v in pool.items()}

        for element, sub in byst.groupby("element", observed=True):
            if len(sub) < MIN_BYSTANDER:
                continue
            n_l, num, w = [], 0.0, 0.0
            for l in lanes:
                s = sub.loc[sub.gem_group == l, score_col].to_numpy()
                if len(s) == 0 or not np.isfinite(base[l]):
                    n_l.append(0)
                    continue
                n_l.append(len(s))
                num += len(s) * (s.mean() - base[l])
                w += len(s)
            if w == 0:
                continue
            obs = num / w
            null = np.empty(B_DRAWS)
            for i in range(B_DRAWS):
                nn = 0.0
                for l, k in zip(lanes, n_l):
                    if k == 0:
                        continue
                    p = pool[l]
                    nn += k * (p[RNG.integers(0, len(p), k)].mean() - base[l])
                null[i] = nn / w
            sd = float(np.std(null))
            rows.append(
                dict(
                    score=score_col,
                    element=element,
                    target=sub.target.iloc[0],
                    is_ntc=bool(sub.is_ntc.iloc[0]),
                    n_bystander=int(len(sub)),
                    delta_isg=obs,
                    null_sd=sd,
                    z_empirical=obs / sd if sd > 0 else np.nan,
                    p_perm=(1 + int(np.sum(np.abs(null) >= abs(obs)))) / (B_DRAWS + 1),
                )
            )
    out = pd.DataFrame(rows)
    out["q_perm"] = np.nan
    for sc_name, idx in out.groupby("score").groups.items():
        out.loc[idx, "q_perm"] = C.bh_fdr(out.loc[idx, "p_perm"].to_numpy())
    out = out.sort_values(["score", "delta_isg"])
    out.to_csv(C.tbl("06_isg_by_element.tsv"), sep="\t", index=False)

    main_tab = out[out.score == "isg_clean"]
    tested = main_tab[~main_tab.is_ntc]
    sig = tested[tested.q_perm < 0.1]
    print(
        f"\nbystander ISG, target-free score: {len(tested)} elements tested, "
        f"{len(sig)} at q<0.1 ({int((sig.delta_isg < 0).sum())} down, "
        f"{int((sig.delta_isg > 0).sum())} up)"
    )
    print("\nlargest decreases:")
    print(
        tested.nsmallest(15, "delta_isg")[
            ["element", "target", "n_bystander", "delta_isg", "z_empirical", "q_perm"]
        ].round(4).to_string(index=False)
    )
    print("\nlargest increases:")
    print(
        tested.nlargest(8, "delta_isg")[
            ["element", "target", "n_bystander", "delta_isg", "z_empirical", "q_perm"]
        ].round(4).to_string(index=False)
    )
    sig.to_csv(C.tbl("06_isg_hits.tsv"), sep="\t", index=False)

    # ------------------------------------------------- against the paper's list
    rows = []
    for tgt in PAPER_ISG_DOWN + PAPER_IFN_PATHWAY:
        for score_col in ("isg", "isg_clean"):
            sub = out[(out.target == tgt) & (out.score == score_col)]
            if sub.empty:
                rows.append(dict(target=tgt, score=score_col, note="not testable"))
                continue
            for _, r in sub.iterrows():
                n = out[(out.score == score_col) & (~out.is_ntc)]
                rows.append(
                    dict(
                        target=tgt,
                        score=score_col,
                        element=r.element,
                        paper=("down" if tgt in PAPER_ISG_DOWN else "pathway group"),
                        n_bystander=int(r.n_bystander),
                        delta_isg=r.delta_isg,
                        q_perm=r.q_perm,
                        rank=int((n.delta_isg < r.delta_isg).sum() + 1),
                        n_tested=len(n),
                    )
                )
    cmp = pd.DataFrame(rows)
    cmp.to_csv(C.tbl("06_isg_paper_comparison.tsv"), sep="\t", index=False)
    print("\nthe paper's interferon results in our ranking:")
    print(cmp.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
