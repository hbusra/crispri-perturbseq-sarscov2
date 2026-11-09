# A guide-aware QC and hit-calling workflow for pooled single-cell perturbation screens

An independently written pipeline for QC'ing and calling hits in a CRISPRi
Perturb-seq screen, validated by reproducing a published screen of SARS-CoV-2
host factors in Calu-3 cells ([Sunshine et al., *Nature Communications* 14, 2023](https://doi.org/10.1038/s41467-023-41788-4);
data [GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240)).

**Three headline findings**

1. **The screen reproduces.** All 21 of the published viral-load calls come back
   in the right direction — the 14 reported hits at q < 0.05 and the 7 reported
   non-hits failing q < 0.1 — using a different statistical approach from the
   original. All four of the paper's novel interferon hits (GPR89A/B, KEAP1,
   SPNS1, BRD2) and the canonical JAK-STAT axis return at q < 0.01.
2. **One published arm rests on a reagent that does not work.** The FURIN
   elements show **no knockdown of FURIN** (+6% to −19%, against a library median
   of 88%), so the single-gene FURIN arm of the protease-redundancy argument is
   uninterpretable rather than negative. The same holds for MEPCE and LARP1.
3. **The screen's internal reproducibility is more modest than it looks.** The 12
   genes carried by two separately cloned vectors agree at r = 0.55 on knockdown.
   Pooling them with pairs that share guide design inflates this to r = 0.74,
   which is the kind of number worth not quoting.

![Figure 4](results/figures/fig4_concordance.png)

Bar colour is what the *paper* claims; fill is what *we* call. Agreeing that a
target has **no** effect is a reproduction too, so the grey bars count.

---

## The question

> How reliably can we identify perturbations that change an infection phenotype
> in a pooled single-cell screen, and how do we know the hits are real?

Three things have to be true before a hit list means anything, and each gets its
own script:

| | Question | Script | Answer here |
|---|---|---|---|
| 1 | Did the perturbation work? | `04_guide_qc.py` | median 88% on-target knockdown; 76% of element–gene pairs clear 75%; **FURIN does not knock down at all** |
| 2 | Is the phenotype real, or a measurement artefact? | `03_infection.py` | three independent classifiers agree ≥ 99.3%; calls stable to 10× downsampling |
| 3 | Is the effect bigger than control-to-control noise? | `07_hit_calling.py` | lane-matched resampling null from the non-targeting controls; 0 of 20 control elements reach p < 0.1 |

## Methods in brief

- **Infection state** from a two-component mixture on log10 viral fraction
  (29.2% infected, 61.1% bystander, 9.7% borderline and excluded), cross-checked
  against three other rules including the paper's own.
- **Two estimands, kept apart.** Infection *establishment* is a Mantel–Haenszel
  odds ratio stratified by 10x channel (30 of 164 elements at q < 0.1). Viral
  *burden given infection* is a negative-binomial GLM with `log(host UMI)` as an
  offset (5 of 149). The offset is never total UMI, which contains the outcome.
- **Blocking.** 10x channels 5–8 carry 1.4× the infection rate of channels 1–4,
  so lane is a stratum in every model and every null.
- **Significance** from a control-based resampling null: 2,000 lane-matched
  pseudo-elements drawn from the pooled non-targeting cells. This is not a
  permutation of perturbation labels, and
  [docs/statistical_methods.md](docs/statistical_methods.md) says exactly what
  it is.
- **Knockdown** measured in bystander cells only, pseudobulk, because shutoff in
  infected cells flatters every ratio.

**Full results, figures and discrepancies:** [docs/results.md](docs/results.md).
**Estimands, nulls, testing families and sensitivity analyses:** [docs/statistical_methods.md](docs/statistical_methods.md).

## What this design cannot claim

The eight `gem_group` values are 10x channels from a **single infection
experiment**, not eight independent infections. The resampling null controls for
cell-level and channel-level sampling noise; it does **not** establish
reproducibility across independent biological experiments, and no p-value here
should be read that way. 47 of 237 elements fall below the 55-cell cut, mostly
essential genes whose knockdown kills cells, so their absence from the hit list
is a systematic gap rather than evidence of no effect.

---

## Repository layout

```
scripts/
  00_download.sh              GEO download, checksums, unpacking
  01_build_anndata.py         AnnData with host/viral/lenti split, guide calls, lanes
  02_qc.py                    cell QC, lane structure, element representation
  03_infection.py             infection classification (4 rules) + depth sensitivity
  04_guide_qc.py              on-target knockdown, construct concordance
  05_viral_load.py            descriptive reproduction of the published Fig. 3A
  06_isg.py                   reproduction of the published Fig. 5
  07_hit_calling.py           MH odds ratios, NB GLM, control resampling null, FDR
  08_representation.py        subsampling curves, resolution, underpowered elements
  09_infected_vs_bystander.py the extension
  10_figures.py               every figure, drawn from the tables only
  _common.py                  paths, gene sets, BH-FDR
  _hitcalling_core.py         the estimators and nulls shared by 07 and 08
tests/
  test_estimators.py          12 tests: MH vs statsmodels, NB recovery, null calibration
docs/
  results.md                  every figure, the discrepancies, the limitations
  statistical_methods.md      estimands, nulls, and what each one licenses
  RUNBOOK.md                  every command, in order, with expected output
results/tables/               every number in this repository, as .tsv
results/figures/              the figures
```

Analysis and plotting are kept apart: `10_figures.py` reads only from
`results/tables/`, so any figure can be regenerated or redrawn without rerunning
a single statistic.

## Running it

```bash
conda env create -f environment.yml && conda activate perturbseq-cov2
bash scripts/00_download.sh data/raw            # ~1.2 GB from GEO
python scripts/01_build_anndata.py --raw data/raw --out data
for s in 02_qc 03_infection 04_guide_qc 05_viral_load 06_isg \
         07_hit_calling 08_representation 09_infected_vs_bystander 10_figures; do
    python scripts/$s.py
done
python -m pytest tests -q                       # 12 estimator tests, ~10 s
```

About 14 minutes end to end on 4 cores, of which 11 are the one-off parse of the
337-million-entry matrix in `01`; that result is cached, so every rerun of the
analysis is under 3 minutes. Peak memory is about 6 GB, also in `01`. See
[docs/RUNBOOK.md](docs/RUNBOOK.md) for the step-by-step version with the expected
output at each stage.

## Data

Raw and processed data are from GEO [GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240),
deposited by the authors of the original study. Nothing in `data/` is committed.
The authors' own per-cell guide assignments (`cell_identities.csv`) are used as
the starting point, so this is a reanalysis from the processed matrix rather than
from FASTQs.

If you use the data, cite the original study:

> Sunshine S, Puschnik AS, Replogle JM, Laurie MT, Liu J, Zha BS, Nuñez JK,
> Byrum JR, McMorrow AH, Frieman MB, Winkler J, Qiu X, Rosenberg OS, Leonetti MD,
> Ye CJ, Weissman JS, DeRisi JL, Hein MY.
> Systematic functional interrogation of SARS-CoV-2 host factors using Perturb-seq.
> *Nature Communications* 14 (2023). doi:10.1038/s41467-023-41788-4

## Licence

Code in this repository is MIT licensed (see [LICENSE](LICENSE)). That covers the
code only: the GSE208240 data and the original publication's figures remain under
their own terms, and the reproduction above is of published results that belong
to their authors.
