# CRISPRi Perturb-seq: QC and hit calling on a SARS-CoV-2 host-factor screen

[![tests](https://github.com/hbusra/crispri-perturbseq-sarscov2/actions/workflows/tests.yml/badge.svg)](https://github.com/hbusra/crispri-perturbseq-sarscov2/actions/workflows/tests.yml)

I reanalysed the Calu-3 CRISPRi Perturb-seq screen of SARS-CoV-2 host factors
from [Sunshine et al., *Nature Communications* 14 (2023)](https://doi.org/10.1038/s41467-023-41788-4),
starting from the processed count matrix in [GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240)
and writing the QC and hit-calling code from scratch. Two reasons: to check
whether the published calls hold up under a different statistical approach, and
to build the per-element QC I would want on any pooled screen before trusting
its hit list.

27,882 cells with a guide assignment, 237 library elements (22 non-targeting),
eight 10x channels, one infection experiment.

## What I found

**1. The published viral-load calls reproduce.** All 21 calls in the paper's
Fig. 3A come back in the right direction: the 14 reported hits at q < 0.05 and
the 7 reported non-hits failing q < 0.1. The interferon panel (paper Fig. 5A)
reproduces too — all four novel hits at q < 0.01 (KEAP1, SPNS1, BRD2, and
GPR89A/B through the GPR89A and GPR89A+GPR89B elements) along with the JAK-STAT
axis. One element-level exception: GPR89B on its own does nothing here (rank 131
of 190), while its paralogue GPR89A is rank 9. The paper reports the pair
together.

![Figure 1](results/figures/fig1_concordance.png)

**2. The FURIN element never knocks FURIN down.** Across the three elements
carrying FURIN, target mRNA changes by +6% to −19% against a library median
knockdown of 88%. The paper uses FURIN-alone as a negative result in its
protease-redundancy argument; on this data that arm is an uninterpretable
reagent rather than a negative. MEPCE and LARP1 have the same problem. This is
the first thing I check on a screen and the reason I check it.

**3. Internal reproducibility is weaker than a pooled number suggests.** Twelve
genes in this library are carried by two separately cloned single-guide vectors.
Those 12 pairs agree at r = 0.55 on knockdown. Pooling them with pairs that
share guide design gives r = 0.74, which measures the same reagent twice as much
as it measures reproducibility, so I report the 12.

## How hits are called

Three questions had to be answered before I would use any hit list from this
screen:

| | Question | Script | Answer |
|---|---|---|---|
| 1 | Did the perturbation work? | `04_guide_qc.py` | median 88% on-target knockdown; 76% of element–gene pairs clear 75%; FURIN, MEPCE and LARP1 do not knock down |
| 2 | Is the phenotype real or a depth artefact? | `03_infection.py` | three threshold-free classifiers agree ≥ 99.3%; calls stable to 10× downsampling |
| 3 | Is the effect bigger than control-to-control noise? | `07_hit_calling.py` | lane-matched resampling null from the non-targeting cells; controls checked per phenotype (infection 0/20 at p < 0.1, burden 3/21, ISG 3/20, against 10% expected) |

The analysis choices that matter most:

- **Infection state** comes from a two-component mixture on log10 viral fraction
  (29.2% infected, 61.1% bystander, 9.7% borderline and dropped). Three other
  rules are run alongside it, including the paper's own.
- **Two phenotypes, kept separate.** Whether a cell gets infected at all is a
  Mantel–Haenszel odds ratio stratified by 10x channel (30 of 164 elements at
  q < 0.1). How much virus an infected cell carries is a negative-binomial GLM
  with `log(host UMI)` as the offset (5 of 149). The offset is host UMI, never
  total UMI: total UMI contains the viral counts being modelled and is itself
  pushed down by shutoff, so using it means conditioning on the outcome.
- **Lane is a blocking factor, not noise.** Channels 5–8 carry 1.4× the
  infection rate of channels 1–4, matching the paper's pulsed and non-pulsed
  conditions. Fitting lane as a random intercept instead
  (`11_mixed_effects.py`) gives the same answers (r = 0.9999 on effect sizes,
  27 of 30 hits shared), so the stratified estimator stays primary.
- **Significance** comes from resampling the control cells: 2,000
  lane-matched pseudo-elements per element, drawn from the pooled non-targeting
  pool. It is not a permutation of perturbation labels, and
  [docs/statistical_methods.md](docs/statistical_methods.md) spells out what it
  does instead.
- **Knockdown** is measured in bystander cells only. Infected cells have lost
  most of their host mRNA to NSP1 shutoff, which makes every knockdown ratio
  look better than it is.

Everything else — all seven figures, the discrepancies with the paper, the
sensitivity analyses — is in [docs/results.md](docs/results.md).

## What this screen cannot answer

The eight `gem_group` values are 10x channels from one infection experiment, not
eight infections. The resampling null covers cell- and channel-level sampling
noise and nothing above that, so none of these p-values speak to whether the
screen would replicate in a second infection. Separately, 47 of 237 elements
have fewer than 55 cells, mostly essential genes whose knockdown kills cells;
those targets are missing from the hit list for a systematic reason, not because
they were tested and found null.

## Repository layout

```
scripts/
  00_download.sh              GEO download, checksums, unpacking
  01_build_anndata.py         AnnData with host/viral/lenti split, guide calls, lanes
  02_qc.py                    cell QC, lane structure, element representation
  03_infection.py             infection classification (4 rules) + depth sensitivity
  04_guide_qc.py              on-target knockdown, construct concordance
  05_viral_load.py            descriptive reproduction of the paper's Fig. 3A
  06_isg.py                   reproduction of the paper's Fig. 5A
  07_hit_calling.py           MH odds ratios, NB GLM, control resampling null, FDR
  08_representation.py        subsampling curves, resolution, underpowered elements
  09_infected_vs_bystander.py do perturbation effects survive infection?
  10_figures.py               the figures
  11_mixed_effects.py         lane random-intercept GLMM vs the stratified estimator
  _common.py                  paths, gene sets, BH-FDR
  _hitcalling_core.py         estimators and nulls shared by 07 and 08
tests/
  test_estimators.py          12 tests: MH vs statsmodels, NB recovery, null calibration
docs/
  results.md                  all seven figures, the discrepancies, the limits
  statistical_methods.md      estimands, nulls, testing families, sensitivity checks
  RUNBOOK.md                  every command in order, with the output to expect
results/tables/               every number quoted anywhere here, as .tsv
results/figures/              the figures
```

`10_figures.py` reads only the tables in `results/tables/`, so a figure can be
redrawn without rerunning any statistics.

## Running it

```bash
conda env create -f environment.yml && conda activate perturbseq-cov2
bash scripts/00_download.sh data/raw            # ~1.2 GB from GEO
python scripts/01_build_anndata.py --raw data/raw --out data
for s in 02_qc 03_infection 04_guide_qc 05_viral_load 06_isg \
         07_hit_calling 08_representation 09_infected_vs_bystander 10_figures \
         11_mixed_effects; do
    python scripts/$s.py
done
python -m pytest tests -q                       # 12 estimator tests, ~10 s
```

About 14 minutes on 4 cores, 11 of which are the one-off parse of the
337-million-entry matrix in `01`. That result is cached, so reruns take under 3
minutes. Peak memory ~6 GB, also in `01`.
[docs/RUNBOOK.md](docs/RUNBOOK.md) lists each step with the numbers it should
print.

## Data

Everything comes from GEO [GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240),
deposited by the authors of the original study; nothing in `data/` is committed.
I start from their per-cell guide assignments (`cell_identities.csv`), so this is
a reanalysis from the processed matrix, not from FASTQs.

If you use the data, cite the original study:

> Sunshine S, Puschnik AS, Replogle JM, Laurie MT, Liu J, Zha BS, Nuñez JK,
> Byrum JR, McMorrow AH, Frieman MB, Winkler J, Qiu X, Rosenberg OS, Leonetti MD,
> Ye CJ, Weissman JS, DeRisi JL, Hein MY.
> Systematic functional interrogation of SARS-CoV-2 host factors using Perturb-seq.
> *Nature Communications* 14 (2023). doi:10.1038/s41467-023-41788-4

## Licence

Code here is MIT licensed (see [LICENSE](LICENSE)). That covers the code only:
the GSE208240 data and the original figures stay under their own terms, and the
results being reproduced belong to their authors.
