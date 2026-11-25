# crispri-perturbseq-sarscov2

[![tests](https://github.com/hbusra/crispri-perturbseq-sarscov2/actions/workflows/tests.yml/badge.svg)](https://github.com/hbusra/crispri-perturbseq-sarscov2/actions/workflows/tests.yml)

Reanalysis of the Calu-3 CRISPRi Perturb-seq screen of SARS-CoV-2 host factors
from Sunshine et al., *Nat Commun* 14 (2023), GEO
[GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240),
starting from the processed count matrix: 27,882 cells with a guide assignment,
237 library elements, eight 10x channels, one infection experiment.

I went looking for a Perturb-seq screen in Calu-3 cells because I work on
SARS-CoV-2 in lung organoids and wanted to see whether a published screen could
say anything about the programs I see there. Reproducing the paper's calls came
first; the reagent and phenotype problems below turned up on the way, and the
organoid comparison came last.

Short on time: [docs/overview.md](docs/overview.md) is the same findings with the
figures and links to the notebook each one comes from.

## What came out of it

All 21 viral-load calls in the paper's Fig. 3A come back in the right direction:
14 reported hits at q < 0.05, and 7 reported non-hits that fail q < 0.1 here too.
The interferon panel (Fig. 5A) reproduces as well, with the four novel hits at
q < 0.01. One piece of it does not hold at the element level — GPR89B on its own
shows nothing, and only the GPR89A element and the dual-targeting element carry
the effect. `notebook/03`, `notebook/04`.

The elements carrying FURIN barely move it: +6% to −19%, against a library median
of 88% knockdown. That matters because the paper's protease-redundancy argument
treats the single-gene FURIN arm as a negative result, and a reagent that does
not knock the gene down cannot carry that weight. MEPCE and LARP1 look the same.
`notebook/02`.

Splitting the published viral-load ranking into two phenotypes — whether a cell
gets infected at all, and how much virus it carries once it is — separates
targets the single axis merges. ACE2 and the TMPRSS2 combinations block
establishment so thoroughly that too few infected cells are left to estimate
burden (12 cells for ACE2). EXOC8 moves burden with nothing on establishment
(q = 0.02 against q = 0.94). TMPRSS2 and BRD2 move both. `notebook/03`.

Some genes appear in two elements, which I had taken for two independent
reagents. For 12 of them that is true — separately cloned single-guide vectors —
and those agree at r = 0.55 on knockdown. Pooling them with the pairs that share
guide design gives 0.74, which is the flattering number.
`notebook/02`.

On the organoid question: 17 elements separate from the control band on my drug
signature, and 14 of them sit on the interferon axis — KEAP1 (ρ = 0.52), IRF9 and
IRF9+IRF3, STAT1, STAT2, TYK2, IFNAR1 and IFNAR2, GPR89A, SPNS1, BRD2 —
including all four of the paper's novel ISG regulators. PIP4K2C is not in this
library, so none of this speaks to that target directly, and the
autophagy–endolysosomal elements show no significant match (best q = 0.10, the
strongest being VAC14 at ρ = 0.22). Thin coverage is not the explanation, so
either the effect is small or it is something a transcriptional comparison cannot
see. `notebook/07`.

## Layout

```
notebook/   the analysis, in order; each one executed, with its figures
  01_qc_and_infection_calling   cell QC, lane structure, four infection classifiers
  02_knockdown                  on-target knockdown, construct concordance
  03_hit_calling                MH odds ratios, NB GLM, control-resampling null, FDR
  04_interferon                 bystander ISG scores
  05_coverage_and_resolution    subsampling curves, underpowered elements
  06_infection_state            effects by state; the cells between the two modes
  07_organoid_signatures        elements against the lung-organoid signatures
  08_sensitivity_mixed_effects  lane as a random intercept instead of a stratum
  methods.py                    the few things reused across notebooks: loading,
                                estimators, gene sets, palette
scripts/    00_download.sh, 01_build_anndata.py  (download and the one-off matrix parse)
tests/      test_estimators.py — MH against statsmodels, NB recovery, null calibration
docs/       statistical_methods.md — estimands, nulls, testing families
results/    tables/ (every number quoted above, numbered by the notebook that
            writes it) and figures/
Snakefile   the whole thing as one graph: download, build, the eight notebooks
```

## Running it

```bash
conda env create -f environment.yml && conda activate perturbseq-cov2
snakemake --cores 4                             # download, build, run notebooks
python -m pytest tests -q
```

The `Snakefile` wires the download, the build script and the eight notebooks into
one graph whose edges are the tables each notebook actually reads, so
`snakemake --cores 4 hit_calling` runs notebook 03 and only what it needs, and
editing 03 re-runs 05, 06 and 08 because they read its hit table. `snakemake -n`
says what would run; `snakemake --dag | dot -Tpng > dag.png` draws it. The
notebooks execute in place, so the copies in the repository keep their outputs.

To work through them by hand instead:

```bash
bash scripts/00_download.sh data/raw            # ~1.2 GB from GEO
python scripts/01_build_anndata.py --raw data/raw --out data
jupyter lab notebook/                           # run 01 through 08 in order
```

`01_build_anndata.py` streams the 337-million-entry matrix and caches the result,
so it takes about 11 minutes once and seconds thereafter. The notebooks run in a
few minutes each, apart from `07`, which needs the organoid differential-expression
tables. Those are not public and are not in this repository; point `ALO_LFC` at
your own copy.

## Data

GSE208240, deposited by the authors of the original study. Nothing in `data/` is
committed. I start from their per-cell guide assignments, so this is a reanalysis
from the processed matrix rather than from FASTQs.

## Citation

> Sunshine S, Puschnik AS, Replogle JM, Laurie MT, Liu J, Zha BS, Nuñez JK,
> Byrum JR, McMorrow AH, Frieman MB, Winkler J, Qiu X, Rosenberg OS, Leonetti MD,
> Ye CJ, Weissman JS, DeRisi JL, Hein MY. Systematic functional interrogation of
> SARS-CoV-2 host factors using Perturb-seq. *Nature Communications* 14 (2023).
> doi:10.1038/s41467-023-41788-4

## Licence

MIT, for the code only (see [LICENSE](LICENSE)). The GSE208240 data and the
original figures remain under their own terms.
