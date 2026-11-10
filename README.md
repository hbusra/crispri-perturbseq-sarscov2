# crispri-perturbseq-sarscov2

Reanalysis of the Calu-3 CRISPRi Perturb-seq screen of SARS-CoV-2 host factors
from Sunshine et al., *Nat Commun* 14 (2023), GEO
[GSE208240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208240),
starting from the processed count matrix: 27,882 cells with a guide assignment,
237 library elements, eight 10x channels, one infection experiment.

Three things, in order: reproduce the published calls, test the ones that do not
hold up, and use the perturbation data against signatures from my own
SARS-CoV-2 lung organoid work.

## What came out of it

**The published calls reproduce.** All 21 viral-load calls in the paper's Fig. 3A
come back in the right direction — 14 reported hits at q < 0.05, 7 reported
non-hits failing q < 0.1 — and the interferon panel (Fig. 5A) reproduces with all
four novel hits at q < 0.01. `notebook/03`, `notebook/04`.

**The FURIN element does not knock FURIN down** (+6% to −19% across the three
elements carrying it, library median 88%), so the single-gene FURIN arm of the
paper's protease-redundancy argument is an uninterpretable reagent rather than a
negative result. MEPCE and LARP1 likewise. `notebook/02`.

**Establishment and burden are different phenotypes.** The published single
viral-load axis merges them. Split, ACE2 and the TMPRSS2 combinations block
establishment but leave burden unchanged, while MARK2 and EXOC8 lower burden
without clearing the establishment cut. `notebook/03`.

**Two elements per gene is not two independent reagents.** The 12 genes carried
by separately cloned single-guide vectors agree at r = 0.55 on knockdown; pooling
them with pairs that share guide design gives 0.74. `notebook/02`.

**Against the organoid signatures**, 17 elements separate from the control band
on the drug signature and 15 of them are the interferon axis — KEAP1 (ρ = 0.52),
IRF9+IRF3, STAT2, SPNS1, TYK2, IFNAR1/2, GPR89A/B, STAT1, BRD2 — including all
four of the paper's novel ISG regulators. PIP4K2C is not in this library, so
nothing here validates that target; the autophagy–endolysosomal elements show no
match (best q = 0.10), which is what you expect if the mechanism is
post-transcriptional. `notebook/07`.

## Layout

```
notebook/   the analysis, in order; each one executed, with output
  01_qc_and_infection_calling   cell QC, lane structure, four infection classifiers
  02_knockdown                  on-target knockdown, construct concordance
  03_hit_calling                MH odds ratios, NB GLM, control-resampling null, FDR
  04_interferon                 bystander ISG scores
  05_coverage_and_resolution    subsampling curves, underpowered elements
  06_infection_state            effects by state; the cells between the two modes
  07_organoid_signatures        elements against the lung-organoid signatures
  08_sensitivity_mixed_effects  lane as a random intercept instead of a stratum
  09_figures                    the panels, drawn from results/tables only
src/perturbseq/   estimators, loading, gene sets, palette, figure functions
scripts/    00_download.sh, 01_build_anndata.py  (download and the one-off matrix parse)
tests/      test_estimators.py — MH against statsmodels, NB recovery, null calibration
docs/       statistical_methods.md — estimands, nulls, testing families
results/    tables/ (every number quoted above) and figures/
```

## Running it

```bash
conda env create -f environment.yml && conda activate perturbseq-cov2
bash scripts/00_download.sh data/raw            # ~1.2 GB from GEO
python scripts/01_build_anndata.py --raw data/raw --out data
jupyter lab notebook/                           # run 01 through 09 in order
python -m pytest tests -q
```

`01_build_anndata.py` streams the 337-million-entry matrix and caches the result,
so it takes about 11 minutes once and seconds thereafter; the notebooks take
about 3 minutes end to end, apart from `07`, which needs the organoid
differential-expression tables. Those are not public and are not in this
repository; point `ALO_LFC` at your own copy.

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
