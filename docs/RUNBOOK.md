# Runbook

Every command, in order, with what it should print. Written so the analysis can be
rerun from nothing and the numbers checked against the README.

## 0. Environment

```bash
conda env create -f environment.yml
conda activate perturbseq-cov2
```

Needs about 10 GB of free disk (1.2 GB download, ~5 GB unpacked, ~1.4 GB of derived
objects) and about 8 GB of RAM for step 01.

## 1. Download

```bash
bash scripts/00_download.sh data/raw
```

Pulls two tarballs from the GEO FTP site, writes `md5sums.txt`, unpacks them and
creates the `data/raw/matrix_dir` and `data/raw/viral_references` symlinks that the
rest of the pipeline expects.

Expected:

```
barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz  cell_identities.csv
```

## 2. Build the object

```bash
python scripts/01_build_anndata.py --raw data/raw --out data
```

Streams the 337-million-entry matrix in chunks, keeps the 44,273 barcodes that are
either guide-assigned or guide-free (the ambient reference), and writes
`data/perturbseq.h5ad`. The subset matrix is cached as `data/X_kept.npz`, so a
rerun skips the parse.

The script asserts rather than reports: it stops if the matrix is not
33,549 features, if the exact-match cell count is not 27,882, or if the number of
library elements is not 237. Expected output:

```
features by kind: {'host': 33538, 'viral': 9, 'lenti': 2}
assigned cells      : 27882
wt candidate cells  : 16391
library elements with >=1 assigned cell: 237
non-targeting elements: 22
median_cells_per_element   113
elements_under_55_cells     47
median_host_umi_assigned    18066.5
```

If you wrote the `.h5ad` somewhere other than `data/`, point the rest of the
pipeline at it:

```bash
export PERTURBSEQ_H5AD=/path/to/perturbseq.h5ad
```

## 3. The analysis, in order

Order matters: 03 is upstream of 04 because knockdown is measured in bystander
cells; 07 reads the cell states from 03; 08 reads the hit list from 07; 10 reads
only tables.

```bash
python scripts/02_qc.py
python scripts/03_infection.py
python scripts/04_guide_qc.py
python scripts/05_viral_load.py
python scripts/06_isg.py
python scripts/07_hit_calling.py
python scripts/08_representation.py
python scripts/09_infected_vs_bystander.py
python scripts/10_figures.py
python scripts/11_mixed_effects.py
python -m pytest tests -q
```

### What each step should report

**02_qc** — median host UMI 18,066 in assigned cells; lane median host-UMI spread
1.23×; 237 elements (22 NTC, 193 single-gene, 22 multi-gene); 47 elements under 55
cells. The line to look at is the correlation between total and host UMI: 0.991
overall but 0.757 in cells where the virus takes more than half the library, which
is why host UMI is the offset everywhere downstream.

**03_infection** — 705 ambient-reference cells; mixture modes at viral fraction
0.0023 and 0.100 with a crossover at 0.0088; 29.2% infected, 61.1% bystander, 9.7%
borderline. The three threshold-free classifiers agree at 99.3% or better; the paper's
literal mean + 2 SD rule calls only 6.2% infected on the reconstructable
reference and is reported as a documented limitation of that reference. Downsampling the whole
library to a tenth changes 0.45% of calls.

**04_guide_qc** — 177 of 181 target genes present in the reference (ADCK4,
C18orf8, TCEB1, TCEB2 are absent under those symbols). Median knockdown 88%, 76% of
element–gene pairs above 75%. Construct concordance is reported by pair type: the
12 genes carried by two separately cloned single-guide vectors agree at r = 0.55
(Spearman 0.36, median |difference| 0.148), the 24 single-vs-combination pairs at
r = 0.55, and the 5 genes in three or more elements at r = 0.97. Only the first is
quoted. Four element–gene pairs have negative knockdown; FURIN is the one that
matters.

**05_viral_load** — NTC baseline infected fraction 0.347. The paper's 14 hits
occupy ranks 5–33 of 237 in the expected direction, SEC62 rank 231, and its 7
reported non-hits land between ranks 79 and 221.

**06_isg** — 55 of 55 ISGs present; 12 of them are also library targets and are
dropped from the clean score. Bystander mean ISG 2.41 against 1.04 in infected
cells. 21 of 170 targeting elements significant at q < 0.1.

**07_hit_calling** — 184 elements testable for infection (164 targeting + 20
non-targeting); BH runs over the 164 targeting elements and 30 are significant at
q < 0.1. Burden: 149 targeting elements, 5 significant. MH and logistic-GLM effect
sizes agree at r = 0.9999. NB dispersion α = 1.66. **21 of 21 published calls
reproduced.** Zero of the 20 testable non-targeting elements reach p < 0.1 for
infection establishment (median p = 0.44). Offset sensitivity: Spearman 0.967 vs
detected-gene count, 0.870 vs total UMI, 0.559 for depth-banding with no offset.

**08_representation** — smallest detectable odds ratio 2.42× at 25 cells, 1.73× at
60, 1.39× at 160, 1.31× at 240. 73 of 237 elements fail at least one phenotype's coverage
requirement.

**09_infected_vs_bystander** — 133 elements have ≥ 25 cells in both states.
Program-level effect correlations between states: ISG 0.74, MHC-I 0.54, NF-κB 0.34,
unfolded protein 0.20.

**10_figures** — seven PNGs in `results/figures/`.

**11_mixed_effects** — lane random-intercept GLMM for the 184 elements over the
55-cell cut, about 20 s. Effect sizes correlate with the Mantel-Haenszel
estimates at r = 0.9999 (median |difference| 0.006 log OR); 28 hits against 30,
27 shared, Jaccard 0.87; fitted lane SD 0.33 on the log-odds scale.

**11_mixed_effects** — lane random-intercept GLMM, ~20 s; r = 0.9999 with the
Mantel-Haenszel estimates.

**12_borderline_cells** — the 2,696 cells between the two viral-load modes:
entry-factor knockdowns deplete them, ORF1ab:N ratio 0.100 vs 0.068 in infected
cells, no shutoff at matched depth. The published interferon result is unchanged
by how they are classified.

**13_infection_programs** — autophagy/lysosome programme scores, raw and after
downsampling every cell to 1,500 host UMI, against ten random gene sets.

**14_organoid_signature_match** — needs the organoid DE tables, which are not
public; 17 elements match the drug signature, 0 match the cross-run infection
signature.

**pytest** — `python -m pytest tests -q`: 12 estimator tests, ~10 s.

## 4. Runtime

Measured on 4 cores:

| step | time |
|---|---|
| 00 download | 2–10 min, network bound |
| 01 build (first run, parses the 337M-entry matrix) | ~11 min |
| 01 build (rerun, uses the cache) | ~1 min |
| 02, 03, 04 | 8–9 s each |
| 05 | 1 s |
| 06 isg (2,000 permutations × two scores) | ~60 s |
| 07 hit calling (2,000 permutations × two phenotypes + 334 GLM fits) | ~25 s |
| 08, 09 | 8–10 s each |
| 10 figures | ~5 s |

So about 14 minutes end to end on a cold start, and under 3 minutes for everything
downstream of the matrix parse.

Peak memory is in step 01, at roughly 6 GB while the subset matrix is assembled.

## 5. Determinism

Every random draw is seeded: `numpy.random.default_rng(0)` for permutations and
subsampling, `random_state=0` for the Gaussian mixture and for `score_genes`.
Rerunning on the same package versions reproduces the tables exactly.

Permutation p-values are floored at 1/2001, so many hits report the same p. Compare
hit *sets* and effect sizes across environments, not p-values at the floor.

## 6. If something fails

- **`could not find the matrix directory`** — `00_download.sh` did not finish
  unpacking. Check that `data/raw/matrix_dir` resolves.
- **An assertion in 01** — the GEO files changed or were partially downloaded. Check
  `data/raw/md5sums.txt` and re-download.
- **`cached matrix shape mismatch`** — `data/X_kept.npz` is stale, e.g. the cell
  filter changed. Delete it and rerun 01.
- **A script cannot find a table** — the steps are not independent. Run them in the
  order in §3.
