# Statistical methods

What is estimated, against what null, and what each result licenses you to say.

---

## 1. The experimental unit

The unit is the **library element**, not the cell and not the gene.

A cell is not an independent replicate: cells sharing an element share a lentiviral
integration history, a selection history and a droplet batch. Treating 27,882 cells
as 27,882 replicates is the commonest way a Perturb-seq analysis manufactures
significance.

The element is not the gene either. The library has 237 elements for 181 target
genes:

| kind | n | note |
|---|---|---|
| single-gene elements | 193 | two sgRNAs in one dual-guide vector |
| two-gene combination elements | 22 | e.g. `CRISPRi_172_ACE2_TMPRSS2` |
| non-targeting elements | 22 | the empirical null |

### 1a. What counts as an independent construct

Because the standard element carries **both** sgRNAs for its gene in a single
dual-guide vector, the usual split-guide concordance check is unavailable for most
of the library: the guide pair *is* the reagent.

Two elements hitting the same gene are therefore not automatically two independent
constructs, and `04_guide_qc.py` keeps the categories apart rather than pooling
them (`results/tables/04_construct_concordance.tsv`):

| category | n genes | what it is | knockdown agreement |
|---|---|---|---|
| **independent vectors** | **12** | elements 216–239: 24 separately cloned single-guide vectors, two per gene, different guide sequences | **Pearson r = 0.55, Spearman ρ = 0.36, median \|difference\| 0.15** |
| single vs combination | 24 | a single-gene element and a two-gene element containing it, which share guide design for the shared gene | r = 0.55, ρ = 0.20 |
| three or more elements | 5 | genes appearing in several elements | r = 0.98 |

**The reproducibility figure this repository quotes is r = 0.55 on the 12-gene
independent-vector panel.** Pooling all 41 multi-element genes into one correlation
gives r = 0.74, but that number is inflated by the third category and mixes
reagent-to-reagent comparison with measuring one reagent twice. It should not be
reported as construct reproducibility, and earlier drafts of this repository did so
in error.

Twelve genes is a small panel, and the honest reading of r = 0.55 with a median
absolute difference of 15 percentage points is that **knockdown efficiency varies
materially between independently cloned guides against the same gene**. That is a
reason to treat element-level, not gene-level, results as primary.

## 2. Blocking: the 10x channel

`gem_group` runs 1–8 and is a real experimental factor, not a label. Channels 5–8
show about 1.4× the infection rate and about 1.8× the median viral fraction of
channels 1–4, reflecting the paper's pulsed and non-pulsed infection conditions,
which it combined after finding no significant difference between them.

Every model carries lane as a stratum or covariate, and every null is matched on
lane. An element over-represented in channels 5–8 would otherwise look proviral.
`tests/test_estimators.py::test_stratification_actually_matters` is a worked case:
with no true effect and the element confined to the high-infection lanes, the
pooled estimator reports a spurious 1.6× odds ratio and the stratified estimator
correctly reports null.

## 3. The phenotypes

### 3a. Infection establishment

**Estimand.** The odds of a cell being infected rather than bystander, for cells
carrying element *e*, relative to cells carrying a non-targeting element, within
10x channel.

**Estimator.** The Mantel–Haenszel odds ratio across the eight channels, with the
Robins–Breslow–Greenland standard error. Verified against
`statsmodels.stats.contingency_tables.StratifiedTable` to 1e-10 on the point
estimate and 1e-8 on the standard error, and checked for recovery of a known odds
ratio under simulation (`tests/test_estimators.py`).

A logistic GLM with lane as a covariate is fitted in parallel; the two agree at
r = 0.9999, so the stratified estimator is used for speed and because it needs no
per-element lane parameters.

Borderline cells (mixture posterior between 0.1 and 0.9, 9.7% of cells) are
excluded. They are exactly the cells whose call would be depth-sensitive.

### 3b. Viral burden given infection

**Estimand — stated precisely, because the denominator is not neutral.** Viral UMIs
per **host transcript captured**, among infected cells. This is *not* absolute viral
RNA abundance per cell. Host UMI is itself lowered by infection (host shutoff cuts
median host UMIs roughly elevenfold) and could in principle be moved by a
perturbation, so the exposure is a measurement quantity rather than a fixed offer.
A perturbation that changed host capture without changing viral replication would
move this estimand.

**Estimator.** A negative-binomial GLM on raw viral counts,

```
viral_umi ~ element + C(gem_group),  offset = log(host_umi)
```

with the dispersion estimated once on the non-targeting cells (α = 1.66) and held
fixed, so the 149 per-element fits are comparable. Recovery of a known rate ratio is
tested under simulation.

**Why not total UMI.** Total UMI contains the viral counts being modelled and is
itself depressed by infection, so conditioning on it means conditioning on the
outcome. `tests/test_estimators.py::test_offset_choice_changes_the_answer`
demonstrates the bias on simulated data where the truth is known.

**Sensitivity to the denominator** (`07_offset_sensitivity_summary.tsv`), as rank
correlation of the per-element effect against the primary specification, and
overlap of the 20 most negative elements:

| specification | n elements | Spearman ρ | top-20 overlap |
|---|---|---|---|
| `log(host UMI)` — primary | 170 | — | — |
| `log(host genes detected)` | 170 | 0.97 | 16/20 |
| `log(total UMI)` — the rejected spec | 170 | 0.87 | 13/20 |
| no offset, infected cells within a 2-fold band of median host UMI | 88 | 0.56 | 3/20 |

The first alternative is the informative one: a capture proxy that saturates
differently from UMI count gives essentially the same answer, so the result is not
an artefact of how host capture is measured. The depth-band row is weak evidence in
either direction and is reported rather than hidden: restricting to a two-fold band
keeps 3,126 of 8,136 infected cells and only 289 non-targeting ones, so roughly half
the elements become unfittable and the rest are estimated on a fraction of the data.
It is underpowered, not contradictory.

Keeping 3a and 3b apart matters: a knockdown that blocks entry moves 3a; one that
limits replication moves 3b. ACE2 and TMPRSS2 move 3a sharply and leave too few
infected cells to test 3b at all, which is itself the result.

### 3c. Interferon tone

**Estimand.** The shift in a cell's ISG score relative to non-targeting cells,
**in bystander cells only**. Infected cells suppress interferon signalling, so a
score on a mixed population is dominated by each element's infected fraction.

**Estimator.** `scanpy.tl.score_genes` over 55 type-I ISGs on host-only CPM (the
denominator excludes viral and lentiviral counts), then a lane-weighted mean
difference against the non-targeting cells of the same lane.

**The self-containment correction.** Twelve of the 55 ISGs are themselves CRISPRi
targets in this library (BST2, DDX58, EIF2AK2, IFIH1, IFITM1/2/3, IRF9, ISG15, LY6E,
STAT1, STAT2). Knocking one down mechanically lowers any score containing it. Every
result is reported with and without them; the repository quotes the target-free
version.

## 4. The null — what it is, exactly

A Wald p-value from any model above assumes cells are independent. They are not
(§1). Significance therefore comes from a **control-based null** built from the 22
non-targeting elements.

**This is not a permutation of perturbation labels.** No label is shuffled. The
procedure constructs synthetic null elements out of control cells. Earlier drafts
of this repository called it a "permutation null"; that was imprecise, and the
column names are now `p_null` / `q_null`.

Two constructions are used, because the two statistics need different ones:

**Phenotype A — lane-matched binomial null** (`_hitcalling_core.mh_null_draws`).
For element *e* with n<sub>ℓ</sub> cells in lane ℓ, each draw generates
a<sub>ℓ</sub> ~ Binomial(n<sub>ℓ</sub>, p̂<sub>ℓ</sub>), where p̂<sub>ℓ</sub> is the
pooled non-targeting infection rate in that lane, and recomputes the MH odds ratio
against the observed non-targeting counts.

**Phenotypes B and 3c — non-parametric bootstrap of control cells**
(`_hitcalling_core.resample_delta_null`). Each draw samples, **with replacement**,
the same number of cells in each lane as the element has, from that lane's pooled
non-targeting cells, and recomputes the lane-weighted mean difference.

Answering the specific questions an interviewer will ask:

| question | answer |
|---|---|
| Sampled with replacement? | Phenotypes B and 3c: yes, `rng.integers` indexes into the lane pool independently each time. Phenotype A draws binomial counts rather than cells, so replacement does not apply. |
| Sampled independently within each lane? | Yes. Each lane is drawn separately and the draws are summed into one statistic. |
| Can the same control cell appear in several pseudo-elements, or twice in one? | Yes to both. The control pool is fixed and reused; this is a bootstrap, not a partition. |
| Plus-one correction? | Yes: p = (1 + #{\|null\| ≥ \|observed\|}) / (B + 1), so p is never 0 and is floored at 1/2001 with B = 2,000. |
| One null per phenotype? | Yes. Each phenotype builds its own null from its own eligible control cells, and its own statistic. |
| FDR across what? | Benjamini–Hochberg within each phenotype, across that phenotype's **targeting** elements only (§5). |
| Is the null re-matched per element? | Yes. Every element gets its own 2,000 draws at its own per-lane cell counts. |

**Calibration.** Running the real non-targeting elements through the pipeline as if
they were perturbations gives 0 of 20 below p = 0.05 and a median p of 0.44. The
Wald p-values are also well behaved (0 below 0.05, median 0.42), so at these cell
counts the independence assumption is not visibly costing anything; the
control-based null is still used because it does not rest on that assumption.
`tests/test_estimators.py` additionally checks both null constructions against
simulated data where the truth is null, and requires a false-positive rate near 5%
and a mean p near 0.5.

**Robustness.** Repeating phenotype A under a different infection classifier gives
effect sizes correlated at r = 0.99, with 23 elements significant under both rules
and 31 under either (Jaccard 0.74). Elements near the q = 0.1 boundary move, which
is what a boundary means.

## 5. Testing families, and why the denominators differ

FDR depends on the family, so each analysis's eligible population is defined here.
Controls are **excluded from the discovery family**: they are not candidate hits,
and including them would both inflate the family and treat a control as a
hypothesis. They are judged on raw p in the calibration table instead.

| analysis | eligibility rule | targeting elements (FDR family) | NTC elements (calibration only) | excluded, and why |
|---|---|---|---|---|
| library as deposited | — | 215 | 22 | — |
| infection establishment | ≥ 55 cells | **164** | 20 | 51 elements with < 55 cells (47 targeting + 2 NTC below the cut, plus 2 NTC absent); mostly essential genes |
| viral burden | ≥ 20 infected cells | **149** | 21 | 67 elements without 20 infected cells, including the entry factors whose knockdown leaves almost no infected cells |
| bystander ISG | ≥ 30 bystander cells | **170** | 20 | 47 elements without 30 bystander cells |
| published-list comparison | the 21 targets named in the paper | 21 rows | — | not an FDR family; it reads q values off the infection table |

The 184 figure that appears next to phenotype A is 164 + 20: the number of elements
the script *evaluates*, controls included. The FDR family is the 164. The
published-list comparison is a lookup of 21 named targets, not a separate test.

## 6. Coverage and resolution

`08_representation.py` reports, from the null alone, the smallest odds ratio
separable from zero at a given cell count: 2.42× at 25 cells, 1.73× at 60, 1.48× at
120, 1.31× at 240.

Subsampling real hits shows that the **point estimate stays close to the full-data
estimate** as cells are removed, while **detection probability falls sharply for
small effects**. Both statements are about behaviour relative to the full-data
estimate on these cells; neither is a claim of unbiasedness with respect to the true
biological effect, which this design cannot assess. Quantitatively, subsampling to
25 cells per element shifts the median point estimate by −0.01 log OR while the
resampling spread around it is wide (median SD 0.49 log OR, against 0.08 at 160
cells), so "stable" means stable in expectation, not precise.

47 of 237 elements fall below the 55-cell cut, mostly essential genes, so their
under-representation is systematic: absence from a hit list carries no information
about those targets.

## 7. What this design cannot support

- **No claim of biological reproducibility.** The eight channels are 10x lanes from
  one infection experiment. The null controls cell- and lane-level sampling noise.
  It says nothing about whether the screen would replicate in a second infection.
- **No mixed-effects model with a biological random effect**, for the same reason:
  no grouping factor above the lane has an estimable variance. Lane enters as a
  fixed effect or a stratum, which is what the design supports.
- **No formal power analysis.** §6 describes this screen's resolution given its
  coverage, not power for a future experiment.
- **No unbiasedness claim** from the subsampling in §6; see the wording there.
- **No causal claim from the extension.** `09_infected_vs_bystander.py` compares
  effect profiles between infection states descriptively. An interaction test needs
  a perturbation-by-state term and a variance component this design does not
  provide.
- **Gene-level conclusions are weaker than element-level ones**, given §1a.

## 8. Verification

`tests/test_estimators.py` (12 tests, run with `python -m pytest tests -q`) checks:
MH against statsmodels and against simulation; that stratification fixes a
confounded design a pooled estimator fails; NB recovery of a known rate ratio; that
the total-UMI offset is measurably the worse estimator; BH against
`statsmodels.multipletests`, including NaN handling; the plus-one correction; and
calibration of both null constructions on data generated under the null.
