# Statistical methods

What is estimated, against what null, and what each result licenses you to say.

---

## 1. The experimental unit

The unit is the **library element**, not the cell and not the gene.

A cell is not an independent replicate: cells sharing an element share a lentiviral
integration history, a selection history and a droplet batch. Treating 27,882 cells
as 27,882 replicates is the single most common way a Perturb-seq analysis
manufactures significance.

The element is not the gene either. The library has 237 elements for 181 target
genes:

| kind | n | note |
|---|---|---|
| single-gene elements | 193 | two sgRNAs in one dual-guide vector |
| two-gene elements | 22 | e.g. `CRISPRi_172_ACE2_TMPRSS2` |
| non-targeting elements | 22 | the empirical null |

Because the standard element carries **both** sgRNAs for its gene in one vector,
the usual Perturb-seq split-guide concordance check is not available for most of
the library — the pair *is* the reagent. The exception is elements 216–239, which
are 12 genes carried by two independent single-guide vectors. Together with genes
appearing in both a single and a dual element, this gives 41 genes with two
independent elements, and they are used as the internal reproducibility check
(`04_guide_qc.py`, knockdown r = 0.74).

## 2. Blocking: the 10x channel

`gem_group` runs 1–8 and is a real experimental factor, not a label. Channels 5–8
show about 1.4× the infection rate and about 1.8× the median viral fraction of
channels 1–4, reflecting the paper's pulsed and non-pulsed infection conditions,
which it combined after finding no significant difference between them.

Every model below carries lane as a stratum or covariate, and every permutation is
matched on lane. An element that happened to be over-represented in channels 5–8
would otherwise look like a proviral hit.

## 3. The phenotypes

### 3a. Infection establishment

**Estimand.** The odds of a cell being infected rather than bystander, for cells
carrying element *e*, relative to cells carrying a non-targeting element, within
10x channel.

**Estimator.** The Mantel–Haenszel odds ratio across the eight channels, with the
Robins–Breslow–Greenland standard error. A logistic GLM with lane as a covariate is
fitted in parallel as a check; the two agree at r = 0.9999, so the stratified
estimator is used for speed and because it needs no per-element lane parameters.

Borderline cells (mixture posterior between 0.1 and 0.9, 9.7% of cells) are
excluded from this contrast. They are exactly the cells whose call would be
depth-sensitive, and keeping them would pull every element toward the mean.

### 3b. Viral burden given infection

**Estimand.** Viral UMIs per host UMI, among infected cells only.

**Estimator.** A negative-binomial GLM on raw viral counts,

```
viral_umi ~ element + C(gem_group),  offset = log(host_umi)
```

with the dispersion estimated once on the non-targeting cells (α = 1.66) and held
fixed, so that 149 per-element fits are comparable to each other. A lane-weighted
difference in mean log₁₀ viral rate is computed alongside and is what the
permutation null is built on; the two agree at r = 0.83.

**Why the offset is host UMI and not total UMI.** Total UMI contains the viral
counts being modelled, and infection itself depresses host transcription (host
shutoff, roughly an 11-fold drop in median host UMIs). Conditioning on total UMI
would therefore condition on the outcome: it would partly absorb the very effect
being estimated and bias every element toward null. Host UMI is a denominator the
perturbation can move but the viral counts do not enter.

This is also why the two phenotypes are kept apart. A knockdown that blocks entry
moves 3a; one that limits replication moves 3b. ACE2 and TMPRSS2 move 3a sharply
and have too few infected cells left to test 3b at all — which is itself the
result.

### 3c. Interferon tone

**Estimand.** The shift in a cell's ISG score relative to non-targeting cells,
**in bystander cells only**.

Infected cells suppress interferon signalling, so a score computed on a mixed
population is dominated by each element's infected fraction rather than by its
effect on interferon. The restriction is the paper's and is kept.

**Estimator.** `scanpy.tl.score_genes` over 55 type-I ISGs, computed on host-only
CPM (the denominator excludes viral and lentiviral counts), then a lane-weighted
mean difference against the non-targeting cells of the same lane.

**The self-containment correction.** Twelve of the 55 ISGs are themselves CRISPRi
targets in this library (BST2, DDX58, EIF2AK2, IFIH1, IFITM1/2/3, IRF9, ISG15,
LY6E, STAT1, STAT2). Knocking one down mechanically lowers any score containing it,
which looks like a signalling phenotype and is not. Every result is reported with
and without those genes; the repository quotes the target-free version.

## 4. The null

A Wald p-value from any of the models above assumes cells are independent. They are
not (§1). So significance comes from an empirical null built out of the 22
non-targeting elements:

> For element *e* and each of B = 2,000 draws, construct a pseudo-element with the
> same number of cells **in each lane** as *e*, drawn from the pooled non-targeting
> cells of that lane, and recompute the statistic. The two-sided p is the fraction
> of draws at least as extreme, floored at 1/(B+1), then Benjamini–Hochberg
> adjusted across the family of elements tested for that phenotype.

The null is matched on size and on lane composition by construction, so an element
with 60 cells is judged against a 60-cell null, and an element concentrated in
channels 5–8 is judged against a null drawn from channels 5–8.

**Calibration.** Running the 22 real non-targeting elements through the pipeline as
if they were perturbations gives 0 elements below p = 0.05 and a median p of 0.44.
The Wald p-values are also well behaved here (0 below 0.05, median 0.42), so at
these cell counts the independence assumption is not visibly costing anything. The
permutation null is still what calls hits, because it does not rest on that
assumption.

**Robustness.** Repeating the infection phenotype under a different infection
classifier (the depth-aware per-gene ambient rule rather than the mixture) gives
effect sizes correlated at r = 0.99 and a hit-set Jaccard of 0.60. The Jaccard is
the honest number: 18 elements are significant under both rules, 30 under either.
Elements near the q = 0.1 boundary move, which is what a boundary means.

## 5. Multiple testing

Benjamini–Hochberg within each phenotype family, over the elements that pass that
phenotype's coverage requirement:

| phenotype | requirement | elements tested |
|---|---|---|
| infection establishment | ≥ 55 cells | 184 (164 targeting + 20 NTC) |
| viral burden | ≥ 20 infected cells | 170 |
| bystander ISG | ≥ 30 bystander cells | 190 |

Families are kept separate rather than pooled, because they answer different
questions and a reader interested in interferon should not pay a multiplicity price
for the viral-load tests.

## 6. Coverage and resolution

`08_representation.py` reports, from the null alone, the smallest odds ratio
separable from zero at a given cell count: 2.42× at 25 cells, 1.75× at 60, 1.48× at
120, 1.31× at 240. Subsampling real hits confirms that effect *estimates* stay
unbiased down to 25 cells while *detection* degrades sharply for small effects.

47 of 237 elements fall below the 55-cell cut. They are mostly essential genes, so
their under-representation is systematic: absence from a hit list carries no
information about those targets.

## 7. What this design cannot support

- **No claim of biological reproducibility.** The eight channels are 10x lanes from
  one infection experiment. The null controls cell- and lane-level sampling noise.
  It says nothing about whether the screen would replicate in a second infection.
- **No mixed-effects model with a biological random effect**, for the same reason:
  there is no grouping factor above the lane whose variance is estimable. Lane
  enters as a fixed effect or a stratum, which is what the design supports.
- **No formal power analysis.** §6 is a statement about this screen's resolution
  given its coverage, not an estimate of power for a future experiment.
- **No causal claim from the extension.** §3 of `09_infected_vs_bystander.py`
  compares effect profiles between infection states descriptively. A real
  interaction test needs a perturbation-by-state term and a variance component this
  design does not provide.
