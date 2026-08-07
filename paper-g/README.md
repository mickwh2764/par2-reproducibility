# Paper G — Fibonacci structure and phase-gated AR(2) dynamics

Reproduction bundle for *A Time-Domain Analogue to Fibonacci Structure via
Phase-Gated AR(2) Dynamics* (The Fibonacci Quarterly, accepted and in press),
a reply to Boman's spatial-recursion argument.

This is an **independent reimplementation in Python**, not the code that
produced the manuscript's tables — those were computed by the platform's
TypeScript. Where the two disagree, this README says so rather than tuning the
Python to match. See [Agreement with the published numbers](#agreement-with-the-published-numbers).

## Contents

```
code/download_data.py   fetch the GEO series matrix and platform annotation
code/geo.py             parse GEO matrices, map probes to gene symbols
code/hierarchy.py       AR(2) fit of every gene in every tissue, category tests
code/fibonacci_null.py  phi-proximity of AR(2) coefficient ratios vs a null
code/permutation.py     seeded permutation tests (clock median, 1/phi zone)
code/theory_checks.py   numerical check of the paper's analytic propositions
code/gene_categories.json  the gene -> functional category assignment used throughout
data/accessions.tsv     public GEO accessions with URLs and sampling metadata
data/checksums.txt      SHA-256 of the exact upstream files used
results/                committed expected outputs
verify.py               rerun everything and diff against results/
tests/                  unit tests, no network access needed
```

Raw data is **not** committed. `download_data.py --verify` fetches it from NCBI
and checks it against `data/checksums.txt`, so the exact upstream release is
pinned without vendoring 27 MB of someone else's data.

## Data

| Accession | Organism | Design | Used for |
| --- | --- | --- | --- |
| [GSE54650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54650) | *Mus musculus* | 12 tissues, 24 timepoints at 2 h, GPL6246 | everything below |

Probes are collapsed to gene symbols by keeping, per symbol, the probe with the
largest interquartile range across samples; probes mapping to multiple symbols
are dropped. That leaves 20,955 genes × 288 samples (12 tissues × 24 timepoints).

## Claim → script → result

| Claim in the paper | Script | Output | What we get |
| --- | --- | --- | --- |
| Boman's spatial matrix *M* is the AR(2) companion matrix at (1,1) | `theory_checks.py` | `results/theory_checks.json` | holds exactly |
| The Fibonacci point (1,1) is non-stationary with spectral radius φ | `theory_checks.py` | " | holds exactly |
| The equal-coefficient ray leaves stationarity at *c* = 1/2, so *c* = 1 is twice as far | `theory_checks.py` | " | holds exactly |
| φ identities and the φ² memory sum (Props S1, S2, Cor S2) | `theory_checks.py` | " | hold exactly |
| The empirical stable band (0.52, 0.72) admits Boman's p=0/1 and p=2 but not p=3 | `theory_checks.py` | " | holds |
| Clock genes are more persistent than their targets | `hierarchy.py` | `results/hierarchy_summary.json` | clock > target in **12/12** tissues, median gap +0.238 |
| Clock > target > background | `hierarchy.py` | " | full ordering in **8/12** tissues (see caveats) |
| Clock persistence is not a property of arbitrary gene sets | `permutation.py` | `results/permutation_tests.json` | significant at *p* < 0.05 in **11/12** tissues |
| φ-proximity is not enriched once the null is stationarity-constrained | `fibonacci_null.py` | `results/fibonacci_null.json` | observed 0–5.3% per tissue vs 4.8% expected — no enrichment |
| The 1/φ zone is not enriched in classified genes | `permutation.py` | " | 43.7% vs 40.8% expected, *p* = 0.22 |

All 7 analytic checks pass. The two φ results are **negative**, and that is the
point of the bundle: the paper states that its φ-zone enrichment should be
treated as an upper bound pending a stringent null, and that the genome-wide
analysis was not significant. This code makes that testable rather than
asserted.

## Eigenvalue hierarchy (GSE54650)

Median |λ| per category, per tissue:

| Tissue | Clock | Target | Background | Gap | Clock > target > background |
| --- | ---: | ---: | ---: | ---: | :---: |
| adrenal | 0.6570 | 0.3573 | 0.4165 | +0.2997 | no |
| aorta | 0.6773 | 0.4363 | 0.4193 | +0.2410 | yes |
| brainstem | 0.6157 | 0.4054 | 0.4395 | +0.2103 | no |
| brown fat | 0.7102 | 0.3895 | 0.4232 | +0.3207 | no |
| cerebellum | 0.5776 | 0.4864 | 0.4404 | +0.0912 | yes |
| heart | 0.7292 | 0.4937 | 0.4376 | +0.2355 | yes |
| hypothalamus | 0.4468 | 0.4379 | 0.4006 | +0.0089 | yes |
| kidney | 0.8195 | 0.5300 | 0.4610 | +0.2895 | yes |
| liver | 0.6467 | 0.5285 | 0.4963 | +0.1182 | yes |
| lung | 0.8591 | 0.5290 | 0.4904 | +0.3301 | yes |
| muscle | 0.6313 | 0.4858 | 0.4529 | +0.1455 | yes |
| white fat | 0.6828 | 0.3974 | 0.4155 | +0.2854 | no |

Per-gene fits are in `results/gene_lambdas_<tissue>.csv.gz`; per-category
statistics, including Mann–Whitney *p* against the genome-wide background with
Benjamini–Hochberg correction within each tissue, are in
`results/category_tests.csv`.

## Caveats

- **The full three-level hierarchy does not hold everywhere.** In adrenal,
  brainstem, brown fat and white fat the *target* median falls **below** the
  genome-wide background, so the ordering is clock > background > target. The
  clock > target part holds in all twelve tissues. The paper's headline liver
  numbers are unaffected.
- **Hypothalamus is the weak tissue.** Its clock median (0.4468) barely exceeds
  its target median (0.4379) and the permutation test is not significant
  (*p* = 0.17). It is the one tissue where the effect is absent.
- **Category membership is a curated list, not an ontology.** 213 of 20,955
  genes are classified (16 clock, 22 target, rest across seven other
  categories); everything unclassified is "background", which therefore
  includes genuine clock-controlled genes. This makes the background a
  conservative comparator, not a clean negative set.
- **|λ| depends on sampling interval.** All numbers here are for 2-hourly
  sampling; they are not comparable to values from a differently sampled
  series without rescaling.
- **The φ results are negative and the null is a modelling choice.** The null
  in `fibonacci_null.py` draws (β₁, β₂) uniformly from [−2, 2] × [−1, 1] and
  keeps stationary draws. A null matched to the empirical coefficient
  distribution would be more stringent still.

## Agreement with the published numbers

Compared against the platform's TypeScript outputs for GSE54650:

- The headline liver values agree **exactly**: clock 0.6467, target 0.5285,
  background 0.4963.
- Gene counts agree in every category and tissue.
- Across all 108 category × tissue medians, 81 agree to four decimal places and
  100 agree to within 0.01. The largest disagreement is 0.0371 (lung,
  chromatin); the others above 0.02 are heart/metabolic, white fat/housekeeping,
  brown fat/housekeeping, kidney/housekeeping and brainstem/target.
- The likely cause is probe-to-symbol collapse: the platform and this code
  break ties between multiple probes for one gene differently, which moves a
  handful of genes in and out of small categories and shifts their median. No
  claim in the paper depends on the affected values.

Nothing in this repository has been tuned to reproduce the platform's output.

## Citation

See [`CITATION.cff`](../CITATION.cff). Please cite both the paper and the
`par2-circadian` package.
