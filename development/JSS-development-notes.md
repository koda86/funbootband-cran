# `funbootband`: JSS Development and Evidence Record

Last updated: 2026-09-18

## Purpose

This file is the working source of truth for development of `funbootband` after
the CRAN 0.3.0 release and for preparation of a Journal of Statistical Software
(JSS) manuscript. It records methodological decisions, implemented behavior,
validation already completed, claims that may safely be made, and work that is
still planned.

This is a development document, not package documentation. Keep it out of the
built CRAN source package by listing `^development$` in `.Rbuildignore`.

## Current release baseline

- CRAN version: 0.3.0.
- Release date: 2026-09-18.
- Git tag: `v0.3.0`.
- Post-release development branch: `develop/jss`.
- Planned development version: 0.3.0.9000.

Version 0.3.0 is the methodological baseline for all subsequent work. Later
changes must be documented here, in `NEWS.md`, in user documentation, and in
tests when they affect observable behavior.

## Statistical scope and current decisions

### Data representation and scope

The package is intended for dense functional observations recorded on a common
grid. The current preprocessing represents each observed curve by a finite
Fourier series and reconstructs it on the observed grid before bootstrap
calibration.

Consequences that must be stated clearly:

- `k.coef` controls the retained Fourier harmonics and therefore the amount of
  high-frequency detail in the reconstructed curves.
- The implementation clamps `k.coef` to the maximum meaningful value for the
  grid.
- A Fourier basis is most natural for periodic or approximately periodic data.
  It is not a universally appropriate smoother for arbitrary functional data.
- The current interface assumes a common, ordered grid. Support for an explicit
  grid is planned below.

### Independent-curve design

For `iid = TRUE`, curves are treated as independent observational units.

- A confidence band targets the population mean function.
- A prediction band targets one future curve from the same curve population.
- Calibration is simultaneous over the evaluation grid, rather than a
  collection of separately calibrated pointwise intervals.

### Clustered design: estimand and weighting

For `iid = FALSE`, the observational hierarchy is subjects (clusters) with one
or more curves per subject.

The current clustered confidence target is the subject-weighted population
mean function. Each subject receives equal weight; within a subject, its curves
receive equal weight. Thus subjects with more recorded curves do not
automatically receive greater population weight.

The current and default clustered prediction target is:

> one future Fourier-reconstructed curve from an independent new subject drawn
> from the target subject population.

This is the **new-subject/new-curve** target. It includes both between-subject
and within-subject curve variation under the population and sampling model
represented by the observed clustered data.

Other prediction targets, such as a new curve from an already observed subject
or a new subject mean curve, are not currently implemented. They are possible
future extensions and must not be implied by the present API.

### Clustered bootstrap

The current default is an intact-subject cluster bootstrap:

1. Sample subjects with replacement.
2. For every selected subject, retain that subject's complete observed set of
   curves as an intact cluster.
3. Apply equal-subject/equal-within-subject weighting in constructing the
   relevant estimates and calibration quantities.

There is no additional within-subject resampling stage in the current default.
This choice preserves observed within-subject dependence and aligns the
resampling unit with the new-subject inferential target.

The package must not be described as exactly implementing the two-stage
hierarchical procedure in Koska et al. (2023). That article motivates the
hierarchical application, but version 0.3.0 deliberately uses a revised,
intact-subject procedure.

### Band calibration

The current prediction calibration retains a curve-level supremum statistic
for each pseudo-future curve. It does not replace those quantities by one
maximum over the entire observed collection. The implementation uses
bootstrap-replicate-specific pointwise scales.

For clustered confidence bands, subject mean curves are the cluster-level
units used to form pointwise standard errors. The bootstrap resamples intact
subjects and uses replicate-specific studentization.

### Interpretation of coverage

The intended coverage is simultaneous across the evaluation grid. A statement
such as “90% prediction band” refers to the probability that an entire future
curve lies within the band over that grid under the stated target and data-
generating assumptions.

Do not describe nominal coverage as a universal finite-sample guarantee. The
procedure is bootstrap-calibrated, and its actual coverage depends on the data
structure, sample size, basis representation, grid, bootstrap size, and the
adequacy of the resampling assumptions. Empirical validation supports only the
specific scenarios studied.

## Evidence and validation completed for 0.3.0

- Unit tests: 53 passed, 0 failed, and 0 warnings in the final local test run.
- Local package checks completed without errors or warnings. A local
  compiler-configuration note about `-mno-omit-leaf-frame-pointer` was external
  to the package source.
- Win-builder R-devel completed with status OK.
- A Win-builder R-release run produced one incoming-feasibility NOTE caused by
  network timeouts while checking GitHub URLs; the same source subsequently
  passed the substantive checks.
- The coverage-validation script used 500 outer Monte Carlo replicates and 999
  bootstrap replicates in balanced and unbalanced clustered scenarios.
- Those simulations examined the new-subject/new-curve prediction target and
  the subject-weighted population-mean confidence target.

Limits of this evidence:

- The simulations validate the stated data-generating scenarios, not all
  possible functional-data distributions.
- The existing simulation grid is not yet a full comparative study.
- Runtime and scaling have not yet been characterized systematically.
- Sensitivity to `k.coef`, bootstrap replicate count, nonperiodicity, grid
  resolution, small subject counts, and stronger within-subject dependence
  needs broader evaluation for the JSS paper.

## Claims discipline for documentation and manuscript

### Claims currently supportable

- `funbootband` computes bootstrap-calibrated simultaneous confidence and
  prediction bands for densely sampled functional data on a common grid.
- It supports independent curves and clustered curves.
- Its clustered default uses intact-subject resampling, equal subject weighting,
  and a new-subject/new-curve prediction target.
- Curves are represented through finite Fourier series before calibration.
- Computational kernels are implemented with Rcpp.

### Claims requiring qualification or more evidence

- “Correct coverage” must be tied to assumptions and empirical scenarios.
- “Fast” or “computationally efficient” requires reproducible benchmarks.
- “Applicable to functional data” must be narrowed to the supported dense,
  common-grid setting and must discuss the Fourier representation.
- “Smoothness” is induced by basis truncation; it should not be presented as an
  automatically validated property of every data set.
- “Periodicity” is a modeling consideration, not a fact inferred by the
  package.
- Comparisons with alternative software require a current, source-supported
  related-software review.

### Claims to avoid

- The package exactly reproduces the hierarchical algorithm in Koska et al.
  (2023).
- The clustered procedure supports every possible future-observation target.
- Nominal simultaneous coverage is guaranteed for arbitrary samples or
  distributions.
- The default value of `k.coef` is optimal without data-specific assessment.

## Planned package-development sequence

### 1. Development infrastructure

- Set the version to 0.3.0.9000.
- Keep this record under `development/` and exclude that directory from builds.
- Record every user-visible change in `NEWS.md`.
- Work in focused branches and require clean tests and checks before merging.

### 2. Stable returned object and S3 interface

- Give results a documented S3 class while retaining the existing list
  components needed for backward compatibility.
- Implement `print()` with the band type, design, nominal level, dimensions,
  and target.
- Implement `summary()` with settings, target/estimand, cluster information,
  weighting, Fourier order, bootstrap count, and useful band-width summaries.
- Implement `plot()` for the band and estimated mean, with optional observed
  curves and suitable labels.
- Add method registration, documentation, examples, and tests.

Design decisions still needed:

- Final class name (`funbootband` versus `funbootband_band`).
- Plot defaults and whether base graphics alone remain sufficient.
- Which diagnostics belong in `meta` and which should be computed by
  `summary()`.

### 3. Explicit grid support

- Add a documented grid argument while preserving current behavior when it is
  omitted.
- Validate length, finiteness, ordering, duplicates, and spacing assumptions.
- Decide whether version 0.3.x supports only equally spaced grids or also a
  well-defined transformation for irregular grids.
- Use the actual grid for plots and return it in the fitted object.
- Explain that Fourier reconstruction presumes an appropriate periodic domain.
- Add unit tests and vignette examples.

### 4. Input and metadata improvements

- Encourage explicit `id` values for clustered data; keep column-name inference
  only if it can be explained and tested unambiguously.
- Make the inferential target, weighting convention, bootstrap unit, effective
  Fourier order, grid, cluster sizes, and randomization settings readily
  inspectable.
- Review validation messages and edge cases, especially one-subject designs,
  singleton clusters, missing values, and degenerate pointwise variability.

### 5. Documentation and teaching examples

- Revise the reference page and vignette after the API is stable.
- Include separate independent and clustered examples.
- Explain confidence versus prediction targets before presenting code.
- Include an unequal-cluster example that makes subject weighting visible.
- Explain selection and sensitivity of `B` and `k.coef` without implying a
  universally optimal choice.
- Add an applied biomechanical example if licensing and reproducibility permit.

### 6. JSS replication and evaluation suite

- Create a standalone replication entry point for every manuscript table and
  figure.
- Record package and platform information and use reproducible seeds.
- Separate quick checks from full computational experiments.
- Add coverage experiments across subject counts, cluster-size balance,
  dependence strength, curve distributions, grid sizes, and Fourier mismatch.
- Add sensitivity analyses for `B` and `k.coef`.
- Add runtime and memory benchmarks over representative dimensions.
- Compare only with methods and packages that target genuinely comparable
  inferential objects; clearly label pointwise or mean-only methods as such.

## Provisional JSS manuscript structure

1. Introduction and statistical problem
2. Inferential targets and assumptions
   1. Independent curves
   2. Clustered curves and subject weighting
   3. Confidence versus prediction bands
3. Methodology
   1. Fourier representation
   2. Simultaneous calibration
   3. Intact-subject bootstrap
4. Software design and implementation
   1. User interface and returned object
   2. Rcpp implementation
   3. Reproducibility and input validation
5. Using `funbootband`
   1. Independent example
   2. Clustered example
   3. Interpretation and diagnostics
6. Empirical evaluation
   1. Coverage study
   2. Sensitivity analysis
   3. Computational performance
7. Applied biomechanical example
8. Related software
9. Discussion, limitations, and future work

## JSS completion requirements to track

- Use the current official JSS manuscript style and author instructions.
- Supply source for the manuscript and software.
- Supply a standalone replication script or reproducibility bundle for
  empirical results.
- Explain the statistical methodology sufficiently for readers to understand
  the implemented targets and assumptions.
- Explain implementation and API design, not only the motivating statistics.
- Compare with related implementations by scope and target rather than by a
  superficial feature list.
- Ensure all manuscript examples run against the archived package version used
  for publication.
- Freeze the public API before finalizing manuscript code and screenshots.

## Open evidence tasks

- Re-audit the current CRAN, Bioconductor, GitHub, JSS, and scholarly software
  landscape before writing the related-software section.
- Verify every description of Lenhoff et al. (1999) and Koska et al. (2023)
  against the primary texts.
- Document the mathematical definition of every estimand and bootstrap
  statistic in notation consistent with the code.
- Expand simulation scenarios and pre-specify reported outcomes.
- Benchmark the Rcpp implementation against a transparent reference
  implementation where feasible.
- Select and document a distributable applied data example.
- Decide whether optional future prediction targets are warranted only after
  the default target is fully documented and validated.

## Change log for this record

- 2026-09-18: Created from the decisions and validation work leading to CRAN
  version 0.3.0 and established the post-release JSS development plan.
