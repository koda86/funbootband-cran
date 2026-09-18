# funbootband (development version)

* Began JSS-oriented development following the CRAN 0.3.0 release.

# funbootband 0.3.0

* Replaced the former curve-by-curve hierarchical sampler with an intact-subject
  cluster bootstrap. Subjects are sampled with replacement and all observed
  curves belonging to a selected subject are retained.
* Defined the clustered prediction target as one future curve from a new
  subject. Subjects are weighted equally and curves are weighted equally within
  subject, including for unequal cluster sizes.
* Prediction calibration now retains one supremum statistic per pseudo-future
  curve instead of maximizing over the whole observed collection, and uses a
  replicate-specific pointwise scale.
* Clustered confidence bands now use subject mean curves for their pointwise
  standard errors, resample intact subjects, and use replicate-specific
  studentization.
* The returned metadata now records the estimand, weighting convention,
  bootstrap unit, cluster sizes, and Fourier-reconstructed curve
  representation.
* Documentation and examples now distinguish the revised intact-subject method
  from the hierarchical procedure described in Koska et al. (2023).
* Added numerical reference tests for the Rcpp kernels, unequal-cluster
  examples, edge-case tests, and a reproducible coverage-validation script.

# funbootband 0.2.0

# funbootband 0.1.1 (2025-09-22)

- Add Fourier preprocessing in `band()` via `k.coef` (default 50) to honor smooth/periodic structure.
- New vignette: “funbootband: Simultaneous Bands for Functional Data”.
- Faster, cleaned bootstrap pipeline (Rcpp); small `B` used in examples/tests for CRAN timing.
- Documentation and tests refreshed.

# funbootband 0.1.0 (2025-08-01)

- Initial release: simultaneous prediction and confidence bands with clustered bootstrap.
