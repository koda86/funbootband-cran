## Submission type

Update from CRAN version 0.2.0 to 0.3.0.

## Major changes

* Corrected clustered resampling so that subjects are sampled with replacement
  and all curves belonging to each selected subject are retained intact.
* Defined clustered prediction as one Fourier-reconstructed future curve from
  an independent new subject.
* Introduced equal-subject/equal-within-subject weighting, including for unequal
  cluster sizes.
* Corrected prediction calibration to retain curve-level supremum statistics
  and use bootstrap-replicate-specific pointwise scales.
* Clustered confidence bands now use subject mean curves for their pointwise
  standard errors and replicate-specific studentization.
* Updated documentation, examples, tests, and validation materials.

## R CMD check results

Before submission, replace this paragraph with the results from checking the
exact 0.3.0 tarball on local R-release/R-devel and Winbuilder. Explain every
remaining NOTE. The intended submission standard is 0 errors, 0 warnings, and
no unexplained notes.

## Reverse dependencies

No reverse dependencies were listed on the CRAN package page when preparing
this release candidate. Recheck immediately before submission.
