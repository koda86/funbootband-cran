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

## R CMD check results

* Win-builder R-devel: 0 errors, 0 warnings, 0 notes.
* Win-builder R-release: 0 errors, 0 warnings, 1 note.

The R-release NOTE resulted from timeouts while Win-builder was checking
valid GitHub URLs listed in DESCRIPTION and README.md. The same source
package passed these checks on R-devel. This was a transient external
network issue.

## Reverse dependencies

No reverse dependencies were listed on the CRAN package page when preparing
this release candidate. Recheck immediately before submission.
