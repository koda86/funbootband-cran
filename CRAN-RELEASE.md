# CRAN release checklist for funbootband 0.3.0

Run all commands from a clean checkout containing only the intended release
changes.

## 1. Regenerate generated files

```r
Rcpp::compileAttributes()
roxygen2::roxygenise()
```

Review the generated diff, especially `R/RcppExports.R`,
`src/RcppExports.cpp`, and `man/band.Rd`.

## 2. Automated and statistical tests

```r
devtools::test()
devtools::check(args = "--as-cran")
```

Run the installed validation study:

```r
Sys.setenv(FUNBOOTBAND_VALIDATION_R = 500,
           FUNBOOTBAND_VALIDATION_B = 999)
source(system.file("validation", "coverage_simulation.R",
                   package = "funbootband"))
```

Retain the output with the release records. Interpret empirical coverage using
its Monte Carlo standard error; investigate material undercoverage before
submission.

## 3. Build and check the exact tarball

```sh
R CMD build .
R CMD check --as-cran funbootband_0.3.0.tar.gz
```

Do not submit the manually assembled development archive. Submit only the
tarball produced by `R CMD build` and checked in this step.

## 4. External checks

* Run Winbuilder with R-release and R-devel.
* Run R-hub on at least Linux, Windows, and macOS where available.
* Check the current CRAN results page for version 0.2.0.
* Recheck reverse dependencies.
* Verify all URLs and the maintainer email.

## 5. Release metadata

* Confirm `DESCRIPTION` reports version 0.3.0.
* Confirm `NEWS.md` describes the changed statistical behaviour.
* Replace the placeholder check-results paragraph in `cran-comments.md`.
* Commit and tag the exact source used to build the submitted tarball.
* Submit through the CRAN web form and confirm the submission email.
