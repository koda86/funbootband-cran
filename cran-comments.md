## Test environments

* Local: Ubuntu 24.04.2 LTS, R 4.5.1 (x86_64), GCC 13.3.0
* No external services; no network access required.
* Examples and vignette use small `B` on CRAN to keep run time short.
* win-builder (Windows Server 2022):
  - R-release: 0 errors, 0 warnings, 1 note (“New submission”)
  - R-devel:   0 errors, 0 warnings, 1 note (“New submission”)

## R CMD check results

0 errors | 0 warnings | 3 notes

* **New submission**
  - This is the first CRAN release of the package.

* **Compilation flags**
  - NOTE: “Compilation used the following non-portable flag(s): ‘-mno-omit-leaf-frame-pointer’”.
  - This flag is injected by Ubuntu’s system-wide `Makeconf` and is **not** set in the package sources (no `src/Makevars*`). CRAN’s build machines will not inherit these local flags.

* **HTML manual checks**
  - NOTE: Local environment lacks `tidy` and `V8`, so HTML validation and math rendering checks were skipped.
  - CRAN machines provide these tools; this NOTE is local-only.

## Additional package notes

* No writing to the file system beyond temporary directories.
* No remote resources are accessed in examples, tests, or vignettes.
* The vignette adapts computation via `NOT_CRAN` to meet CRAN run-time expectations.

## Downstream dependencies

* This is an initial submission; there are no reverse dependencies yet.

