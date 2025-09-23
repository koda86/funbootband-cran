
# funbootband

<!-- badges: start -->
[![R-CMD-check](https://github.com/koda86/funbootband-cran/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/koda86/funbootband-cran/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/koda86/funbootband-cran/graph/badge.svg)](https://app.codecov.io/gh/koda86/funbootband-cran)
<!-- badges: end -->

`funbootband` computes **simultaneous prediction and confidence bands** for dense functional data (e.g., gait curves sampled on a common grid).  
It supports i.i.d. and **clustered** (hierarchical) designs, uses a **finite Fourier** preprocessing step to honor smoothness/periodicity, and a fast **Rcpp** backend for bootstrap calibration.

## Installation

You can install the development version of funbootband from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("koda86/funbootband-cran")
```

## Example

### i.i.d. example with shaded band

This is a basic example which shows you how to solve a common problem:

``` r
library(funbootband)

set.seed(1)
T <- 60; n <- 8
Y <- matrix(rnorm(T * n, sd = 0.25), nrow = T, ncol = n) +
     outer(seq_len(T), rep(1, n), function(i, j) 0.5 * sin(2*pi*i/T))
fit <- band(Y, type = "prediction", alpha = 0.1, iid = TRUE, B = 25L, k.coef = 12L)
x <- seq_len(fit$meta$T)
plot(x, fit$mean, type = "n",
     ylim = range(c(fit$lower, fit$upper)), xlab = "Index", ylab = "Value")
polygon(c(x, rev(x)), c(fit$lower, rev(fit$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.3), border = NA)
lines(x, fit$mean, lwd = 2)
```

### Clustered example with shaded band

This is a basic example which shows you how to solve a common problem:

``` r
library(funbootband)

set.seed(2)
T  <- 80; m <- c(4, 4)                 # two clusters, few curves
t  <- seq(0, 1, length.out = T)
mu <- list(function(x) 0.7 * sin(2*pi*x),
           function(x) 0.6 * cos(2*pi*x))
Bm <- cbind(sin(2*pi*t), cos(2*pi*t))
gen_curve <- function(k) {
  sc <- rnorm(ncol(Bm), sd = c(0.2, 0.15))
  mu[[k]](t) + as.vector(Bm %*% sc) + rnorm(T, sd = 0.12)
}
Ylist <- lapply(seq_along(m), function(k) sapply(seq_len(m[k]), function(i) gen_curve(k)))
Yh    <- do.call(cbind, Ylist)
id    <- rep(seq_along(m), times = m)
fitH  <- band(Yh, type = "prediction", alpha = 0.1,
              iid = FALSE, id = id, B = 25L, k.coef = 12L)
xh <- seq_len(fitH$meta$T)
plot(xh, fitH$mean, type = "n",
     ylim = range(c(Yh, fitH$lower, fitH$upper), finite = TRUE),
     xlab = "Index", ylab = "Value")
polygon(c(xh, rev(xh)), c(fitH$lower, rev(fitH$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.30), border = NA)
lines(xh, fitH$mean,  lwd = 2)
```

