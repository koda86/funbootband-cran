# funbootband

<!-- badges: start -->
[![R-CMD-check](https://github.com/koda86/funbootband-cran/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/koda86/funbootband-cran/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/koda86/funbootband-cran/graph/badge.svg)](https://app.codecov.io/gh/koda86/funbootband-cran)
<!-- badges: end -->

`funbootband` computes simultaneous prediction and confidence bands for dense
functional data observed on a common grid. Curves are represented by finite
Fourier series before bootstrap calibration. The package supports independent
curves and repeated curves nested within subjects.

For clustered data, version 0.3.0 uses an intact-subject bootstrap: subjects are
sampled with replacement and all curves belonging to a selected subject are
retained. Subjects receive equal weight, including when cluster sizes differ.
The clustered prediction target is one Fourier-reconstructed future curve from
an independent new subject.

## Installation

Install the CRAN release with:

```r
install.packages("funbootband")
```

Install the development version with:

```r
# install.packages("pak")
pak::pak("koda86/funbootband-cran")
```

## Independent curves

```r
library(funbootband)

set.seed(1)
T <- 101L
n <- 30L
x <- seq(0, 1, length.out = T)
mu <- 0.7 * sin(2 * pi * x) - 0.2 * cos(4 * pi * x)

Y <- replicate(n, {
  mu + rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.30) * sin(2 * pi * x) +
    rnorm(1, sd = 0.20) * cos(2 * pi * x)
})

fit_pred <- band(Y, type = "prediction", alpha = 0.10,
                 iid = TRUE, B = 1000L, k.coef = 4L)
fit_conf <- band(Y, type = "confidence", alpha = 0.10,
                 iid = TRUE, B = 1000L, k.coef = 4L)
```

## Repeated curves nested within subjects

```r
set.seed(2)
K_subject <- 12L
m <- rep(c(2L, 3L, 4L), length.out = K_subject)
id <- rep(seq_len(K_subject), m)

subject_effect <- sapply(seq_len(K_subject), function(i) {
  rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.30) * sin(2 * pi * x) +
    rnorm(1, sd = 0.20) * cos(2 * pi * x)
})

Y_clustered <- sapply(seq_along(id), function(j) {
  mu + subject_effect[, id[j]] +
    rnorm(1, sd = 0.18) * sin(4 * pi * x) +
    rnorm(1, sd = 0.12) * cos(4 * pi * x)
})

fit_clustered <- band(
  Y_clustered,
  type = "prediction",
  alpha = 0.10,
  iid = FALSE,
  id = id,
  B = 1000L,
  k.coef = 4L
)

fit_clustered$meta[c(
  "target", "weighting", "bootstrap_unit", "n_clusters"
)]
```

The clustered band is marginal over the subject population. It is not a
conditional band for an already observed subject and does not provide joint
coverage for several future curves.

## References

- Lenhoff et al. (1999) <doi:10.1016/S0966-6362(98)00043-5>
- Koska et al. (2023) <doi:10.1016/j.jbiomech.2023.111506>
