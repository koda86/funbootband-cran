#!/usr/bin/env Rscript

# Coverage study for funbootband 0.3.0
# =====================================
#
# This script is installed for reproducibility but is not run during CRAN
# checks. It evaluates the two clustered targets implemented by band():
#   1. one Fourier-reconstructed curve from a new subject;
#   2. the equally subject-weighted population mean curve.
#
# Quick implementation check:
#   source(system.file("validation", "coverage_simulation.R",
#                      package = "funbootband"))
#
# More precise run before publication or release:
#   Sys.setenv(FUNBOOTBAND_VALIDATION_R = 500,
#              FUNBOOTBAND_VALIDATION_B = 999)
#   source(system.file("validation", "coverage_simulation.R",
#                      package = "funbootband"))

if (!requireNamespace("funbootband", quietly = TRUE)) {
  stop("Install funbootband before running this validation study.")
}

outer_R <- as.integer(Sys.getenv("FUNBOOTBAND_VALIDATION_R", "50"))
bootstrap_B <- as.integer(Sys.getenv("FUNBOOTBAND_VALIDATION_B", "199"))
alpha <- 0.10
Tlen <- 61L
x <- seq(0, 1, length.out = Tlen)
mu_true <- 0.7 * sin(2 * pi * x) - 0.2 * cos(4 * pi * x)

scenarios <- data.frame(
  scenario = c("K15_balanced", "K15_unbalanced", "K30_unbalanced"),
  K = c(15L, 15L, 30L),
  size_pattern = c("balanced", "unbalanced", "unbalanced"),
  stringsAsFactors = FALSE
)

simulate_data <- function(K, size_pattern) {
  m <- if (identical(size_pattern, "balanced")) {
    rep(3L, K)
  } else {
    rep(c(2L, 3L, 5L), length.out = K)
  }
  id <- rep(seq_len(K), m)

  subject_effect <- sapply(seq_len(K), function(i) {
    rnorm(1, sd = 0.35) +
      rnorm(1, sd = 0.30) * sin(2 * pi * x) +
      rnorm(1, sd = 0.20) * cos(2 * pi * x)
  })

  within_effect <- function() {
    rnorm(1, sd = 0.18) * sin(4 * pi * x) +
      rnorm(1, sd = 0.12) * cos(4 * pi * x)
  }

  Y <- sapply(seq_along(id), function(j) {
    mu_true + subject_effect[, id[j]] + within_effect()
  })

  new_subject_effect <-
    rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.30) * sin(2 * pi * x) +
    rnorm(1, sd = 0.20) * cos(2 * pi * x)
  future_curve <- mu_true + new_subject_effect + within_effect()

  list(data = Y, id = id, future = future_curve)
}

set.seed(20260917)
results <- vector("list", nrow(scenarios))

for (s in seq_len(nrow(scenarios))) {
  prediction_covered <- logical(outer_R)
  confidence_covered <- logical(outer_R)
  prediction_width <- numeric(outer_R)
  confidence_width <- numeric(outer_R)

  for (r in seq_len(outer_R)) {
    sample <- simulate_data(scenarios$K[s], scenarios$size_pattern[s])

    fit_prediction <- funbootband::band(
      sample$data, type = "prediction", alpha = alpha,
      iid = FALSE, id = sample$id, B = bootstrap_B, k.coef = 4L
    )
    fit_confidence <- funbootband::band(
      sample$data, type = "confidence", alpha = alpha,
      iid = FALSE, id = sample$id, B = bootstrap_B, k.coef = 4L
    )

    prediction_covered[r] <- all(
      sample$future >= fit_prediction$lower &
        sample$future <= fit_prediction$upper
    )
    confidence_covered[r] <- all(
      mu_true >= fit_confidence$lower & mu_true <= fit_confidence$upper
    )
    prediction_width[r] <- mean(fit_prediction$upper - fit_prediction$lower)
    confidence_width[r] <- mean(fit_confidence$upper - fit_confidence$lower)
  }

  results[[s]] <- data.frame(
    scenario = scenarios$scenario[s],
    target = c("new_subject_new_curve", "subject_weighted_population_mean"),
    nominal_coverage = 1 - alpha,
    empirical_coverage = c(mean(prediction_covered),
                           mean(confidence_covered)),
    monte_carlo_se = c(
      sqrt(mean(prediction_covered) * (1 - mean(prediction_covered)) / outer_R),
      sqrt(mean(confidence_covered) * (1 - mean(confidence_covered)) / outer_R)
    ),
    mean_band_width = c(mean(prediction_width), mean(confidence_width)),
    outer_replicates = outer_R,
    bootstrap_replicates = bootstrap_B
  )
}

coverage_results <- do.call(rbind, results)
row.names(coverage_results) <- NULL
print(coverage_results)

message(
  "Interpret coverage together with its Monte Carlo standard error. ",
  "This study validates the stated data-generating scenarios, not every ",
  "possible functional-data distribution."
)
