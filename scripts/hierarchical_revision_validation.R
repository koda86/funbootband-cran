#!/usr/bin/env Rscript

# funbootband: hierarchical-bootstrap revision and validation
# ============================================================
#
# PURPOSE
# This standalone script documents, in executable R, the two urgent changes in
# funbootband 0.3.0. It is deliberately written in transparent base R;
# the package uses equivalent Rcpp kernels for speed.
#
# WHAT CHANGED
# 1. OLD: repeatedly pick a subject, then pick only one of that subject's
#    curves, using curves without replacement until the subject is exhausted.
#    NEW: sample K subjects with replacement and copy ALL observed curves of
#    every selected subject. This is an intact-subject cluster bootstrap.
#
# 2. OLD: the clustered prediction target was not explicitly defined, pooled
#    curve means implicitly over-weighted subjects with more curves, and the
#    calibration maximized over every observed curve at once.
#    NEW: the target is ONE future curve from ONE new subject. First draw a
#    subject from the subject population, then draw a curve conditionally from
#    that subject. Subjects are equally weighted; curves are equally weighted
#    within subject. Calibration retains one curve-level supremum statistic.
#
# 3. NEW SUPPORTING FIXES:
#    - use the bootstrap-replicate-specific pointwise prediction scale;
#    - base clustered confidence standard errors on subject mean curves;
#    - expose target, weighting, and bootstrap unit in fit$meta.
#
# TARGET AND DERIVATION
# Let Y_ij(t) = mu(t) + A_i(t) + E_ij(t), where subjects i = 1,...,K are
# independent and subject i supplies m_i repeated curves. Here Y denotes the
# representation actually calibrated by the algorithm: in the package this is
# the finite-Fourier reconstruction, not necessarily the raw noisy measurement.
# The intended future quantity is
#
#   Y_new(t) = mu(t) + A_new(t) + E_new(t),
#
# for an independent new subject and one new curve from that subject. The
# subject-first empirical distribution is
#
#   P_hat = (1/K) sum_i (1/m_i) sum_j delta_{Y_ij}.
#
# Hence observed curve (i,j) has future weight q_ij = 1/(K*m_i), and
#
#   mu_hat(t)     = sum_ij q_ij Y_ij(t),
#   sigma_hat^2(t)= sum_ij q_ij {Y_ij(t)-mu_hat(t)}^2.
#
# In bootstrap replicate b, draw K subject labels with replacement. If subject
# i occurs c_bi times, every curve from that subject receives training weight
#
#   w_bij = c_bi/(K*m_i).
#
# Compute mu_b and sigma_b from w_bij. For every empirical pseudo-future curve,
# retain its curve-level simultaneous statistic
#
#   M_bij = max_t |Y_ij(t)-mu_b(t)| / sigma_b(t).
#
# The critical value is the (1-alpha) weighted empirical quantile of all M_bij,
# using q_ij/B. The final band is mu_hat(t) +/- c * sigma_hat(t).
#
# This estimates the new-subject/new-curve marginal prediction target. It does
# NOT target another curve from an already observed subject, a subject-specific
# conditional band, or joint coverage of several future curves. Those are later
# extensions and require different resampling/prediction constructions.

.cluster_curve_weights_reference <- function(id) {
  id <- as.integer(as.factor(id))
  sizes <- tabulate(id)
  1 / (length(sizes) * sizes[id])
}

.intact_subject_boot_weights_reference <- function(id, B) {
  id <- as.integer(as.factor(id))
  groups <- split(seq_along(id), id)
  K <- length(groups)
  n <- length(id)
  out <- matrix(0, nrow = B, ncol = n)

  selected <- matrix(sample.int(K, B * K, replace = TRUE), nrow = B)
  for (b in seq_len(B)) {
    counts <- tabulate(selected[b, ], nbins = K)
    for (i in which(counts > 0L)) {
      # All curves are retained; there is no within-subject second-stage draw.
      out[b, groups[[i]]] <- counts[i] / (K * length(groups[[i]]))
    }
  }
  out
}

.weighted_ecdf_quantile_reference <- function(x, w, probability) {
  stopifnot(length(x) == length(w), probability > 0, probability < 1)
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  x[which(cumsum(w) >= probability)[1L]]
}

prediction_band_new_subject_reference <- function(data, id, alpha = 0.10,
                                                  B = 499L,
                                                  ridge = 1e-12) {
  data <- as.matrix(data)
  id <- as.integer(as.factor(id))
  stopifnot(is.numeric(data), ncol(data) == length(id),
            length(unique(id)) >= 2L, B >= 2L)

  future_weights <- .cluster_curve_weights_reference(id)
  mu_hat <- as.numeric(data %*% future_weights)
  residuals <- sweep(data, 1L, mu_hat, FUN = "-")
  sigma_hat <- sqrt(as.numeric((residuals^2) %*% future_weights))
  sigma_hat <- pmax(sigma_hat, ridge)

  boot_weights <- .intact_subject_boot_weights_reference(id, B)
  B <- nrow(boot_weights)
  n <- ncol(data)
  curve_maxima <- matrix(NA_real_, nrow = B, ncol = n)

  for (b in seq_len(B)) {
    w <- boot_weights[b, ]
    mu_b <- as.numeric(data %*% w)
    residuals_b <- sweep(data, 1L, mu_b, FUN = "-")
    sigma_b <- sqrt(as.numeric((residuals_b^2) %*% w))
    sigma_b <- pmax(sigma_b, ridge)

    # One maximum per curve: deliberately no max over the curve dimension.
    standardized <- abs(sweep(residuals_b, 1L, sigma_b, FUN = "/"))
    curve_maxima[b, ] <- apply(standardized, 2L, max)
  }

  critical_value <- .weighted_ecdf_quantile_reference(
    as.vector(curve_maxima),
    rep(future_weights, each = B),
    1 - alpha
  )

  list(
    lower = mu_hat - critical_value * sigma_hat,
    mean = mu_hat,
    upper = mu_hat + critical_value * sigma_hat,
    critical_value = critical_value,
    meta = list(
      target = "new_subject_new_curve",
      weighting = "equal_subject_then_equal_curve_within_subject",
      bootstrap_unit = "intact_subject",
      curve_representation = "input_to_reference_function",
      B = B,
      alpha = alpha
    )
  )
}

run_structural_validation <- function() {
  message("1/3 Checking intact-subject bootstrap weights ...")
  id <- c(1, 1, 2, 2, 2, 3, 3, 3, 3)
  set.seed(20260907)
  W <- .intact_subject_boot_weights_reference(id, B = 200L)
  K <- length(unique(id))

  stopifnot(max(abs(rowSums(W) - 1)) < 1e-12)
  for (i in unique(id)) {
    # A selected subject enters with all curves at the same within-subject
    # weight; an unselected subject has zero weight for every one of its curves.
    stopifnot(all(apply(W[, id == i, drop = FALSE], 1L,
                        function(z) length(unique(z)) == 1L)))
    subject_totals <- rowSums(W[, id == i, drop = FALSE])
    stopifnot(max(abs(K * subject_totals - round(K * subject_totals))) < 1e-12)
  }

  message("2/3 Checking the subject-weighted target with unequal cluster sizes ...")
  Tlen <- 11L
  id2 <- c(1, 1, 2, 2, 2, 2)
  Y <- matrix(rep(c(0, 0, 12, 12, 12, 12), each = Tlen), nrow = Tlen)
  fit <- prediction_band_new_subject_reference(Y, id2, B = 99L)
  # Equal subjects => (0 + 12)/2 = 6. A pooled curve mean would incorrectly be 8.
  stopifnot(max(abs(fit$mean - 6)) < 1e-12)

  message("3/3 Checking package metadata when the revised package is installed ...")
  if (requireNamespace("funbootband", quietly = TRUE) &&
      utils::packageVersion("funbootband") >= "0.3.0") {
    set.seed(20260907)
    package_fit <- funbootband::band(
      Y, type = "prediction", iid = FALSE, id = id2,
      alpha = 0.10, B = 99L, k.coef = 0L
    )
    stopifnot(identical(package_fit$meta$target, "new_subject_new_curve"),
              identical(package_fit$meta$bootstrap_unit, "intact_subject"),
              max(abs(package_fit$mean - 6)) < 1e-10)
  } else {
    message("    Revised package not installed; skipped package-level comparison.")
  }

  message("Structural validation passed.")
  invisible(TRUE)
}

.simulate_hierarchical_dataset <- function(K = 25L, m = 4L, Tlen = 41L) {
  x <- seq(0, 1, length.out = Tlen)
  mu <- 0.5 * sin(2 * pi * x)
  curves <- vector("list", K * m)
  id <- rep(seq_len(K), each = m)
  cursor <- 1L

  for (i in seq_len(K)) {
    subject_effect <-
      rnorm(1, sd = 0.35) +
      rnorm(1, sd = 0.25) * sin(2 * pi * x) +
      rnorm(1, sd = 0.15) * cos(4 * pi * x)
    for (j in seq_len(m)) {
      curve_effect <-
        rnorm(1, sd = 0.22) * cos(2 * pi * x) +
        rnorm(1, sd = 0.12) * sin(4 * pi * x) +
        rnorm(Tlen, sd = 0.05)
      curves[[cursor]] <- mu + subject_effect + curve_effect
      cursor <- cursor + 1L
    }
  }

  # Independent subject effect plus independent within-subject curve effect:
  # exactly the new-subject/new-curve target used in the derivation above.
  new_subject_effect <-
    rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.25) * sin(2 * pi * x) +
    rnorm(1, sd = 0.15) * cos(4 * pi * x)
  new_curve_effect <-
    rnorm(1, sd = 0.22) * cos(2 * pi * x) +
    rnorm(1, sd = 0.12) * sin(4 * pi * x) +
    rnorm(Tlen, sd = 0.05)

  list(
    data = do.call(cbind, curves),
    id = id,
    future = mu + new_subject_effect + new_curve_effect
  )
}

.wilson_interval <- function(successes, trials, level = 0.95) {
  p <- successes / trials
  z <- stats::qnorm(1 - (1 - level) / 2)
  denominator <- 1 + z^2 / trials
  centre <- (p + z^2 / (2 * trials)) / denominator
  half <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) / denominator
  c(lower = centre - half, estimate = p, upper = centre + half)
}

run_monte_carlo_validation <- function(full = FALSE) {
  # Quick mode is an implementation smoke test. For manuscript evidence, run:
  #   Sys.setenv(FUNBOOTBAND_FULL_VALIDATION = "true")
  #   source("scripts/hierarchical_revision_validation.R")
  outer_replicates <- if (full) 300L else 50L
  bootstrap_replicates <- if (full) 999L else 199L
  nominal <- 0.90
  covered <- logical(outer_replicates)
  mean_width <- numeric(outer_replicates)

  set.seed(20260907)
  message(sprintf(
    "Monte Carlo validation: R=%d, B=%d, nominal=%.2f ...",
    outer_replicates, bootstrap_replicates, nominal
  ))
  for (r in seq_len(outer_replicates)) {
    sample <- .simulate_hierarchical_dataset()
    fit <- prediction_band_new_subject_reference(
      sample$data, sample$id, alpha = 1 - nominal,
      B = bootstrap_replicates
    )
    covered[r] <- all(sample$future >= fit$lower & sample$future <= fit$upper)
    mean_width[r] <- mean(fit$upper - fit$lower)
  }

  interval <- .wilson_interval(sum(covered), outer_replicates)
  result <- data.frame(
    target = "new_subject_new_curve",
    nominal_coverage = nominal,
    empirical_coverage = unname(interval["estimate"]),
    wilson_95_lower = unname(interval["lower"]),
    wilson_95_upper = unname(interval["upper"]),
    mean_band_width = mean(mean_width),
    outer_replicates = outer_replicates,
    bootstrap_replicates = bootstrap_replicates
  )
  print(result, row.names = FALSE)
  message(
    "Interpretation: this checks the implemented estimand under one model; ",
    "it is not a general proof of validity or a substitute for the manuscript study."
  )
  invisible(result)
}

run_structural_validation()
full_validation <- identical(
  tolower(Sys.getenv("FUNBOOTBAND_FULL_VALIDATION", "false")), "true"
)
validation_result <- run_monte_carlo_validation(full = full_validation)

# PACKAGE CHECKS TO RUN FROM THE REVISED PACKAGE ROOT
# ---------------------------------------------------
# install.packages(c("Rcpp", "roxygen2", "testthat", "devtools"))
# Rcpp::compileAttributes()
# roxygen2::roxygenise()
# devtools::test()
# devtools::check(args = "--as-cran")
#
# The package-level tests additionally verify:
# - unequal cluster sizes produce the subject-weighted rather than pooled mean;
# - bootstrap weight rows retain whole subjects and sum to one;
# - the prediction kernel returns a B x n matrix of curve-level maxima;
# - fit$meta reports the prediction target and intact-subject sampling unit.
