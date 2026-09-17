# Tests for funbootband::band()

.prediction_kernel_reference <- function(data, boot_weights, ridge = 1e-12) {
  B <- nrow(boot_weights)
  n <- ncol(data)
  out <- matrix(NA_real_, nrow = B, ncol = n)
  for (b in seq_len(B)) {
    mu_b <- as.numeric(data %*% boot_weights[b, ])
    residuals <- sweep(data, 1L, mu_b, FUN = "-")
    sd_b <- sqrt(as.numeric((residuals^2) %*% boot_weights[b, ]))
    sd_b <- pmax(sd_b, ridge)
    out[b, ] <- apply(abs(sweep(residuals, 1L, sd_b, FUN = "/")), 2L, max)
  }
  out
}

.confidence_kernel_reference <- function(unit_data, mu_hat, boot_weights,
                                         ridge = 1e-12) {
  U <- ncol(unit_data)
  apply(boot_weights, 1L, function(w) {
    mu_b <- as.numeric(unit_data %*% w)
    residuals <- sweep(unit_data, 1L, mu_b, FUN = "-")
    se_b <- sqrt(as.numeric((residuals^2) %*% w) / (U - 1L))
    se_b <- pmax(se_b, ridge)
    max(abs((mu_b - mu_hat) / se_b))
  })
}

test_that("band() runs for valid inputs (i.i.d., prediction & confidence)", {
  set.seed(1)
  T <- 40; n <- 20
  Y <- matrix(rnorm(T * n, sd = 0.3), nrow = T, ncol = n)

  # prediction band
  fitP <- band(Y, type = "prediction", alpha = 0.10, iid = TRUE,
               B = 120, k.coef = 8L)
  expect_type(fitP, "list")
  expect_equal(sort(names(fitP)), sort(c("lower","mean","upper","meta")))
  expect_equal(length(fitP$lower), T)
  expect_equal(length(fitP$upper), T)
  expect_true(all(is.finite(fitP$lower)))
  expect_true(all(is.finite(fitP$upper)))
  expect_true(all(fitP$upper >= fitP$lower))

  # confidence band
  fitC <- band(Y, type = "confidence", alpha = 0.10, iid = TRUE,
               B = 120, k.coef = 8L)
  expect_type(fitC, "list")
  expect_equal(length(fitC$lower), T)
  expect_true(all(fitC$upper >= fitC$lower))
})

test_that("band() works for clustered data via explicit id", {
  set.seed(2)
  T <- 30; n <- 24
  id <- rep(1:6, each = 4)                      # 6 clusters, each size 4
  Y  <- matrix(rnorm(T * n, sd = 0.25), nrow = T, ncol = n)

  fit <- band(Y, type = "prediction", alpha = 0.10, iid = FALSE,
              id = id, B = 100, k.coef = 8L)
  expect_equal(length(fit$lower), T)
  expect_true(all(is.finite(fit$lower)))
  expect_true(all(is.finite(fit$upper)))
  expect_true(all(fit$upper >= fit$lower))
  expect_identical(fit$meta$target, "new_subject_new_curve")
  expect_identical(fit$meta$bootstrap_unit, "intact_subject")
  expect_equal(fit$meta$n_clusters, 6L)
})

test_that("clustered centre is subject-weighted when cluster sizes differ", {
  # Subject 1 has two zero curves; subject 2 has four curves equal to 12.
  # A pooled curve mean is 8, whereas the stated new-subject target has mean 6.
  T <- 8
  id <- c(1, 1, 2, 2, 2, 2)
  Y <- matrix(rep(c(0, 0, 12, 12, 12, 12), each = T), nrow = T)

  set.seed(20)
  fit <- band(Y, type = "prediction", alpha = 0.10, iid = FALSE,
              id = id, B = 50, k.coef = 0)
  fit_conf <- band(Y, type = "confidence", alpha = 0.10, iid = FALSE,
                   id = id, B = 50, k.coef = 0)

  expect_equal(fit$mean, rep(6, T), tolerance = 1e-10)
  expect_equal(fit_conf$mean, rep(6, T), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fit$mean, rep(rowMeans(Y)[1], T))))
  expect_identical(fit$meta$weighting,
                   "equal_subject_then_equal_curve_within_subject")
  expect_identical(fit_conf$meta$target, "subject_weighted_population_mean")
})

test_that("cluster bootstrap copies selected subjects intact", {
  id <- c(1, 1, 2, 2, 2, 3, 3, 3, 3)
  set.seed(21)
  W <- funbootband:::.bootstrap_weight_matrix(
    n = length(id), B = 80, iid = FALSE, id = id
  )

  expect_equal(rowSums(W), rep(1, nrow(W)), tolerance = 1e-12)
  for (g in unique(id)) {
    # All curves from a selected subject have the same positive weight; all
    # curves from an unselected subject have zero weight together.
    expect_true(all(apply(W[, id == g, drop = FALSE], 1, function(z) {
      length(unique(z)) == 1L
    })))
  }

  # Subject totals must be 0, 1/K, 2/K, ... according to its selection count.
  K <- length(unique(id))
  totals <- vapply(unique(id), function(g) rowSums(W[, id == g, drop = FALSE]),
                   numeric(nrow(W)))
  expect_equal(totals * K, round(totals * K), tolerance = 1e-12)
})

test_that("prediction kernel retains one maximum per pseudo-future curve", {
  Y <- rbind(c(0, 1, 2), c(0, 2, 4), c(0, 1, 3))
  W <- rbind(c(1 / 3, 1 / 3, 1 / 3), c(1 / 2, 1 / 2, 0))
  M <- funbootband:::prediction_curve_max_dev_weighted_cpp(Y, W, 1e-12)
  expect_equal(dim(M), c(2L, 3L))
  expect_true(all(is.finite(M)))
})

test_that("Rcpp kernels agree with transparent R reference calculations", {
  Y <- rbind(
    c(-1.0, 0.5, 2.0, 1.5),
    c( 0.0, 1.0, 3.0, 2.0),
    c( 1.0, 0.0, 2.5, 4.0)
  )
  W <- rbind(
    c(0.25, 0.25, 0.25, 0.25),
    c(0.50, 0.00, 0.25, 0.25),
    c(0.00, 0.50, 0.50, 0.00)
  )

  pred_cpp <- funbootband:::prediction_curve_max_dev_weighted_cpp(Y, W, 1e-12)
  pred_r <- .prediction_kernel_reference(Y, W)
  expect_equal(pred_cpp, pred_r, tolerance = 1e-12)

  mu_hat <- rowMeans(Y)
  conf_cpp <- funbootband:::confidence_max_dev_studentized_cpp(
    Y, mu_hat, W, 1e-12
  )
  conf_r <- .confidence_kernel_reference(Y, mu_hat, W)
  expect_equal(conf_cpp, conf_r, tolerance = 1e-12)
})

test_that("band() is reproducible when the random seed is reset", {
  set.seed(22)
  Y <- matrix(rnorm(20 * 12), nrow = 20)
  id <- rep(seq_len(4), each = 3)

  set.seed(987)
  fit1 <- band(Y, type = "prediction", iid = FALSE, id = id,
               B = 60, k.coef = 3)
  set.seed(987)
  fit2 <- band(Y, type = "prediction", iid = FALSE, id = id,
               B = 60, k.coef = 3)

  expect_equal(fit1, fit2, tolerance = 0)
})

test_that("constant curves produce finite zero-width bands", {
  Y <- matrix(3, nrow = 15, ncol = 8)

  set.seed(23)
  fitP <- band(Y, type = "prediction", iid = TRUE, B = 40, k.coef = 0)
  set.seed(23)
  fitC <- band(Y, type = "confidence", iid = TRUE, B = 40, k.coef = 0)

  expect_true(all(is.finite(c(fitP$lower, fitP$upper,
                              fitC$lower, fitC$upper))))
  expect_equal(fitP$lower, rep(3, nrow(Y)), tolerance = 1e-12)
  expect_equal(fitP$upper, rep(3, nrow(Y)), tolerance = 1e-12)
  expect_equal(fitC$lower, rep(3, nrow(Y)), tolerance = 1e-12)
  expect_equal(fitC$upper, rep(3, nrow(Y)), tolerance = 1e-12)
})

test_that("band() infers clusters from column-name prefixes when id is missing", {
  set.seed(3)
  T <- 25; n <- 12
  Y <- matrix(rnorm(T * n, sd = 0.3), nrow = T, ncol = n)
  # prefixes define clusters: subj1_* , subj2_* , subj3_*
  colnames(Y) <- c(paste0("subj1_rep", 1:4),
                   paste0("subj2_rep", 1:4),
                   paste0("subj3_rep", 1:4))
  expect_no_error(
    band(Y, type = "confidence", alpha = 0.10, iid = FALSE,
         B = 80, k.coef = 8L)
  )
})

test_that("excessive k.coef is clamped with a warning", {
  set.seed(10)
  Y <- matrix(rnorm(10 * 6), nrow = 10)

  expect_warning(
    fit <- band(
      Y,
      type = "prediction",
      iid = TRUE,
      B = 20L,
      k.coef = 50L
    ),
    "exceeds maximum 4"
  )

  expect_identical(fit$meta$k.coef, 4L)
})

test_that("Invalid inputs raise informative errors", {
  T <- 10; n <- 5
  Y <- matrix(rnorm(T * n), nrow = T)

  # invalid type
  expect_error(band(Y, type = "not-a-type", alpha = 0.1, iid = TRUE, B = 50))

  # alpha out of bounds
  expect_error(band(Y, type = "prediction", alpha = -0.1, iid = TRUE, B = 50),
               "alpha.*\\(0,1\\)", ignore.case = TRUE)
  expect_error(band(Y, type = "prediction", alpha = 1.1, iid = TRUE, B = 50),
               "alpha.*\\(0,1\\)", ignore.case = TRUE)

  # too few rows / cols
  expect_error(band(matrix(1, nrow = 1, ncol = 4), B = 10),
               "at least 2 time points", ignore.case = TRUE)
  expect_error(band(matrix(1, nrow = 4, ncol = 1), B = 10),
               "at least 2.*curves", ignore.case = TRUE)

  # iid = FALSE but no usable id or names
  Y2 <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colnames(Y2) <- NULL
  expect_error(
    band(Y2, type = "prediction", alpha = 0.1, iid = FALSE, B = 20),
    "supply `id` or.*column names", ignore.case = TRUE
  )

  expect_error(band(Y, B = 1), "B.*>= 2")
  expect_error(band(Y, B = 2.5), "B.*integer")
  expect_error(band(Y, k.coef = 1.5), "k.coef.*integer")
  expect_error(band(Y, alpha = c(0.05, 0.10)), "alpha.*one finite")

  id_na <- c(1, 1, 2, 2, NA)
  expect_error(
    band(Y, iid = FALSE, id = id_na, B = 20),
    "id.*missing"
  )
})

test_that("Non-numeric input and NA values are handled with errors", {
  T <- 12; n <- 6
  Y <- matrix(rnorm(T * n), nrow = T, ncol = n)

  # introduce NA -> underlying quantile/SD should complain
  Yna <- Y
  Yna[1, 1] <- NA
  expect_error(
    band(Yna, type = "prediction", alpha = 0.1, iid = TRUE, B = 40),
    "missing|NA", ignore.case = TRUE
  )

  # data.frame with non-numeric column -> coercion loses numeric type, should error
  DF <- as.data.frame(Y)
  DF$non_numeric <- rep("a", T)
  expect_error(
    band(DF, type = "prediction", alpha = 0.1, iid = TRUE, B = 40),
    "must be a numeric matrix", ignore.case = TRUE
  )
})
