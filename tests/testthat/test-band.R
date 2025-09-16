# Tests for funbootband::band()

test_that("band() runs for valid inputs (i.i.d., prediction & confidence)", {
  set.seed(1)
  T <- 40; n <- 20
  Y <- matrix(rnorm(T * n, sd = 0.3), nrow = T, ncol = n)

  # prediction band
  fitP <- band(Y, type = "prediction", alpha = 0.10, iid = TRUE, B = 120)
  expect_type(fitP, "list")
  expect_equal(sort(names(fitP)), sort(c("lower","mean","upper","meta")))
  expect_equal(length(fitP$lower), T)
  expect_equal(length(fitP$upper), T)
  expect_true(all(is.finite(fitP$lower)))
  expect_true(all(is.finite(fitP$upper)))
  expect_true(all(fitP$upper >= fitP$lower))

  # confidence band
  fitC <- band(Y, type = "confidence", alpha = 0.10, iid = TRUE, B = 120)
  expect_type(fitC, "list")
  expect_equal(length(fitC$lower), T)
  expect_true(all(fitC$upper >= fitC$lower))
})

test_that("band() works for clustered data via explicit id", {
  set.seed(2)
  T <- 30; n <- 24
  id <- rep(1:6, each = 4)                      # 6 clusters, each size 4
  Y  <- matrix(rnorm(T * n, sd = 0.25), nrow = T, ncol = n)

  fit <- band(Y, type = "prediction", alpha = 0.10, iid = FALSE, id = id, B = 100)
  expect_equal(length(fit$lower), T)
  expect_true(all(is.finite(fit$lower)))
  expect_true(all(is.finite(fit$upper)))
  expect_true(all(fit$upper >= fit$lower))
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
    band(Y, type = "confidence", alpha = 0.10, iid = FALSE, B = 80)
  )
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

