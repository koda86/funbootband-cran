test_that("print() reports the main inferential settings", {
  set.seed(101)
  Y <- matrix(rnorm(24 * 10), nrow = 24)
  fit <- band(Y, type = "prediction", alpha = 0.10,
              iid = TRUE, B = 30L, k.coef = 4L)

  printed <- capture.output(returned <- print(fit))
  printed_text <- paste(printed, collapse = "\n")

  expect_identical(returned, fit)
  expect_match(printed_text, "90% simultaneous prediction band")
  expect_match(printed_text, "independent curves")
  expect_match(printed_text, "one future independent curve")
  expect_match(printed_text, "24 grid points, 10 curves")
  expect_match(printed_text, "30 replicates")
})

test_that("summary() returns settings and band-width summaries", {
  set.seed(102)
  id <- rep(seq_len(4L), each = 3L)
  Y <- matrix(rnorm(20 * length(id)), nrow = 20)
  fit <- band(Y, type = "confidence", alpha = 0.05,
              iid = FALSE, id = id, B = 30L, k.coef = 3L)

  out <- summary(fit)

  expect_s3_class(out, "summary.funbootband")
  expect_identical(out$type, "confidence")
  expect_equal(out$level, 0.95)
  expect_identical(out$design, "clustered_curves")
  expect_identical(out$target, "subject_weighted_population_mean")
  expect_equal(out$n_clusters, 4L)
  expect_equal(out$cluster_sizes, rep(3L, 4L))
  expect_named(out$band_width, c("minimum", "median", "mean", "maximum"))
  expect_true(all(is.finite(out$band_width)))
  expect_lte(out$band_width[["minimum"]], out$band_width[["maximum"]])

  printed <- capture.output(returned <- print(out))
  printed_text <- paste(printed, collapse = "\n")
  expect_identical(returned, out)
  expect_match(printed_text, "95% simultaneous confidence band")
  expect_match(printed_text, "4 subjects")
  expect_match(printed_text, "subject-weighted population mean function")
})

test_that("plot() draws a band and accepts an explicit display grid", {
  set.seed(103)
  Y <- matrix(rnorm(18 * 8), nrow = 18)
  fit <- band(Y, type = "prediction", alpha = 0.10,
              iid = TRUE, B = 30L, k.coef = 3L)
  grid <- seq(0, 1, length.out = nrow(Y))
  plot_file <- tempfile(fileext = ".pdf")

  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  expect_invisible(plot(fit, grid = grid, xlab = "Normalized time"))
  grDevices::dev.off()

  expect_true(file.exists(plot_file))
  expect_gt(file.info(plot_file)$size, 0)
})

test_that("plot() rejects invalid display grids", {
  set.seed(104)
  Y <- matrix(rnorm(16 * 8), nrow = 16)
  fit <- band(Y, B = 20L, k.coef = 3L)

  expect_error(plot(fit, grid = 1:3), "one value per grid point")
  expect_error(plot(fit, grid = rev(seq_len(nrow(Y)))), "strictly increasing")
})
