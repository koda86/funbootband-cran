#' Methods for `funbootband` objects
#'
#' Convenient display, summary, and plotting methods for objects returned by
#' [band()]. The underlying object remains a list, so the `lower`, `mean`,
#' `upper`, and `meta` components remain directly accessible.
#'
#' @param x A `funbootband` object.
#' @param object A `funbootband` object.
#' @param ... Additional arguments. For `plot()`, these are passed to
#'   [graphics::plot.default()]. They are currently ignored by `print()` and
#'   `summary()`.
#' @param grid Optional numeric vector giving the horizontal coordinates. By
#'   default, the method uses `x$meta$grid` when available and otherwise the
#'   integer sequence along the band.
#' @param ylim Optional numeric vector of length two giving the vertical plot
#'   limits. By default, the range of the lower and upper bands is used.
#' @param xlab,ylab Axis labels.
#' @param main Plot title. By default, the title reports the nominal coverage
#'   and band type.
#' @param band.col Fill color for the band.
#' @param border Border color for the band polygon. The default is no border.
#' @param mean.col Color for the estimated mean curve.
#' @param mean.lwd Line width for the estimated mean curve.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `summary.funbootband` containing the inferential target, design and
#'   computation settings, and summaries of band width. `plot()` returns `x`
#'   invisibly after drawing the band. The print method for a summary object
#'   returns that object invisibly.
#'
#' @name funbootband-methods
NULL

#' @rdname funbootband-methods
#' @export
print.funbootband <- function(x, ...) {
  .validate_funbootband(x)
  meta <- x$meta

  cat("<funbootband>\n")
  cat("  Band:      ", .band_label(meta), "\n", sep = "")
  cat("  Design:    ", .design_label(meta), "\n", sep = "")
  cat("  Target:    ", .target_label(meta$target), "\n", sep = "")
  cat("  Data:      ", meta$T, " grid points, ", meta$n, " curves\n",
      sep = "")
  cat("  Bootstrap: ", meta$B, " replicates\n", sep = "")

  invisible(x)
}

#' @rdname funbootband-methods
#' @export
summary.funbootband <- function(object, ...) {
  .validate_funbootband(object)
  meta <- object$meta
  widths <- object$upper - object$lower

  out <- list(
    type = meta$type,
    level = 1 - meta$alpha,
    design = if (isTRUE(meta$iid)) "independent_curves" else "clustered_curves",
    target = meta$target,
    weighting = meta$weighting,
    bootstrap_unit = meta$bootstrap_unit,
    n_curves = meta$n,
    n_grid_points = meta$T,
    n_clusters = meta$n_clusters,
    cluster_sizes = meta$cluster_sizes,
    bootstrap_replicates = meta$B,
    k.coef = meta$k.coef,
    curve_representation = meta$curve_representation,
    band_width = c(
      minimum = min(widths),
      median = stats::median(widths),
      mean = mean(widths),
      maximum = max(widths)
    )
  )
  class(out) <- c("summary.funbootband", "list")
  out
}

#' @rdname funbootband-methods
#' @export
print.summary.funbootband <- function(x, ...) {
  cat("Summary of <funbootband>\n")
  cat("  Band:      ", sprintf("%g%% simultaneous %s band",
                               100 * x$level, x$type), "\n", sep = "")
  cat("  Design:    ", if (identical(x$design, "independent_curves")) {
    "independent curves"
  } else {
    paste0("clustered curves (", x$n_clusters, " subjects)")
  }, "\n", sep = "")
  cat("  Target:    ", .target_label(x$target), "\n", sep = "")
  cat("  Weighting: ", .weighting_label(x$weighting), "\n", sep = "")
  cat("  Data:      ", x$n_grid_points, " grid points, ", x$n_curves,
      " curves\n", sep = "")
  if (!is.null(x$cluster_sizes)) {
    cat("  Clusters:  sizes ", paste(x$cluster_sizes, collapse = ", "),
        "\n", sep = "")
  }
  cat("  Fourier:   ", x$k.coef, " harmonics\n", sep = "")
  cat("  Bootstrap: ", x$bootstrap_replicates, " replicates\n", sep = "")
  cat("  Band width:\n")
  print(round(x$band_width, 6L))

  invisible(x)
}

#' @rdname funbootband-methods
#' @export
plot.funbootband <- function(x,
                             grid = NULL,
                             ylim = NULL,
                             xlab = NULL,
                             ylab = "Value",
                             main = NULL,
                             band.col = grDevices::adjustcolor(
                               "steelblue", alpha.f = 0.25
                             ),
                             border = NA,
                             mean.col = "black",
                             mean.lwd = 2,
                             ...) {
  .validate_funbootband(x)

  if (is.null(grid)) {
    grid <- x$meta$grid
    if (is.null(grid)) grid <- seq_along(x$mean)
  }
  if (!is.numeric(grid) || length(grid) != length(x$mean) ||
      any(!is.finite(grid))) {
    stop("`grid` must be a finite numeric vector with one value per grid point.")
  }
  if (is.unsorted(grid, strictly = TRUE)) {
    stop("`grid` must be strictly increasing.")
  }

  ylim_supplied <- !is.null(ylim)
  if (!ylim_supplied) {
    ylim <- range(c(x$lower, x$upper), finite = TRUE)
    if (ylim[1L] == ylim[2L]) {
      padding <- if (ylim[1L] == 0) 1 else 0.04 * abs(ylim[1L])
      ylim <- ylim + c(-padding, padding)
    }
  }
  if (!is.numeric(ylim) || length(ylim) != 2L || any(!is.finite(ylim)) ||
      ylim[1L] >= ylim[2L]) {
    stop("`ylim` must contain two finite increasing values.")
  }
  if (is.null(xlab)) {
    xlab <- if (is.null(x$meta$grid)) "Grid index" else "Grid"
  }
  if (is.null(main)) main <- .band_label(x$meta)

  graphics::plot.default(
    grid, x$mean,
    type = "n",
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  graphics::polygon(
    c(grid, rev(grid)),
    c(x$lower, rev(x$upper)),
    col = band.col,
    border = border
  )
  graphics::lines(grid, x$mean, col = mean.col, lwd = mean.lwd)

  invisible(x)
}

# Internal validation and labels shared by the public methods.
.validate_funbootband <- function(x) {
  required <- c("lower", "mean", "upper", "meta")
  if (!is.list(x) || !all(required %in% names(x)) || !is.list(x$meta)) {
    stop("`x` is not a valid `funbootband` object.")
  }

  n_grid <- length(x$mean)
  if (n_grid < 1L || length(x$lower) != n_grid || length(x$upper) != n_grid ||
      !all(vapply(x[c("lower", "mean", "upper")], is.numeric, logical(1L))) ||
      any(!is.finite(c(x$lower, x$mean, x$upper)))) {
    stop("`x` is not a valid `funbootband` object.")
  }

  required_meta <- c(
    "type", "alpha", "iid", "B", "n", "T", "k.coef", "n_clusters",
    "cluster_sizes", "target", "weighting", "bootstrap_unit",
    "curve_representation"
  )
  if (!all(required_meta %in% names(x$meta)) ||
      !identical(as.integer(x$meta$T), as.integer(n_grid))) {
    stop("`x` is not a valid `funbootband` object.")
  }

  invisible(x)
}

.band_label <- function(meta) {
  sprintf("%g%% simultaneous %s band", 100 * (1 - meta$alpha), meta$type)
}

.design_label <- function(meta) {
  if (isTRUE(meta$iid)) {
    "independent curves"
  } else {
    paste0("clustered curves (", meta$n_clusters, " subjects)")
  }
}

.target_label <- function(target) {
  switch(
    target,
    new_iid_curve = "one future independent curve",
    new_subject_new_curve = "one future curve from a new subject",
    iid_population_mean = "population mean function",
    subject_weighted_population_mean = "subject-weighted population mean function",
    target
  )
}

.weighting_label <- function(weighting) {
  switch(
    weighting,
    equal_curve = "equal weight per curve",
    equal_subject_then_equal_curve_within_subject =
      "equal subjects, then equal curves within subject",
    weighting
  )
}
