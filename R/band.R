#' Simultaneous Bands for Functional Data
#'
#' Create simultaneous bootstrap bands for dense functional data
#' (rows are time points, columns are curves). For clustered designs, subjects
#' are treated as the independent sampling units and are resampled intact.
#'
#' @param data Numeric matrix with T rows (time) and n columns (curves).
#'   A data.frame of numeric columns is also accepted and coerced to a matrix.
#' @param type Character, either "prediction" or "confidence".
#' @param alpha Numeric in (0, 1). Use 0.05 for 95% bands.
#' @param iid Logical; if FALSE, use an intact-cluster bootstrap (requires `id`
#'   or infers clusters from column-name prefixes). The clustered prediction
#'   target is one future curve from a new subject. Subjects are weighted
#'   equally, and curves are weighted equally within subject.
#' @param id Optional integer/factor vector of length ncol(data) giving a cluster id
#'   for each curve (used when `iid = FALSE`). If NULL and `iid = FALSE`, clusters
#'   are inferred from column names by prefix (up to the first underscore, hyphen, or dot).
#' @param B Integer, number of bootstrap iterations (e.g., 1000 for final results;
#'   use smaller values in examples/tests).
#' @param k.coef Integer; number of Fourier harmonics (default 50).
#'   Automatically clamped to \eqn{\lfloor (T-2)/2 \rfloor} based on the grid
#'   length. This keeps sine/cosine harmonics in complete pairs on the
#'   periodic grid. Larger values fit more high-frequency detail; smaller
#'   values smooth more.
#'
#' @details
#' For `iid = FALSE`, the clustered prediction band is marginal over the
#' subject population. It targets one curve from an independent new subject;
#' it is not conditional on an already observed subject and does not target
#' joint coverage of several future curves. The construction assumes
#' independent subjects and exchangeable repeated curves within subject. Since
#' calibration follows Fourier preprocessing, the formal target is the future
#' curve in the same Fourier-reconstructed representation, rather than raw
#' pointwise measurement noise that the chosen basis does not retain.
#'
#' The i.i.d. calibration follows the curve-level functional-bootstrap target
#' of Lenhoff et al. (1999). The clustered construction implemented here is an
#' intact-subject bootstrap with subject-first empirical weighting. It is a
#' revision of, rather than a literal implementation of, the hierarchical
#' resampling description in Koska et al. (2023).
#' Clustered confidence inference treats subject mean curves as the independent
#' units and studentizes every bootstrap replicate with its own pointwise
#' standard error.
#'
#' @return An object of class `funbootband`, implemented as a list with elements
#'   `lower`, `mean`, `upper` (each of length T) and `meta`. Existing code can
#'   continue to access these components with `$`. For clustered prediction,
#'   `meta$target` records the estimand `"new_subject_new_curve"` and
#'   `meta$weighting` records the subject-first weighting convention.
#'
#' @example inst/examples/iid_example.R
#' @example inst/examples/clustered_example.R
#'
#' @references
#' Koska, D., Oriwol, D., & Maiwald, C. (2023).
#' Comparison of statistical models for characterizing continuous differences
#' between two biomechanical measurement systems.
#' *Journal of Biomechanics*, 149, 111506.
#' <doi:10.1016/j.jbiomech.2023.111506>
#'
#' Lenhoff, M. W., Santner, T. J., Otis, J. C., Peterson, M. G. E., Williams, B. J., & Backus, S. I. (1999).
#' Bootstrap prediction and confidence bands: a superior statistical method for analysis of gait data.
#' *Gait & Posture*, 9(1), 10–17.
#' <doi:10.1016/S0966-6362(98)00043-5>
#'
#' Davison, A. C., & Hinkley, D. V. (1997).
#' *Bootstrap Methods and Their Application*.
#' Cambridge University Press.
#' <doi:10.1017/cbo9780511802843>
#'
#' @export
band <- function(data,
                 type  = c("prediction","confidence"),
                 alpha = 0.05,
                 iid   = TRUE,
                 id    = NULL,
                 B     = 1000L,
                 k.coef = 50L) {

  type  <- match.arg(type)

  # ---- Input normalization & checks ----
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data)) stop("`data` must be a numeric matrix [T x n].")
  if (any(!is.finite(data))) stop("`data` must not contain NA/NaN/Inf.")
  Tlen <- nrow(data); ncur <- ncol(data)
  if (Tlen < 2L || ncur < 2L) stop("`data` must have at least 2 time points (rows) and 2 curves (cols).")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      !(alpha > 0 && alpha < 1)) {
    stop("`alpha` must be one finite number in (0,1).")
  }
  alpha <- as.numeric(alpha)
  if (!is.logical(iid) || length(iid) != 1L) stop("`iid` must be a logical scalar.")
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 2 || B != floor(B)) {
    stop("`B` must be one integer >= 2.")
  }
  B <- as.integer(B)

  # clusters
  if (!iid) {
    if (!is.null(id)) {
      if (length(id) != ncur) stop("`id` must have length ncol(data).")
      if (anyNA(id)) stop("`id` must not contain missing values.")
      id <- as.integer(as.factor(id))
    } else {
      cn <- colnames(data)
      if (is.null(cn) || any(!nzchar(cn)) || anyNA(cn)) {
        stop("For iid = FALSE, supply `id` or meaningful column names to infer clusters.")
      }
      id <- tolower(trimws(sub("(_|-|\\.).*$", "", cn)))
      id <- as.integer(as.factor(id))
    }
    cluster_sizes <- tabulate(id)
    if (length(cluster_sizes) < 2L) {
      stop("Clustered inference requires at least two subjects/clusters.")
    }
    if (!any(cluster_sizes > 1L)) {
      stop("Cluster structure has no repeated curves; use `iid = TRUE` or supply repeated curves per subject.")
    }
  } else {
    id <- NULL
    cluster_sizes <- NULL
  }

  # ---- Fourier preprocessing (Lenhoff-style) ----
  if (!is.numeric(k.coef) || length(k.coef) != 1L || !is.finite(k.coef) ||
      k.coef < 0 || k.coef != floor(k.coef)) {
    stop("`k.coef` must be one nonnegative integer.")
  }
  k.coef <- as.integer(k.coef)
  # The first and last grid points have the same Fourier phase, leaving T-1
  # distinct phases. With odd T, harmonic (T-1)/2 has an identically zero
  # sine column (the Nyquist frequency). Restrict to complete, independent
  # sine/cosine pairs; for even T this leaves the previous cap unchanged.
  maxK <- as.integer(floor((Tlen - 2L) / 2L))
  if (k.coef > maxK) {
    warning("`k.coef` = ", k.coef, " exceeds maximum ", maxK,
            " for T = ", Tlen, ". Using ", maxK, " instead.")
    k.coef <- maxK
  }
  fit <- fit_fourier(data, K = k.coef)
  data <- fit$fitted   # replace raw with Fourier-reconstructed curves

  # ---- Estimand and empirical distribution ----
  # IID: every curve has weight 1/n.
  # Clustered: first select a subject uniformly, then a curve uniformly within
  # that subject. Thus curve (i,j) has weight 1 / (K * m_i). This defines both
  # the subject-weighted mean and the new-subject/new-curve prediction target.
  future_weights <- if (iid) {
    rep.int(1 / ncur, ncur)
  } else {
    .cluster_curve_weights(id)
  }

  mu_hat <- as.numeric(data %*% future_weights)
  centered <- sweep(data, 1L, mu_hat, FUN = "-")
  sd_hat <- sqrt(as.numeric((centered^2) %*% future_weights))
  sd_hat <- .ridge_scale(sd_hat)

  if (type == "prediction") {
    # Bootstrap weights encode the training sample. In the clustered case, K
    # subjects are sampled with replacement and every selected subject is
    # copied intact. No second-stage curve resampling is performed.
    boot_weights <- .bootstrap_weight_matrix(
      n = ncur, B = B, iid = iid, id = id
    )
    # Each row is one bootstrap training sample; each column is one possible
    # future curve. Keeping curve-specific maxima (rather than maxing over all
    # curves) calibrates coverage for ONE future curve.
    M <- prediction_curve_max_dev_weighted_cpp(data, boot_weights, 1e-12)
    calibration_weights <- rep(future_weights, each = B)
    c_p <- .weighted_quantile(as.vector(M), calibration_weights, 1 - alpha)
    lower <- mu_hat - c_p * sd_hat
    upper <- mu_hat + c_p * sd_hat

  } else { # confidence
    if (iid) {
      independent_units <- data
    } else {
      # Subject mean curves, rather than individual repeated curves, are the
      # independent units for clustered confidence inference.
      independent_units <- .cluster_means(data, id)
    }
    n_units <- ncol(independent_units)
    se_hat <- apply(independent_units, 1L, stats::sd) / sqrt(n_units)
    se_hat <- .ridge_scale(se_hat)

    # Resample the independent units and studentize each replicate with its
    # own pointwise standard error. For clustered data, an independent unit is
    # a complete subject represented by its subject mean curve.
    unit_boot_weights <- .bootstrap_weight_matrix(
      n = n_units, B = B, iid = TRUE, id = NULL
    )
    C <- confidence_max_dev_studentized_cpp(
      independent_units, mu_hat, unit_boot_weights, 1e-12
    )
    c_c <- stats::quantile(C, probs = 1 - alpha, names = FALSE, type = 7)
    lower <- mu_hat - c_c * se_hat
    upper <- mu_hat + c_c * se_hat
  }

  out <- list(
    lower = as.numeric(lower),
    mean  = as.numeric(mu_hat),
    upper = as.numeric(upper),
    meta  = list(
      type = type,
      alpha = alpha,
      iid = iid,
      B = B,
      n = ncur,
      T = Tlen,
      k.coef = k.coef,
      n_clusters = if (iid) NA_integer_ else length(cluster_sizes),
      cluster_sizes = if (iid) NULL else as.integer(cluster_sizes),
      target = if (type == "prediction") {
        if (iid) "new_iid_curve" else "new_subject_new_curve"
      } else {
        if (iid) "iid_population_mean" else "subject_weighted_population_mean"
      },
      weighting = if (iid) {
        "equal_curve"
      } else {
        "equal_subject_then_equal_curve_within_subject"
      },
      bootstrap_unit = if (iid) "curve" else "intact_subject",
      curve_representation = "finite_fourier_reconstruction",
      engine = "cpp"
    )
  )

  # Keep the established list structure and component names while adding an
  # S3 class for print(), summary(), and plot() methods.
  class(out) <- c("funbootband", "list")
  out
}

# ----- internal helpers (do not export) -----

# Finite Fourier design, T x (2K+1)
fourier_design <- function(Tlen, K) {
  stopifnot(Tlen >= 2L, K >= 0L)
  t_idx <- 0:(Tlen - 1L)
  denom <- (Tlen - 1L)   # Lenhoff: '-1' for periodic closure
  X <- cbind(1, matrix(NA_real_, nrow = Tlen, ncol = 2L * K))
  col <- 2L
  if (K > 0L) {
    for (k in seq_len(K)) {
      X[, col] <- cos(2 * pi * k * t_idx / denom); col <- col + 1L
      X[, col] <- sin(2 * pi * k * t_idx / denom); col <- col + 1L
    }
  }
  X
}

# Fit Fourier series for all curves
fit_fourier <- function(data, K) {
  if (is.data.frame(data)) data <- as.matrix(data)
  stopifnot(is.matrix(data), is.numeric(data))
  Tlen <- nrow(data)
  X <- fourier_design(Tlen, K)
  qrX <- qr(X)                         # stable LS
  coef_mat <- qr.coef(qrX, data)       # (2K+1) x n
  if (is.null(dim(coef_mat))) coef_mat <- matrix(coef_mat, ncol = 1L)
  fitted <- X %*% coef_mat             # T x n
  list(X = X, coef = coef_mat, fitted = fitted)
}

# Equal probability for subjects, followed by equal probability for curves
# within subject. These weights sum to one even when cluster sizes differ.
.cluster_curve_weights <- function(id) {
  sizes <- tabulate(id)
  1 / (length(sizes) * sizes[id])
}

# T x K matrix of subject-specific mean curves.
.cluster_means <- function(data, id) {
  groups <- split(seq_len(ncol(data)), id)
  out <- vapply(groups, function(j) rowMeans(data[, j, drop = FALSE]),
                numeric(nrow(data)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  out
}

# B x n matrix of empirical bootstrap weights.
#
# IMPORTANT REVISION: clustered rows are generated by drawing K subjects with
# replacement and retaining every curve belonging to each selected subject.
# If subject i is selected c_i times, each of its curves receives weight
# c_i / (K * m_i). This is an intact-cluster bootstrap and handles unequal m_i.
.bootstrap_weight_matrix <- function(n, B, iid, id) {
  out <- matrix(0, nrow = B, ncol = n)

  if (iid || is.null(id)) {
    draws <- matrix(sample.int(n, B * n, replace = TRUE), nrow = B)
    for (b in seq_len(B)) out[b, ] <- tabulate(draws[b, ], nbins = n) / n
    return(out)
  }

  groups <- split(seq_len(n), id)
  K <- length(groups)
  draws <- matrix(sample.int(K, B * K, replace = TRUE), nrow = B)
  for (b in seq_len(B)) {
    counts <- tabulate(draws[b, ], nbins = K)
    for (g in which(counts > 0L)) {
      out[b, groups[[g]]] <- counts[g] / (K * length(groups[[g]]))
    }
  }
  out
}

.ridge_scale <- function(x, ridge = 1e-12) {
  x[!is.finite(x) | x < ridge] <- ridge
  x
}

# Left-continuous inverse of a weighted empirical CDF. This matches the
# coverage equation P(M <= c) >= p and avoids interpolation between observed
# curve-level maximum deviations.
.weighted_quantile <- function(x, w, p) {
  if (length(x) != length(w)) stop("Internal error: `x` and `w` lengths differ.")
  if (any(!is.finite(x)) || any(!is.finite(w)) || any(w < 0)) {
    stop("Internal error: invalid weighted-quantile inputs.")
  }
  keep <- w > 0
  x <- x[keep]
  w <- w[keep]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  x[which(cumsum(w) >= p)[1L]]
}
