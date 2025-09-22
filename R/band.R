#' Simultaneous Bands for Functional Data
#'
#' Create simultaneous bootstrap bands for dense functional data
#' (rows are time points, columns are curves). Supports clustered designs
#' via a simple cluster bootstrap when `iid = FALSE`.
#'
#' @param data Numeric matrix with T rows (time) and n columns (curves).
#'   A data.frame of numeric columns is also accepted and coerced to a matrix.
#' @param type Character, either "prediction" or "confidence".
#' @param alpha Numeric in (0, 1). Use 0.05 for 95% bands.
#' @param iid Logical; if FALSE, use a cluster bootstrap (requires `id` or
#'   infers clusters from column-name prefixes).
#' @param id Optional integer/factor vector of length ncol(data) giving a cluster id
#'   for each curve (used when `iid = FALSE`). If NULL and `iid = FALSE`, clusters
#'   are inferred from column names by prefix (up to the first underscore, hyphen, or dot).
#' @param B Integer, number of bootstrap iterations (e.g., 1000 for final results;
#'   use smaller values in examples/tests).
#'
#' @return A list with elements `lower`, `mean`, `upper` (each of length T) and `meta`
#'   (a list with settings such as type, alpha, iid, B, n, T).
#'
#' @examples
#' ## i.i.d. example (small, fast) with shaded band
#' set.seed(1)
#' T <- 60; n <- 8
#' Y <- matrix(rnorm(T * n, sd = 0.25), nrow = T, ncol = n) +
#'      outer(seq_len(T), rep(1, n), function(i, j) 0.5 * sin(2*pi*i/T))
#' fit <- band(Y, type = "prediction", alpha = 0.1, iid = TRUE, B = 25L, k.coef = 12L)
#' x <- seq_len(fit$meta$T)
#' plot(x, fit$mean, type = "n",
#'      ylim = range(c(fit$lower, fit$upper)), xlab = "Index", ylab = "Value")
#' polygon(c(x, rev(x)), c(fit$lower, rev(fit$upper)),
#'         col = grDevices::adjustcolor("steelblue", alpha.f = 0.3), border = NA)
#' lines(x, fit$mean, lwd = 2)
#'
#' ## clustered example (also kept small and fast)
#' \donttest{
#' set.seed(2)
#' T  <- 80; m <- c(4, 4)                 # two clusters, few curves
#' t  <- seq(0, 1, length.out = T)
#' mu <- list(function(x) 0.7 * sin(2*pi*x),
#'            function(x) 0.6 * cos(2*pi*x))
#' Bm <- cbind(sin(2*pi*t), cos(2*pi*t))
#' gen_curve <- function(k) {
#'   sc <- rnorm(ncol(Bm), sd = c(0.2, 0.15))
#'   mu[[k]](t) + as.vector(Bm %*% sc) + rnorm(T, sd = 0.12)
#' }
#' Ylist <- lapply(seq_along(m), function(k) sapply(seq_len(m[k]), function(i) gen_curve(k)))
#' Yh    <- do.call(cbind, Ylist)
#' id    <- rep(seq_along(m), times = m)
#' fitH  <- band(Yh, type = "prediction", alpha = 0.1,
#'               iid = FALSE, id = id, B = 25L, k.coef = 12L)
#' xh <- seq_len(fitH$meta$T)
#' plot(xh, fitH$mean, type = "n",
#'      ylim = range(c(Yh, fitH$lower, fitH$upper), finite = TRUE),
#'      xlab = "Index", ylab = "Value")
#' polygon(c(xh, rev(xh)), c(fitH$lower, rev(fitH$upper)),
#'         col = grDevices::adjustcolor("steelblue", alpha.f = 0.30), border = NA)
#' lines(xh, fitH$mean,  lwd = 2)
#' }

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
  alpha <- as.numeric(alpha)

  # ---- Input normalization & checks ----
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data)) stop("`data` must be a numeric matrix [T x n].")
  if (any(!is.finite(data))) stop("`data` must not contain NA/NaN/Inf.")
  Tlen <- nrow(data); ncur <- ncol(data)
  if (Tlen < 2L || ncur < 2L) stop("`data` must have at least 2 time points (rows) and 2 curves (cols).")
  if (!(alpha > 0 && alpha < 1)) stop("`alpha` must be in (0,1).")
  if (!is.logical(iid) || length(iid) != 1L) stop("`iid` must be a logical scalar.")

  # clusters
  if (!iid) {
    if (!is.null(id)) {
      if (length(id) != ncur) stop("`id` must have length ncol(data).")
      id <- as.integer(as.factor(id))
    } else {
      cn <- colnames(data)
      if (is.null(cn) || any(!nzchar(cn)) || anyNA(cn)) {
        stop("For iid = FALSE, supply `id` or meaningful column names to infer clusters.")
      }
      id <- tolower(trimws(sub("(_|-|\\.).*$", "", cn)))
      id <- as.integer(as.factor(id))
    }
    if (!any(tabulate(id) > 1L)) stop("Cluster structure inferred/provided has no cluster with size > 1.")
  } else {
    id <- NULL
  }

  # ---- Fourier preprocessing (Lenhoff-style) ----
  k.coef <- as.integer(k.coef)
  if (k.coef < 0L) stop("`k.coef` must be nonnegative.")
  # practical ceiling for periodic Fourier basis
  maxK <- floor((Tlen - 1L) / 2L)
  if (k.coef > maxK) {
    warning("`k.coef` = ", k.coef, " exceeds maximum ", maxK,
            " for T = ", Tlen, ". Using ", maxK, " instead.")
    k.coef <- maxK
  }
  fit <- fit_fourier(data, K = k.coef)
  data <- fit$fitted   # replace raw with Fourier-reconstructed curves

  # ---- Basic estimates ----
  mu_hat <- rowMeans(data)                      # mean curve, length T
  sd_hat <- apply(data, 1L, stats::sd)          # sample sd per time
  sd_hat[sd_hat == 0] <- 1e-12                  # ridge for stability

  if (type == "prediction") {
    resid_mat <- data - matrix(mu_hat, nrow = Tlen, ncol = ncur, byrow = FALSE)
    idx_mat   <- .resample_idx_mat_two_stage(n = ncur, B = B, iid = iid, id = id)
    # pass the whole matrix (B x n)
    M <- prediction_max_dev_cpp(resid_mat, mu_hat, sd_hat, idx_mat)
    c_p <- stats::quantile(M, probs = 1 - alpha, names = FALSE, type = 7)
    lower <- mu_hat - c_p * sd_hat
    upper <- mu_hat + c_p * sd_hat

  } else { # confidence
    se_hat <- sd_hat / sqrt(ncur)
    idx_mat <- .resample_idx_mat_two_stage(n = ncur, B = B, iid = iid, id = id)
    C <- confidence_max_dev_cpp(data, mu_hat, se_hat, idx_mat)
    c_c <- stats::quantile(C, probs = 1 - alpha, names = FALSE, type = 7)
    lower <- mu_hat - c_c * se_hat
    upper <- mu_hat + c_c * se_hat
  }

  list(
    lower = as.numeric(lower),
    mean  = as.numeric(mu_hat),
    upper = as.numeric(upper),
    meta  = list(type = type, alpha = alpha, iid = iid, B = as.integer(B),
                 n = ncur, T = Tlen, engine = "cpp")
  )
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

# One curve per column draw, but clusters first (Stage 1), then curve within cluster (Stage 2)
# B x n matrix for bootstrap means (confidence) *and* prediction
.resample_idx_mat_two_stage <- function(n, B, iid, id) {
  if (iid || is.null(id)) {
    matrix(sample.int(n, B * n, replace = TRUE), nrow = B, ncol = n)
  } else {
    cl <- split(seq_len(n), id)
    cl_ids <- seq_along(cl)
    out <- matrix(NA_integer_, nrow = B, ncol = n)
    for (b in seq_len(B)) {
      # within a replicate, draw within clusters w/o replacement until a cluster is exhausted,
      # then continue with replacement
      avail <- lapply(cl, identity)
      for (j in seq_len(n)) {
        k <- sample(cl_ids, 1L)                      # Stage 1: pick a cluster
        if (length(avail[[k]]) == 0L) avail[[k]] <- cl[[k]]  # replenish -> with replacement
        pos <- sample.int(length(avail[[k]]), 1L)    # Stage 2: pick one curve within cluster
        out[b, j] <- avail[[k]][[pos]]
        avail[[k]] <- avail[[k]][-pos]
      }
    }
    out
  }
}

