#' Simultaneous Bands for Functional Data (Prediction or Confidence)
#'
#' @description
#' Create **simultaneous** bootstrap bands for dense functional data (rows = time, cols = curves).
#' Supports **hierarchical/clustered** designs via a simple cluster bootstrap.
#'
#' @param data Numeric matrix \code{[T x n]} with rows = time points and columns = curves.
#'             A data.frame of numeric columns is also accepted and will be coerced to a matrix.
#' @param type One of \code{c("prediction","confidence")}.
#'   - \code{"prediction"}: band for a **new individual** curve.
#'   - \code{"confidence"}: band for the **mean function** \eqn{\mu(t)}.
#' @param alpha Numeric in (0,1). Use \code{0.05} for 95\% bands.
#' @param iid Logical; if \code{FALSE}, use a **cluster bootstrap** (needs \code{id} or
#'            infers clusters from column-name prefixes).
#' @param id Optional integer/factor vector of length \code{ncol(data)} giving a cluster id
#'           for each curve (used when \code{iid = FALSE}). If \code{NULL}, and \code{iid = FALSE},
#'           clusters are inferred from column names by prefix (up to first `_`/`-`/`.`).
#' @param B Integer, number of bootstrap iterations (e.g. 1000 for final results; use smaller in examples/tests).
#' @return A list with
#'   \item{lower}{numeric vector length \code{T} (lower band).}
#'   \item{mean}{numeric vector length \code{T} (estimated mean).}
#'   \item{upper}{numeric vector length \code{T} (upper band).}
#'   \item{meta}{list with settings (type, alpha, iid, B, n, T).}
#'
#' @details
#' **Simultaneity** is enforced via the max-deviation approach:
#' - For \code{type="prediction"}: compute \eqn{M_i = \max_t |Y_i^\*(t) - \hat\mu(t)| / \hat\sigma(t)} and take
#'   \eqn{c_p = Q_{1-\alpha}(M_i)}; then band is \eqn{\hat\mu(t) \pm c_p\,\hat\sigma(t)}.
#' - For \code{type="confidence"}: bootstrap the mean curve \eqn{\hat\mu_i^\*(t)}, form
#'   \eqn{C_i = \max_t |\hat\mu_i^\*(t) - \hat\mu(t)| / \widehat{\mathrm{se}}(t)} with
#'   \eqn{\widehat{\mathrm{se}}(t) = \hat\sigma(t)/\sqrt{n}}, then \eqn{c_c = Q_{1-\alpha}(C_i)} and band is
#'   \eqn{\hat\mu(t) \pm c_c\,\widehat{\mathrm{se}}(t)}.
#'
#' Cluster bootstrap (\code{iid = FALSE}): for the **confidence** band we resample **clusters** with replacement
#' and include all member curves of each drawn cluster until we reach \code{n} curves; for the **prediction** band
#' we draw one cluster (uniformly over clusters) and then one of its member curves for the residual.
#'
#' @section Input checks:
#' The function errors if \code{data} is not numeric, has <2 rows or <2 columns, or if \code{iid = FALSE}
#' but no cluster structure can be determined.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' T <- 60; n <- 30
#' Y <- matrix(rnorm(T*n, sd = .3), nrow = T, ncol = n)
#'
#' # Prediction band (i.i.d.)
#' fitP <- band(Y, type = "prediction", alpha = 0.1, iid = TRUE, B = 200)
#'
#' # Confidence band (i.i.d.)
#' fitC <- band(Y, type = "confidence", alpha = 0.1, iid = TRUE, B = 200)
#' }
#'
#' @importFrom stats quantile sd
#' @export
band <- function(data,
                 type  = c("prediction","confidence"),
                 alpha = 0.05,
                 iid   = TRUE,
                 id    = NULL,
                 B     = 1000L) {

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

  # ---- Basic estimates ----
  mu_hat <- rowMeans(data)                      # mean curve, length T
  sd_hat <- apply(data, 1L, stats::sd)          # sample sd per time
  sd_hat[sd_hat == 0] <- 1e-12                  # ridge for stability

  if (type == "prediction") {
    # residuals matrix (T x n)
    resid_mat <- data - matrix(mu_hat, nrow = Tlen, ncol = ncur, byrow = FALSE)
    # pre-sample indices in R for determinism; C++ consumes them
    pick_idx  <- .resample_prediction_idx(n = ncur, B = B, iid = iid, id = id)
    M <- prediction_max_dev_cpp(resid_mat, mu_hat, sd_hat, pick_idx)
    c_p <- stats::quantile(M, probs = 1 - alpha, names = FALSE, type = 7)
    lower <- mu_hat - c_p * sd_hat
    upper <- mu_hat + c_p * sd_hat

  } else { # confidence
    se_hat <- sd_hat / sqrt(ncur)
    idx_mat <- .resample_confidence_idx_mat(n = ncur, B = B, iid = iid, id = id)
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

# pick one residual curve per bootstrap draw (length B)
.resample_prediction_idx <- function(n, B, iid, id) {
  if (iid || is.null(id)) {
    sample.int(n, B, replace = TRUE)
  } else {
    cl <- split(seq_len(n), id)
    cl_ids <- seq_along(cl)
    vapply(seq_len(B), function(i) {
      k <- sample(cl_ids, 1L)
      sample(cl[[k]], 1L)
    }, integer(1L))
  }
}

# return a B x n matrix of curve indices for bootstrap means
.resample_confidence_idx_mat <- function(n, B, iid, id) {
  if (iid || is.null(id)) {
    matrix(sample.int(n, B * n, replace = TRUE), nrow = B, ncol = n)
  } else {
    cl <- split(seq_len(n), id)
    cl_ids <- seq_along(cl)
    out <- matrix(NA_integer_, nrow = B, ncol = n)
    for (b in seq_len(B)) {
      idx <- integer(0L)
      while (length(idx) < n) {
        k <- sample(cl_ids, 1L)
        idx <- c(idx, cl[[k]])
      }
      out[b, ] <- idx[seq_len(n)]
    }
    out
  }
}

