#' Simultaneous Bands for Functional Functional Data (Prediction or Confidence)
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
  Tlen <- nrow(data); ncur <- ncol(data)
  if (Tlen < 2L || ncur < 2L) stop("`data` must have at least 2 time points (rows) and 2 curves (cols).")
  if (!(alpha > 0 && alpha < 1)) stop("`alpha` must be in (0,1).")
  if (!is.logical(iid) || length(iid) != 1L) stop("`iid` must be a logical scalar.")
  if (!iid) {
    # infer or validate clusters
    if (!is.null(id)) {
      if (length(id) != ncur) stop("`id` must have length ncol(data).")
      id <- as.integer(as.factor(id))
    } else {
      cn <- colnames(data)
      if (is.null(cn) || any(!nzchar(cn)) || anyNA(cn)) {
        stop("For iid = FALSE, supply `id` or meaningful column names to infer clusters.")
      }
      # cluster prefix up to first _ / - / .
      id <- tolower(trimws(sub("(_|-|\\.).*$", "", cn)))
      id <- as.integer(as.factor(id))
    }
    # sanity: at least one cluster with size > 1
    sizes <- tabulate(id)
    if (!any(sizes > 1L)) stop("Cluster structure inferred/provided has no cluster with size > 1.")
  }

  # ---- Basic estimates ----
  mu_hat <- rowMeans(data)                      # mean curve, length T
  # sample sd per time (avoid zero with small ridge)
  sd_hat <- apply(data, 1L, stats::sd)
  sd_hat[sd_hat == 0] <- 1e-12

  # helper: draw one residual curve index for prediction band
  draw_one_index <- function() {
    if (iid) {
      sample.int(ncur, 1L)
    } else {
      # cluster then member (uniform over clusters)
      cl_ids <- unique(id)
      cl <- sample(cl_ids, 1L)
      which(id == cl)[ sample.int(sum(id == cl), 1L) ]
    }
  }

  # helper: resample indices for mean bootstrap (confidence band)
  resample_indices_conf <- function() {
    if (iid) {
      sample.int(ncur, ncur, replace = TRUE)
    } else {
      # cluster bootstrap: sample clusters with replacement, include all their members
      out <- integer(0L)
      cl_ids <- unique(id)
      sizes  <- tabulate(id)
      # keep adding whole clusters until reaching >= ncur, then trim
      while (length(out) < ncur) {
        cl <- sample(cl_ids, 1L)
        out <- c(out, which(id == cl))
      }
      out[seq_len(ncur)]
    }
  }

  if (type == "prediction") {
    # M_i = max_t |Y*_i(t) - mu_hat(t)| / sd_hat(t), where Y*_i uses a resampled residual curve
    M <- numeric(B)
    # precompute residuals (T x n)
    resid_mat <- data - matrix(mu_hat, nrow = Tlen, ncol = ncur, byrow = FALSE)
    for (i in seq_len(B)) {
      pick <- draw_one_index()
      Y_star <- mu_hat + resid_mat[, pick]
      M[i] <- max(abs(Y_star - mu_hat) / sd_hat)
    }
    c_p <- stats::quantile(M, probs = 1 - alpha, names = FALSE, type = 7)
    lower <- mu_hat - c_p * sd_hat
    upper <- mu_hat + c_p * sd_hat

  } else { # confidence
    # bootstrap the mean curve; standard error se_hat = sd_hat / sqrt(n)
    se_hat <- sd_hat / sqrt(ncur)
    C <- numeric(B)
    for (i in seq_len(B)) {
      idx <- resample_indices_conf()
      mu_star <- rowMeans(data[, idx, drop = FALSE])
      C[i] <- max(abs(mu_star - mu_hat) / se_hat)
    }
    c_c <- stats::quantile(C, probs = 1 - alpha, names = FALSE, type = 7)
    lower <- mu_hat - c_c * se_hat
    upper <- mu_hat + c_c * se_hat
  }

  list(
    lower = as.numeric(lower),
    mean  = as.numeric(mu_hat),
    upper = as.numeric(upper),
    meta  = list(type = type, alpha = alpha, iid = iid, B = as.integer(B),
                 n = ncur, T = Tlen)
  )
}

