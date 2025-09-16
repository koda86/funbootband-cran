# Pre-sample indices in R (deterministic with set.seed), then let C++ "consume" them.

# For prediction: pick one residual curve per bootstrap draw
.resample_prediction_idx <- function(n, B, iid, id) {
  if (iid || is.null(id)) {
    sample.int(n, B, replace = TRUE)
  } else {
    # cluster-then-member
    cl <- split(seq_len(n), id)
    cl_ids <- seq_along(cl)
    vapply(seq_len(B), function(i) {
      k <- sample(cl_ids, 1L)
      sample(cl[[k]], 1L)
    }, integer(1L))
  }
}

# For confidence: B x n matrix of curve indices per draw (resample curves; clustered = clusters with replacement)
.resample_confidence_idx_mat <- function(n, B, iid, id) {
  if (iid || is.null(id)) {
    # i.i.d.: each row draws n curves with replacement
    matrix(sample.int(n, B * n, replace = TRUE), nrow = B, ncol = n)
  } else {
    # cluster bootstrap: sample clusters with replacement, include all members, trim to n
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

