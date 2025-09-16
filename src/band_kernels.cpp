#include <Rcpp.h>
using namespace Rcpp;

// Compute max deviation for prediction band
// [[Rcpp::export]]
NumericVector prediction_max_dev_cpp(const NumericMatrix resid_mat, // T x n
                                     const NumericVector mu_hat,    // length T
                                     const NumericVector sd_hat,    // length T (ridge applied in R)
                                     const IntegerVector pick_idx)  // length B, 1..n (R indices)
{
  const int T = resid_mat.nrow();
  const int n = resid_mat.ncol();
  const int B = pick_idx.size();

  if (mu_hat.size() != T) stop("mu_hat length must match nrow(resid_mat).");
  if (sd_hat.size() != T) stop("sd_hat length must match nrow(resid_mat).");

  NumericVector M(B);

  for (int b = 0; b < B; ++b) {
    int j = pick_idx[b] - 1; // 1-based -> 0-based
    if (j < 0 || j >= n) stop("pick_idx out of bounds.");

    double m = 0.0;
    for (int t = 0; t < T; ++t) {
      double ystar = mu_hat[t] + resid_mat(t, j);
      double z = std::fabs((ystar - mu_hat[t]) / sd_hat[t]);
      if (z > m) m = z;
    }
    M[b] = m;
  }
  return M;
}

// Compute max deviation for confidence band
// [[Rcpp::export]]
NumericVector confidence_max_dev_cpp(const NumericMatrix data,    // T x n
                                     const NumericVector mu_hat,  // length T
                                     const NumericVector se_hat,  // length T
                                     const IntegerMatrix idx_mat) // B x n, 1..n
{
  const int T = data.nrow();
  const int n = data.ncol();
  const int B = idx_mat.nrow();

  if (mu_hat.size() != T) stop("mu_hat length must match nrow(data).");
  if (se_hat.size() != T) stop("se_hat length must match nrow(data).");

  NumericVector C(B);

  // For each bootstrap replicate: mu_star = rowMeans(data[, idx, drop=FALSE])
  for (int b = 0; b < B; ++b) {
    double maxdev = 0.0;

    for (int t = 0; t < T; ++t) {
      double rowsum = 0.0;
      for (int k = 0; k < n; ++k) {
        int j = idx_mat(b, k) - 1; // 1-based -> 0-based
        if (j < 0 || j >= n) stop("idx_mat out of bounds.");
        rowsum += data(t, j);
      }
      double mu_star_t = rowsum / static_cast<double>(n);
      double z = std::fabs((mu_star_t - mu_hat[t]) / se_hat[t]);
      if (z > maxdev) maxdev = z;
    }

    C[b] = maxdev;
  }

  return C;
}

