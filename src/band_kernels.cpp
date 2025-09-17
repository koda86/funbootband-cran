#include <Rcpp.h>
using namespace Rcpp;

// Compute max deviation for prediction band
// [[Rcpp::export]]
Rcpp::NumericVector prediction_max_dev_cpp(const Rcpp::NumericMatrix resid_mat, // T x n, = data - mu_hat
                                           const Rcpp::NumericVector mu_hat,    // length T
                                           const Rcpp::NumericVector sd_hat,    // length T (already ridged in R)
                                           const Rcpp::IntegerMatrix idx_mat)   // B x n, 1..n (R indices)
{
  const int T = resid_mat.nrow();
  const int n = resid_mat.ncol();
  const int B = idx_mat.nrow();

  if (mu_hat.size() != T) Rcpp::stop("mu_hat length must match nrow(resid_mat).");
  if (sd_hat.size() != T) Rcpp::stop("sd_hat length must match nrow(resid_mat).");

  Rcpp::NumericVector M(B);

  // For each bootstrap replicate:
  // 1) form mu_star(t) = mean over columns idx_mat[b, ]
  // 2) form delta_resid(t) = mu_star(t) - mu_hat(t)
  // 3) for ALL curves k, compute max_t | (resid_mat(t,k) - delta_resid(t)) / sd_hat[t] |
  for (int b = 0; b < B; ++b) {
    // step 1: mu_star
    std::vector<double> mu_star(T, 0.0);
    for (int j = 0; j < n; ++j) {
      int col = idx_mat(b, j) - 1;               // R (1-based) -> C++ (0-based)
      if (col < 0 || col >= n) Rcpp::stop("idx_mat out of bounds.");
      for (int t = 0; t < T; ++t) {
        mu_star[t] += resid_mat(t, col) + mu_hat[t]; // resid + mu_hat = original data
      }
    }
    for (int t = 0; t < T; ++t) mu_star[t] /= static_cast<double>(n);

    // step 2: delta_resid = mu_star - mu_hat
    std::vector<double> delta_resid(T);
    for (int t = 0; t < T; ++t) delta_resid[t] = mu_star[t] - mu_hat[t];

    // step 3: max over k,t
    double maxdev = 0.0;
    for (int k = 0; k < n; ++k) {
      for (int t = 0; t < T; ++t) {
        double z = std::fabs((resid_mat(t, k) - delta_resid[t]) / sd_hat[t]);
        if (z > maxdev) maxdev = z;
      }
    }
    M[b] = maxdev;
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

