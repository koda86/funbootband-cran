#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Curve-specific simultaneous prediction errors under weighted bootstrap
// training distributions.
//
// IMPORTANT REVISION:
// - boot_weights[b, ] describes the whole bootstrap training distribution.
// - clustered weights are constructed in R by resampling subjects and copying
//   all their curves intact;
// - the return value is B x n, preserving one maximum per possible future
//   curve. The old kernel returned one maximum over all n curves per bootstrap
//   sample, which targeted simultaneous coverage of the complete observed set
//   rather than coverage of one future curve.
// [[Rcpp::export]]
Rcpp::NumericMatrix prediction_curve_max_dev_weighted_cpp(
    const Rcpp::NumericMatrix data,          // T x n reconstructed curves
    const Rcpp::NumericMatrix boot_weights,  // B x n, each row sums to one
    const double ridge = 1e-12) {

  const int T = data.nrow();
  const int n = data.ncol();
  const int B = boot_weights.nrow();

  if (boot_weights.ncol() != n) {
    Rcpp::stop("boot_weights must have ncol(data) columns.");
  }
  if (!R_finite(ridge) || ridge <= 0.0) {
    Rcpp::stop("ridge must be finite and positive.");
  }

  Rcpp::NumericMatrix out(B, n);
  std::vector<double> mu_star(T);
  std::vector<double> sd_star(T);

  for (int b = 0; b < B; ++b) {
    double weight_sum = 0.0;
    for (int j = 0; j < n; ++j) {
      const double w = boot_weights(b, j);
      if (!R_finite(w) || w < 0.0) {
        Rcpp::stop("boot_weights must be finite and nonnegative.");
      }
      weight_sum += w;
    }
    if (std::fabs(weight_sum - 1.0) > 1e-8) {
      Rcpp::stop("Each row of boot_weights must sum to one.");
    }

    // Bootstrap centre and pointwise plug-in scale.
    for (int t = 0; t < T; ++t) {
      double mu = 0.0;
      for (int j = 0; j < n; ++j) mu += boot_weights(b, j) * data(t, j);
      mu_star[t] = mu;

      double variance = 0.0;
      for (int j = 0; j < n; ++j) {
        const double difference = data(t, j) - mu;
        variance += boot_weights(b, j) * difference * difference;
      }
      sd_star[t] = std::max(std::sqrt(variance), ridge);
    }

    // One curve-level supremum for each empirical pseudo-future curve.
    for (int j = 0; j < n; ++j) {
      double maxdev = 0.0;
      for (int t = 0; t < T; ++t) {
        const double z = std::fabs((data(t, j) - mu_star[t]) / sd_star[t]);
        if (z > maxdev) maxdev = z;
      }
      out(b, j) = maxdev;
    }
  }

  return out;
}

// Studentized simultaneous confidence statistic.
//
// unit_data contains independent sampling units: individual curves for an
// i.i.d. analysis and subject mean curves for a clustered analysis. Each
// bootstrap replicate is studentized by its own pointwise standard error.
// [[Rcpp::export]]
Rcpp::NumericVector confidence_max_dev_studentized_cpp(
    const Rcpp::NumericMatrix unit_data,        // T x U independent units
    const Rcpp::NumericVector mu_hat,           // length T
    const Rcpp::NumericMatrix boot_weights,     // B x U, rows sum to one
    const double ridge = 1e-12)
{
  const int T = unit_data.nrow();
  const int U = unit_data.ncol();
  const int B = boot_weights.nrow();

  if (U < 2) Rcpp::stop("At least two independent units are required.");
  if (boot_weights.ncol() != U) {
    Rcpp::stop("boot_weights must have ncol(unit_data) columns.");
  }
  if (mu_hat.size() != T) {
    Rcpp::stop("mu_hat length must match nrow(unit_data).");
  }
  if (!R_finite(ridge) || ridge <= 0.0) {
    Rcpp::stop("ridge must be finite and positive.");
  }

  Rcpp::NumericVector out(B);
  std::vector<double> mu_star(T);
  std::vector<double> se_star(T);

  for (int b = 0; b < B; ++b) {
    double weight_sum = 0.0;
    for (int j = 0; j < U; ++j) {
      const double w = boot_weights(b, j);
      if (!R_finite(w) || w < 0.0) {
        Rcpp::stop("boot_weights must be finite and nonnegative.");
      }
      weight_sum += w;
    }
    if (std::fabs(weight_sum - 1.0) > 1e-8) {
      Rcpp::stop("Each row of boot_weights must sum to one.");
    }

    for (int t = 0; t < T; ++t) {
      double mean = 0.0;
      for (int j = 0; j < U; ++j) {
        mean += boot_weights(b, j) * unit_data(t, j);
      }
      mu_star[t] = mean;

      // With weights count_j / U, this equals sample_variance / U,
      // i.e. the squared standard error of the bootstrap mean.
      double weighted_variance = 0.0;
      for (int j = 0; j < U; ++j) {
        const double difference = unit_data(t, j) - mean;
        weighted_variance += boot_weights(b, j) * difference * difference;
      }
      se_star[t] = std::max(
        std::sqrt(weighted_variance / static_cast<double>(U - 1)), ridge
      );
    }

    double maxdev = 0.0;
    for (int t = 0; t < T; ++t) {
      const double z = std::fabs((mu_star[t] - mu_hat[t]) / se_star[t]);
      if (z > maxdev) maxdev = z;
    }
    out[b] = maxdev;
  }

  return out;
}
