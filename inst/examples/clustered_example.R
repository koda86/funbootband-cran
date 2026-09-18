## Clustered example: repeated curves nested within subjects

set.seed(2)
T <- 101L
x <- seq(0, 1, length.out = T)

# Twelve independent subjects contribute unequal numbers of repeated curves.
K_subject <- 12L
m <- rep(c(2L, 3L, 4L), length.out = K_subject)
id <- rep(seq_len(K_subject), m)

mu_true <- 0.7 * sin(2 * pi * x) - 0.2 * cos(4 * pi * x)

# Smooth subject-specific deviations from the population mean.
subject_effect <- sapply(seq_len(K_subject), function(i) {
  rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.30) * sin(2 * pi * x) +
    rnorm(1, sd = 0.20) * cos(2 * pi * x)
})

# Smooth curve-to-curve deviations within a subject.
within_subject_effect <- function() {
  rnorm(1, sd = 0.18) * sin(4 * pi * x) +
    rnorm(1, sd = 0.12) * cos(4 * pi * x)
}

Y <- sapply(seq_along(id), function(j) {
  mu_true + subject_effect[, id[j]] + within_subject_effect()
})

trial <- ave(id, id, FUN = seq_along)
colnames(Y) <- paste0("subject", id, "_trial", trial)

# The prediction target is one Fourier-reconstructed curve from a new subject.
fit_pred <- band(Y, type = "prediction", alpha = 0.10,
                 iid = FALSE, id = id, B = 500L, k.coef = 4L)

# The confidence target is the equally subject-weighted population mean curve.
fit_conf <- band(Y, type = "confidence", alpha = 0.10,
                 iid = FALSE, id = id, B = 500L, k.coef = 4L)

fit_pred$meta[c("target", "weighting", "bootstrap_unit", "n_clusters")]

# Plot the results.
ylim <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)
plot(x, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Normalized time", ylab = "Value",
     main = "Simultaneous bands (clustered)")
matlines(x, Y, col = grDevices::adjustcolor("gray40", 0.20), lty = 1)
polygon(c(x, rev(x)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
polygon(c(x, rev(x)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("darkorange", 0.30), border = NA)
lines(x, fit_pred$mean, lwd = 2)
lines(x, mu_true, col = "red", lwd = 2, lty = 2)
