library(funbootband)

set.seed(1)
T <- 200; n <- 10
x <- seq(0, 1, length.out = T)

# Simulate smooth Gaussian-process-like curves
mu  <- 10 * sin(2 * pi * x)
ell <- 0.12; sig <- 3
Kmat <- outer(x, x, function(s, t) sig^2 * exp(-(s - t)^2 / (2 * ell^2)))
ev <- eigen(Kmat + 1e-8 * diag(T), symmetric = TRUE)
Z  <- matrix(rnorm(T * n), T, n)
Y  <- mu + ev$vectors %*% (sqrt(pmax(ev$values, 0)) * Z)
Y  <- Y + matrix(rnorm(T * n, sd = 0.2), T, n)

# Fit simultaneous prediction and confidence bands
fit_pred <- band(Y, type = "prediction", alpha = 0.11, iid = TRUE, B = 200L, k.coef = 20L)
fit_conf <- band(Y, type = "confidence", alpha = 0.11, iid = TRUE, B = 200L, k.coef = 20L)

# Plotting
x <- seq_len(fit_pred$meta$T)

out_path <- "~/funbootband-cran/man/figures/README_iid_plot.png"

png(filename = out_path, width = 2000, height = 1400, res = 300)

# Plot grey curves + shaded ribbons
plot(x, fit_pred$mean, type = "n",
     ylim = range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE),
     xlab = "Index (Time)", ylab = "Amplitude")

matlines(x, Y, col = "grey70", lty = 1, lwd = 1)
polygon(c(x, rev(x)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x, rev(x)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = adjustcolor("gray40", alpha.f = 0.3), border = NA)
lines(x, fit_pred$mean, lwd = 1, col = "black")

dev.off()
