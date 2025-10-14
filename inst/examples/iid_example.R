set.seed(1)
T <- 200
n <- 10
x <- seq(0, 1, length.out = T)

# Simulate smooth Gaussian-process-like curves of equal length
mu  <- 10 * sin(2 * pi * x)
ell <- 0.12; sig <- 3
Kmat <- outer(x, x, function(s, t) sig^2 * exp(-(s - t)^2 / (2 * ell^2)))
ev <- eigen(Kmat + 1e-8 * diag(T), symmetric = TRUE)
Z  <- matrix(rnorm(T * n), T, n)
Y  <- mu + ev$vectors %*% (sqrt(pmax(ev$values, 0)) * Z)
Y  <- Y + matrix(rnorm(T * n, sd = 0.2), T, n) # observation noise

# Fit prediction and confidence bands
fit_pred <- band(Y, type = "prediction", alpha = 0.11, iid = TRUE, B = 1000L, k.coef = 50L)
fit_conf <- band(Y, type = "confidence", alpha = 0.11, iid = TRUE, B = 1000L, k.coef = 50L)

# Plot the results
x_idx <- seq_len(fit_pred$meta$T)
ylim  <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)

plot(x_idx, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Index (Time)", ylab = "Amplitude",
     main = "Simultaneous bands (i.i.d.)")

matlines(x_idx, Y, col = "gray70", lty = 1, lwd = 1)
polygon(c(x_idx, rev(x_idx)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x_idx, rev(x_idx)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("gray40", alpha.f = 0.3), border = NA)
lines(x_idx, fit_pred$mean, col = "black", lwd = 1)

