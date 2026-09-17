## Independent-curve example

set.seed(1)
T <- 101L
n <- 30L
x <- seq(0, 1, length.out = T)
mu_true <- 0.7 * sin(2 * pi * x) - 0.2 * cos(4 * pi * x)

generate_curve <- function() {
  mu_true +
    rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.30) * sin(2 * pi * x) +
    rnorm(1, sd = 0.20) * cos(2 * pi * x) +
    rnorm(1, sd = 0.15) * sin(4 * pi * x)
}

Y <- replicate(n, generate_curve())

fit_pred <- band(Y, type = "prediction", alpha = 0.10,
                 iid = TRUE, B = 500L, k.coef = 4L)
fit_conf <- band(Y, type = "confidence", alpha = 0.10,
                 iid = TRUE, B = 500L, k.coef = 4L)

ylim <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)
plot(x, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Normalized time", ylab = "Value",
     main = "Simultaneous bands (i.i.d.)")
matlines(x, Y, col = grDevices::adjustcolor("gray40", 0.25), lty = 1)
polygon(c(x, rev(x)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
polygon(c(x, rev(x)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("darkorange", 0.30), border = NA)
lines(x, fit_pred$mean, lwd = 2)
lines(x, mu_true, col = "red", lwd = 2, lty = 2)
