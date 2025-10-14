## clustered (hierarchical) example

set.seed(2)
T <- 200
m <- c(5, 5)
x <- seq(0, 1, length.out = T)

# Cluster-specific means
mu <- list(
  function(z) 8 * sin(2 * pi * z),
  function(z) 8 * cos(2 * pi * z)
)

# Generate curves with smooth within-cluster variation
Bm <- cbind(sin(2 * pi * x), cos(2 * pi * x))
gen_curve <- function(k) {
  sc <- rnorm(ncol(Bm), sd = c(2.0, 1.5))
  mu[[k]](x) + as.vector(Bm %*% sc)
}

Ylist <- lapply(seq_along(m), function(k) {
  sapply(seq_len(m[k]), function(i) gen_curve(k) + rnorm(T, sd = 0.6))
})
Y <- do.call(cbind, Ylist)
colnames(Y) <- unlist(mapply(
  function(k, mk) paste0("C", k, "_", seq_len(mk)),
  seq_along(m), m
))


# Fit prediction and confidence bands
fit_pred <- band(Y, type = "prediction", alpha = 0.11, iid = FALSE, B = 1000L, k.coef = 50L)
fit_conf <- band(Y, type = "confidence", alpha = 0.11, iid = FALSE, B = 1000L, k.coef = 50L)

# Plot the results
x_idx <- seq_len(fit_pred$meta$T)
ylim   <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)

plot(x_idx, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Index (Time)", ylab = "Amplitude",
     main = "Simultaneous bands (clustered)")

matlines(x_idx, Y, col = "gray70", lty = 1, lwd = 1)
polygon(c(x_idx, rev(x_idx)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x_idx, rev(x_idx)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("gray40", alpha.f = 0.3), border = NA)
lines(x_idx, fit_pred$mean, col = "black", lwd = 1)
