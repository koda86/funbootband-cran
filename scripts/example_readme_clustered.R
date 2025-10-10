library(funbootband)

set.seed(2)
T <- 200
m <- c(5, 5)                          # curves per cluster
x <- seq(0, 1, length.out = T)

# Simulate smooth cluster-specific curves
mu <- list(
  function(z) 8 * sin(2 * pi * z),
  function(z) 8 * cos(2 * pi * z)
)

# Smooth within-cluster variation (two harmonic modes)
Bm <- cbind(sin(2 * pi * x), cos(2 * pi * x))

gen_curve <- function(k) {
  sc <- rnorm(ncol(Bm), sd = c(2.0, 1.5))
  mu[[k]](x) + as.vector(Bm %*% sc)
}

Ylist <- lapply(seq_along(m), function(k) {
  sapply(seq_len(m[k]), function(i) gen_curve(k) + rnorm(T, sd = 0.6))
})
Y <- do.call(cbind, Ylist)

# Name columns so band() can infer clusters from the prefix (iid = FALSE)
colnames(Y) <- unlist(mapply(
  function(k, mk) paste0("C", k, "_", seq_len(mk)),
  seq_along(m), m
))

# Fit simultaneous prediction and confidence bands (cluster-aware)
fit_pred <- band(Y, type = "prediction", alpha = 0.11, iid = FALSE, B = 200L, k.coef = 20L)
fit_conf <- band(Y, type = "confidence", alpha = 0.11, iid = FALSE, B = 200L, k.coef = 20L)

# Plotting
x_idx   <- seq_len(fit_pred$meta$T)
out_path <- "man/figures/README_clustered_plot.png"   # <- adjust if you prefer another path

png(filename = out_path, width = 2000, height = 1400, res = 300)

plot(x_idx, fit_pred$mean, type = "n",
     ylim = range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE),
     xlab = "Index (Time)", ylab = "Amplitude")

matlines(x_idx, Y, col = "grey70", lty = 1, lwd = 1)
polygon(c(x_idx, rev(x_idx)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x_idx, rev(x_idx)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = adjustcolor("gray40", alpha.f = 0.30), border = NA)
lines(x_idx, fit_pred$mean, lwd = 1, col = "black")

dev.off()
