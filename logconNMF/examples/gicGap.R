# One sample with two Gaussian peaks and a constant background in every ion.
t <- seq_len(80)
H <- cbind(exp(-((t - 25) / 5)^2), exp(-((t - 53) / 4)^2))
S <- rbind(c(100, 30, 0), c(0, 20, 80))
Y <- round(2 + H %*% S)
X <- convMatrix(convKernel("gaussian", 15L, list(sigma = 2)), length(t))

# Estimate an ion-specific background level in both fits.
B <- matrix(1, nrow = length(t), ncol = 1L)
fit <- logconNMF(Y, X, time = t, background = B, truth.profiles = H,
                 ncores = 1L, correlation = "none")
gap <- gicGap(fit)
round(c(default = gap$fitted_gic, truth_initialised = gap$reference_gic), 3)
plot(gap, type = "GIC", main = "Two simulated peaks with a constant background")
