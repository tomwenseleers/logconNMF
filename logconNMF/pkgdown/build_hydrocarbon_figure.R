# Refit the corrected gcmsData1 matrix and draw the six-component homepage
# figure. Run from the repository root:
#
# Rscript logconNMF/pkgdown/build_hydrocarbon_figure.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package to rebuild this figure.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
package <- if (file.exists("logconNMF/DESCRIPTION")) "logconNMF" else "."
output <- if (length(args)) args[[1L]] else
  file.path(package, "man", "figures", "gcms-data1-decomposition.png")

pkgload::load_all(package, quiet = TRUE)
data(gcmsData1, package = "logconNMF")

stopifnot(
  identical(dim(gcmsData1), c(119L, 670L)),
  identical(colnames(gcmsData1), as.character(30:699))
)

Y <- gcmsData1
time <- as.numeric(rownames(Y))
mass <- as.numeric(colnames(Y))
# A sensitivity grid over Gaussian widths and precluster counts was filtered to
# six-component solutions with chemically plausible apex positions.  Sigma = 9
# and seven preclusters gave the best balance of GIC, apex accuracy, compact
# profiles and recovery of the six compounds' diagnostic ions.
family <- quasipoisson("identity")
kernel <- convKernel("gaussian", 65L, list(sigma = 9))
# Use the same fixed-PSF construction as GCMSprime::gcmsDeco().  In particular,
# estimatePSF() applies the package's boundary-aware column projection to the
# raw convolution matrix.  Feeding convMatrix(kernel, nrow(Y)) directly would
# therefore describe a slightly different model at the chromatogram edges.
psf <- estimatePSF(
  Y,
  family = family,
  time = time,
  standard = NULL,
  nwindows = 1L,
  precluster = 1L,
  ncores = 8L,
  control = list(kernel_est = kernel, maxiter = 0L)
)
stopifnot(identical(psf$method, "fixed"))
X <- psf$X
B <- matrix(1, nrow(Y), 1L)
fit <- logconNMF(
  Y,
  X = X,
  background = B,
  family = family,
  precluster = 7L,
  correlation = "auto",
  ncores = 8L,
  control = list(initial_relaxed = FALSE, fusion_passes = 0L)
)
stopifnot(fit$rank == 6L)

background <- as.matrix(fit$background %*% fit$alpha)
corrected <- pmax(Y - background, 0)
profiles <- sweep(as.matrix(fit$H), 2L, apply(fit$S, 1L, max), "*")
spectra <- sweep(as.matrix(fit$S), 1L, apply(fit$S, 1L, max), "/") * 100
apex <- fit$apex_scan
# Match GCMSprime::plot.GCMSdeco() exactly so component colours remain stable
# between the compact overview and the detailed profile/spectrum panels.
colours <- grDevices::hcl(
  h = (seq_len(fit$rank) * 137.508) %% 360,
  c = 75,
  l = 48
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
grDevices::png(
  output,
  width = 2100,
  height = 1400,
  res = 180,
  bg = "white",
  type = if (.Platform$OS.type == "windows") "windows" else "cairo"
)
on.exit(grDevices::dev.off(), add = TRUE)

graphics::layout(
  rbind(c(1L, 1L, 1L), c(2L, 3L, 4L), c(5L, 6L, 7L)),
  heights = c(1.85, 1, 1)
)

# Top panel: every background-corrected ion trace in grey, with fitted
# base-peak-channel contributions on top. The square-root transform keeps weak
# ions visible while axis labels remain on the original intensity scale.
top_max <- max(corrected, profiles)
raw_ticks <- pretty(c(0, top_max), n = 5L)
raw_ticks <- raw_ticks[raw_ticks >= 0 & raw_ticks <= top_max]
graphics::par(mar = c(5.2, 6.0, 3.4, 1.2), mgp = c(3.6, .9, 0), las = 1)
graphics::matplot(
  time,
  sqrt(corrected),
  type = "l",
  lty = 1,
  lwd = .45,
  col = grDevices::adjustcolor("grey45", alpha.f = .16),
  xlab = "Retention time (min)",
  ylab = "Intensity",
  ylim = sqrt(c(0, top_max * 1.05)),
  yaxt = "n",
  main = paste("gcmsData1 | square-root display | six-component",
    "log-concave NMF | base-peak-scaled profiles")
)
graphics::axis(2, at = sqrt(raw_ticks), labels = format(raw_ticks, scientific = TRUE), las = 1)
graphics::abline(h = 0, col = "grey70", lty = 3)
for (k in seq_len(fit$rank)) {
  graphics::lines(time, sqrt(profiles[, k]), col = colours[k], lwd = 2.3)
  graphics::text(
    time[apex[k]],
    sqrt(profiles[apex[k], k]),
    labels = k,
    col = colours[k],
    font = 2,
    pos = 3,
    cex = .85,
    offset = .35
  )
}

# Bottom panels: base-peak-normalised resolved spectra, retaining the corrected
# nominal mass axis. Suppress numerical fitting dust below 0.1% of the base peak.
for (k in seq_len(fit$rank)) {
  s <- spectra[k, ]
  keep <- which(s >= .1)
  graphics::par(mar = c(3.8, 3.8, 3.1, .8), mgp = c(2.25, .72, 0), las = 1)
  graphics::plot(
    range(mass),
    c(0, 100),
    type = "n",
    xlab = "m/z",
    ylab = "Relative intensity (%)",
    main = sprintf("Component %d | apex %.3f min", k, time[apex[k]]),
    cex.main = .83,
    cex.axis = .78,
    cex.lab = .82,
    xaxs = "i"
  )
  graphics::segments(mass[keep], 0, mass[keep], s[keep], col = colours[k], lwd = 1.15)
}

grDevices::dev.off()
on.exit(NULL, add = FALSE)

cat("Wrote", normalizePath(output, winslash = "/"), "\n")
cat("Rank:", fit$rank, "\n")
cat("Quasi-GIC:", format(GIC(fit), digits = 10L), "\n")
cat("Total fitting seconds:", fit$seconds[["total"]], "\n")
cat("Apices (min):", paste(format(time[apex], digits = 7L), collapse = ", "), "\n")
