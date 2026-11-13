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
X <- convMatrix(convKernel("gaussian", 65L, list(sigma = 8)), nrow(Y))
B <- matrix(1, nrow(Y), 1L)
fit <- logconNMF(
  Y,
  X = X,
  background = B,
  family = quasipoisson("identity"),
  precluster = 8L,
  correlation = "auto",
  ncores = 8L
)
stopifnot(fit$rank == 6L)

background <- as.matrix(fit$background %*% fit$alpha)
corrected <- pmax(Y - background, 0)
profiles <- sweep(as.matrix(fit$H), 2L, apply(fit$S, 1L, max), "*")
spectra <- sweep(as.matrix(fit$S), 1L, apply(fit$S, 1L, max), "/") * 100
apex <- fit$apex_scan
colours <- grDevices::hcl.colors(fit$rank, "Dark 3")

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
grDevices::png(
  output,
  width = 2100,
  height = 1300,
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
graphics::par(mar = c(4.2, 6.0, 3.4, 1.2), mgp = c(3.9, .9, 0), las = 1)
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
  main = "gcmsData1 | square-root display | six-component log-concave NMF"
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
  graphics::par(mar = c(3.6, 3.8, 2.5, .8), mgp = c(2.2, .72, 0), las = 1)
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
