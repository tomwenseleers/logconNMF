args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript build_benchmark_figures.R <benchmark-results> <output-dir>")
}

results <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], mustWork = TRUE)

for (package in c("GCMSbenchmark", "ggplot2", "scales")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required.")
  }
}

read_table <- function(file) {
  utils::read.csv(file.path(results, file), check.names = FALSE)
}

navy <- "#14273b"
steelblue <- "#4682b4"
light_grey <- "#d9d9d9"
muted <- "#526579"
grid_colour <- "#dfe8eb"

matches <- read_table("matches.csv")
matches <- matches[matches$method == "GCMSprime", , drop = FALSE]
truth <- read_table("truth-peaks.csv")
summary <- read_table("benchmark.csv")
summary <- summary[summary$method == "GCMSprime", , drop = FALSE]
hits <- read_table("spectral-hits.csv")
hits <- hits[hits$method == "GCMSprime" & hits$present, , drop = FALSE]

area_regression <- utils::getFromNamespace("area_regression", "GCMSbenchmark")
area <- area_regression(matches)
grid <- seq(min(area$data$log_true), max(area$data$log_true), length.out = 200L)
interval <- stats::predict(
  area$fit,
  newdata = data.frame(log_true = grid),
  interval = "prediction"
)
limits <- range(c(area$data$log_true, area$data$log_estimated, interval))
padding <- 0.04 * diff(limits)
limits <- limits + c(-padding, padding)

grDevices::png(
  file.path(output_dir, "simultaneous-benchmark-area-accuracy.png"),
  width = 1600L, height = 1080L, res = 160L
)
graphics::par(mar = c(5.2, 5.5, 5.0, 1.4), mgp = c(3.25, 0.85, 0),
  las = 1, cex.axis = 1.0)
graphics::plot(
  area$data$log_true, area$data$log_estimated,
  type = "n", xlim = limits, ylim = limits,
  xlab = expression(log[10] * "(true peak area, integrated counts)"),
  ylab = expression(log[10] * "(estimated peak area, integrated counts)")
)
graphics::grid(col = grid_colour, lty = 1)
graphics::polygon(
  c(grid, rev(grid)), c(interval[, "lwr"], rev(interval[, "upr"])),
  col = light_grey, border = NA
)
graphics::lines(
  grid, interval[, "fit"],
  col = "black", lwd = 3
)
graphics::abline(
  0, 1, col = "black", lty = 2, lwd = 2
)
graphics::points(
  area$data$log_true, area$data$log_estimated,
  pch = 16, cex = 0.95,
  col = grDevices::adjustcolor(steelblue, alpha.f = 0.5)
)
graphics::title("Peak-area accuracy of the simultaneous logconNMF fit",
  col.main = navy, font.main = 2, cex.main = 1.45)
graphics::mtext(
  sprintf(
    "Robust adjusted R-squared %.3f | median absolute area error %.2f%% | %d matched occurrences",
    area$adjusted_R2, 100 * summary$area_median_error_factor, nrow(area$data)
  ),
  side = 3, line = 1.0, col = muted, cex = 0.96
)
graphics::legend(
  "topleft", bty = "n", inset = 0.015, cex = 0.95,
  legend = c("Matched occurrence", "Robust log-log fit",
    "95% robust prediction band", "1:1 reference"),
  col = c(grDevices::adjustcolor(steelblue, alpha.f = 0.5),
    "black", light_grey, "black"),
  pch = c(16, NA, 15, NA), lty = c(NA, 1, NA, 2),
  lwd = c(NA, 3, NA, 2), pt.cex = c(1.35, 1, 1.5, 1)
)
grDevices::dev.off()

distribution_plot <- function(scores, ticks, file, title, axis, statistic,
    statistic_name, caption, cap = NULL, reference = NULL) {
  scores <- scores[is.finite(scores)]
  display <- function(x) {
    if (!is.null(cap)) x <- pmin(x, cap)
    x
  }
  span <- range(stats::qlogis(display(c(scores, ticks))), finite = TRUE)
  transform <- function(x) {
    x <- display(x)
    answer <- stats::qlogis(x)
    answer[x == 0] <- span[[1L]] - 1.4
    answer[x == 1] <- span[[2L]] + 1.4
    answer
  }
  data <- data.frame(
    model = "Simultaneous logconNMF",
    score = scores,
    position = transform(scores)
  )
  quartiles <- stats::quantile(scores, c(0.25, 0.75), names = FALSE)
  centre <- statistic(scores)
  summary_data <- data.frame(
    model = "Simultaneous logconNMF",
    lower = transform(quartiles[[1L]]),
    upper = transform(quartiles[[2L]]),
    centre = transform(centre)
  )
  tick_labels <- vapply(ticks, function(x) {
    if (!is.null(cap) && x >= cap) "1" else
      format(x, trim = TRUE, scientific = FALSE)
  }, character(1L))

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = model, y = position)) +
    ggplot2::geom_violin(fill = light_grey, colour = "grey45", scale = "width",
      width = 0.72, adjust = 0.75, linewidth = 0.6) +
    ggplot2::geom_point(colour = steelblue, alpha = 0.55, size = 1.15,
      position = ggplot2::position_jitter(width = 0.20, height = 0,
        seed = 23L)) +
    ggplot2::geom_linerange(data = summary_data,
      ggplot2::aes(x = model, ymin = lower, ymax = upper),
      inherit.aes = FALSE, colour = "black", linewidth = 1.15) +
    ggplot2::geom_point(data = summary_data,
      ggplot2::aes(x = model, y = centre), inherit.aes = FALSE,
      shape = 21, fill = steelblue, colour = "black", stroke = 1.0,
      size = 3.2) +
    ggplot2::scale_y_continuous(
      breaks = transform(ticks), labels = tick_labels,
      limits = range(transform(c(scores, ticks))) + c(-0.25, 0.25),
      oob = scales::oob_keep, expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::labs(
      x = NULL, y = axis, title = title,
      subtitle = sprintf("%s %.6f | %d correctly matched true spectra | bar: interquartile range",
        statistic_name, centre, length(scores)),
      caption = caption
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", colour = navy, size = 20),
      plot.subtitle = ggplot2::element_text(colour = muted,
        margin = ggplot2::margin(4, 0, 9, 0)),
      plot.caption = ggplot2::element_text(hjust = 0, colour = muted),
      axis.text.y = ggplot2::element_text(face = "bold", colour = navy),
      plot.margin = ggplot2::margin(16, 18, 12, 16)
    )
  if (!is.null(reference)) {
    plot <- plot + ggplot2::geom_hline(yintercept = transform(reference),
      linetype = "dotted", colour = muted, linewidth = 0.55)
  }
  ggplot2::ggsave(file.path(output_dir, file), plot,
    width = 10.0, height = 5.5, dpi = 160, bg = "white")
}

distribution_plot(
  hits$spectral_cosine,
  ticks = c(0.5, 0.7, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999),
  file = "simultaneous-benchmark-spectral-accuracy.png",
  title = "Accuracy of recovered spectra",
  axis = "Spectral cosine similarity (logit scale)",
  statistic = stats::median,
  statistic_name = "Median",
  caption = "Reciprocal matches to compounds present in the sample; values at or above 0.999999 are displayed at 1.",
  cap = 0.999999, reference = 0.9
)

distribution_plot(
  hits$spectral_F1,
  ticks = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1),
  file = "simultaneous-benchmark-spectral-sparsity.png",
  title = "Accuracy of spectral sparsity",
  axis = "Within-spectrum ion-support F1 (logit scale)",
  statistic = mean,
  statistic_name = "Mean",
  caption = "Reference and recovered ion support use the same cutoff: 0.1% of the spectrum's base-peak intensity."
)

profiles <- readRDS(file.path(results, "profiles.rds"))
profiles <- profiles[profiles$method %in% c("Ground truth", "GCMSprime"), , drop = FALSE]
quality <- read_table("profile-quality.csv")
quality <- quality[quality$method == "GCMSprime", , drop = FALSE]
profile_colours <- grDevices::colorRampPalette(
  c("#144bb8", "#189ad5", "#f5c343", "#d82b28")
)(16L)

draw_curves <- function(data, line_type, line_width) {
  for (rank in seq_len(16L)) {
    pieces <- split(
      data[data$sample_rank == rank, , drop = FALSE],
      data$component[data$sample_rank == rank]
    )
    for (curve in pieces) {
      curve <- curve[order(curve$ri), , drop = FALSE]
      graphics::lines(curve$ri, curve$intensity,
        col = profile_colours[[rank]], lty = line_type, lwd = line_width)
    }
  }
}

grDevices::png(
  file.path(output_dir, "simultaneous-benchmark-elution-profiles.png"),
  width = 1900L, height = 950L, res = 160L
)
graphics::par(mfrow = c(1L, 2L), mar = c(4.6, 5.2, 5.2, 0.8),
  oma = c(2.3, 0.2, 1.5, 0.2), mgp = c(3.0, 0.8, 0), las = 1)
for (kind in c("Regular", "Acid")) {
  data <- profiles[profiles$kind == kind, , drop = FALSE]
  known <- data[data$method == "Ground truth", , drop = FALSE]
  window <- range(known$ri[known$intensity > max(known$intensity) * 1e-4])
  window <- window + c(-0.10, 0.10) * diff(window)
  data <- data[data$ri >= window[[1L]] & data$ri <= window[[2L]], , drop = FALSE]
  known <- data[data$method == "Ground truth", , drop = FALSE]
  fitted <- data[data$method == "GCMSprime", , drop = FALSE]
  ymax <- 1.06 * max(data$intensity)
  q <- quality[quality$kind == kind, , drop = FALSE]
  main <- if (kind == "Regular") "Regular EMG-like component" else
    "Shark-fin-shaped acid component"
  graphics::plot(NA, xlim = window, ylim = c(0, ymax),
    xlab = "Retention index", ylab = "Total-ion intensity", main = main,
    col.main = navy, font.main = 2, cex.main = 1.25)
  graphics::grid(col = grid_colour, lty = 1)
  draw_curves(fitted, 1, 1.8)
  draw_curves(known, 2, 1.35)
  graphics::mtext(sprintf("Mean profile cosine %.3f | %d/%d occurrences recovered",
    q$mean_profile_cosine, q$recovered, q$present),
    side = 3, line = 0.55, cex = 0.88, col = muted)
}
graphics::mtext(
  "Solid: inferred profiles | dashed: ground truth | blue: smallest true area | red: largest true area",
  side = 1, outer = TRUE, line = 0.65, cex = 0.98, col = muted
)
graphics::mtext("Representative profile recovery across the 16-sample simultaneous fit",
  side = 3, outer = TRUE, line = 0.4, cex = 1.35, font = 2, col = navy)
grDevices::dev.off()

invisible(TRUE)
