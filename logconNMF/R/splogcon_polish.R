# Score the factors themselves, including every original response channel.
slcn_sync <- function(f, Y) {
  f$mu <- as.matrix(f$H %*% f$S) + f$background %*% f$alpha + f$offset
  f$rank <- ncol(f$H)
  f$loglik <- slcn_loglik(Y, f$mu, f$spec$kind, f$spec$dispersion,
    f$spec$theta, as.numeric(f$weights), as.integer(f$ncores), f$control$censor)
  f$objective <- -2 * f$loglik + f$lambda * f$rank
  f$ic_df <- ncol(Y) * f$rank + if (is.null(f$nuisance_df)) 0L else f$nuisance_df
  f
}

slcn_contained_pairs <- function(f, broad = FALSE, psf_width = f$control$psf_width) {
  K <- ncol(f$H)
  empty <- data.frame(fragment = integer(), anchor = integer(), similarity = numeric())
  if (K < 2L) return(empty)
  H <- f$H;S <- f$S;area <- rowSums(S);geometry <- slcn_profile_geometry(f)
  cosine <- tcrossprod(S) / sqrt(rowSums(S^2) %o% rowSums(S^2))
  if (!broad) hc <- as.matrix(crossprod(H)) / sqrt(colSums(H^2) %o% colSums(H^2))
  parents <- if (broad) which(geometry$width > 4 * psf_width) else seq_len(K)
  out <- list()
  for (b in parents) {
    h <- as.numeric(H[, b]);at <- which(h >= max(h) * .01)
    inside <- colSums(H[min(at):max(at), , drop = FALSE])
    pool <- which(seq_len(K) != b & is.finite(cosine[, b]) &
      inside >= if (broad) .99 else .98)
    pool <- pool[cosine[pool, b] >= if (broad) .97 else .9]
    pool <- pool[area[pool] <= area[b] * if (broad) .15 else .3]
    if (!broad) pool <- pool[hc[pool, b] >= .8]
    if (length(pool)) out[[length(out) + 1L]] <- data.frame(fragment = pool,
      anchor = b, similarity = cosine[pool, b])
  }
  if (!length(out)) return(empty)
  out <- do.call(rbind, out)
  out[order(out$similarity, decreasing = TRUE), , drop = FALSE]
}

slcn_local_block <- function(f, group, anchor, threshold = .01) {
  h <- as.numeric(f$H[, anchor]);at <- which(h > max(h) * 1e-5)
  sort(unique(c(group, which(colSums(f$H[at, , drop = FALSE]) > threshold))))
}

# Lineage preserves proposals when joint refinement redistributes loadings.
slcn_retained_pairs <- function(f, seed) {
  proposals <- slcn_contained_pairs(f)
  original <- slcn_contained_pairs(seed)
  locate <- function(ids) which(vapply(f$members, function(x) any(ids %in% x), TRUE))
  for (i in seq_len(nrow(original))) {
    a <- locate(seed$members[[original$fragment[i]]])
    b <- locate(seed$members[[original$anchor[i]]])
    if (length(a) == 1L && length(b) == 1L && a != b)
      proposals <- rbind(proposals, data.frame(fragment = a, anchor = b,
        similarity = original$similarity[i]))
  }
  proposals[!duplicated(proposals[c("fragment", "anchor")]), , drop = FALSE]
}

slcn_contained_narrow <- function(f, Y, seed, verbose = FALSE) {
  audit <- list()
  for (pass in 1:6) {
    proposals <- head(slcn_retained_pairs(f, seed), 12L)
    if (!nrow(proposals)) break
    best <- NULL;value <- f$objective - 1e-6
    for (i in seq_len(nrow(proposals))) for (shape in c("anchor", "hull")) {
      group <- c(proposals$anchor[i], proposals$fragment[i])
      proposal <- list(group = group, anchor = group[1L], shape = shape)
      seconds <- slcn_time(z <- slcn_merge_anchor(f, Y, proposal,
        cycles = 6L, profile_steps = 2L, neighbor_mass = .05))[["elapsed"]]
      z <- slcn_sync(z, Y)
      audit[[length(audit) + 1L]] <- data.frame(
        source = paste(sort(unlist(f$members[group])), collapse = ";"),
        shape = shape, seconds = seconds, delta = z$objective - f$objective,
        accepted = FALSE)
      if (is.finite(z$objective) && z$objective < value) {
        best <- z;value <- z$objective;chosen <- length(audit)
      }
    }
    if (is.null(best)) break
    audit[[chosen]]$accepted <- TRUE;f <- best
    if (verbose) {
      cat("contained fragments:", ncol(f$H), "components\n");flush.console()
    }
  }
  f$contained_fusion$narrow <- if (length(audit)) do.call(rbind, audit) else data.frame()
  f
}

slcn_contained_broad <- function(f, Y, psf_width, verbose = FALSE) {
  proposals <- slcn_contained_pairs(f, broad = TRUE, psf_width = psf_width)
  parents <- unique(proposals$anchor)
  if (length(parents)) parents <- parents[
    order(slcn_profile_geometry(f)$width[parents], decreasing = TRUE)]
  references <- lapply(parents, function(b) list(anchor = f$members[[b]],
    donors = f$members[proposals$fragment[proposals$anchor == b]]))
  audit <- list()
  for (reference in references) {
    locate <- function(ids) which(vapply(f$members, function(x) any(ids %in% x), TRUE))
    anchor <- locate(reference$anchor)
    if (length(anchor) != 1L) next
    donors <- unique(unlist(lapply(reference$donors, locate)))
    group <- unique(c(anchor, donors))
    if (length(group) < 2L) next
    members <- sort(unlist(f$members[group]))
    block <- slcn_local_block(f, group, anchor)
    seconds_control <- slcn_time(control <- slcn_refine_block(f, Y,
      block = block, cycles = 6L, profile_steps = 4L, adaptive = FALSE,
      free = block, knots = 256L))[["elapsed"]]
    control <- slcn_sync(control, Y)
    best <- if (is.finite(control$objective) && control$objective <= f$objective) control else f
    chosen <- NULL
    for (shape in c("anchor", "hull")) {
      seconds <- slcn_time({
        z <- slcn_merge_anchor(f, Y, list(group = group, anchor = anchor, shape = shape),
          cycles = 0L, neighbor_mass = .05)
        merged <- ncol(z$H);block <- slcn_local_block(z, merged, merged)
        z <- slcn_refine_block(z, Y, block = block, cycles = 6L, profile_steps = 4L,
          adaptive = FALSE, free = block, knots = 256L)
      })[["elapsed"]]
      z <- slcn_sync(z, Y)
      audit[[length(audit) + 1L]] <- data.frame(source = paste(members, collapse = ";"),
        shape = shape, seconds = seconds, control_seconds = seconds_control,
        delta = z$objective - control$objective, accepted = FALSE)
      if (is.finite(z$objective) && z$objective < best$objective - 1e-6) {
        best <- z;chosen <- length(audit)
      }
    }
    if (!is.null(chosen)) audit[[chosen]]$accepted <- TRUE
    f <- best
    if (verbose) {
      cat("broad contained fragments:", ncol(f$H), "components\n");flush.console()
    }
  }
  f$contained_fusion$broad <- if (length(audit)) do.call(rbind, audit) else data.frame()
  slcn_sync(f, Y)
}

slcn_polish <- function(f, Y, X, verbose = FALSE) {
  if (is.null(f$members)) f$members <- lapply(seq_len(ncol(f$H)), identity)
  f <- slcn_sync(f, Y);seed <- f
  advance <- function(label, expression) {
    seconds <- slcn_time(candidate <- force(expression))[["elapsed"]]
    candidate <- slcn_sync(candidate, Y)
    delta <- candidate$objective - f$objective
    ok <- is.finite(delta) && delta <= 1e-6
    history <- f$history
    if (ok) f <<- candidate
    f$history <<- rbind(history, data.frame(stage = label, seconds = seconds,
      delta = delta, accepted = ok, compounds = ncol(candidate$H)))
    if (verbose) {
      cat(label, round(seconds, 2), "s; change", signif(delta, 7),
        if (ok) "accepted;" else "rejected;", ncol(f$H), "components\n")
      flush.console()
    }
  }
  if (ncol(f$H) && f$control$profile_cycles)
    advance("fine_profiles", slcn_refine_block(f, Y, cycles = 4L, profile_steps = 2L,
      adaptive = FALSE, free = seq_len(ncol(f$H))))
  if (f$control$original_cleanup) advance("original_cleanup", slcn_cleanup(f, Y, X))
  if (ncol(f$H) && f$control$profile_cycles)
    advance("final_fine_profiles", slcn_refine_block(f, Y, cycles = 2L, profile_steps = 2L,
      adaptive = FALSE, free = seq_len(ncol(f$H))))
  if (f$control$original_cleanup) advance("final_prune", slcn_cleanup_prune(f, Y)$fit)
  if (ncol(f$H) > 1L && f$control$profile_fusion && f$control$profile_cycles)
    advance("contained_narrow_fusion", slcn_contained_narrow(f, Y, seed, verbose))
  f
}
