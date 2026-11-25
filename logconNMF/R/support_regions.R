# Connected intervals in each sample, joined across samples by shared identity.
cc_support_region <- function(ranges, seeds) {
  ids <- seeds
  groups <- lapply(ranges, function(b) {
    live <- which(is.finite(b[, 1L]))
    if (!length(live)) return(list())
    ordered <- live[order(b[live, 1L])]
    begins <- b[ordered, 1L] > c(-Inf, utils::head(cummax(b[ordered, 2L]), -1L))
    split(ordered, cumsum(begins))
  })
  repeat {
    previous <- length(ids)
    for (sample in groups) for (group in sample) {
      if (any(group %in% ids)) ids <- union(ids, group)
    }
    if (length(ids) == previous) return(sort(ids))
  }
}

cc_region_state <- function(z, ids) {
  a <- new.env(parent = emptyenv())
  a$p <- z$p[ids, , drop = FALSE]
  rownames(a$p) <- as.character(ids)
  a$bands <- cbind(rep(Inf, length(ids)), rep(-Inf, length(ids)))
  a$price <- z$price
  a$threads <- z$threads
  a$trace <- FALSE
  a$history <- list()
  other <- setdiff(seq_len(nrow(z$p)), ids)
  a$data <- lapply(z$data, function(x) {
    s <- new.env(parent = emptyenv())
    s$t <- x$t
    s$warp <- x$warp
    s$phi <- x$phi
    s$y <- x$y
    s$z <- x$z
    s$alpha <- x$alpha + 0
    h <- as(x$h, "dgCMatrix")
    live <- which(diff(h@p) > 0L)
    s$grid <- if (length(live)) conditional_windows(h[, live, drop = FALSE]) else matrix(0L, 0L, 2L)
    s$width <- if (length(x$t) > 1L) max(cc_width(s), 3 * median(diff(s$warp))) else x$width
    ix <- which(ids %in% live)
    if (length(ix)) {
      bounds <- s$grid[match(ids[ix], live), , drop = FALSE]
      a$bands[ix, 1L] <- pmin(a$bands[ix, 1L], s$warp[bounds[, 1L]])
      a$bands[ix, 2L] <- pmax(a$bands[ix, 2L], s$warp[bounds[, 2L]])
    }
    s$h <- x$h[, ids, drop = FALSE]
    s$fixed <- if (length(other)) as.matrix(h[, other, drop = FALSE] %*% z$p[other, , drop = FALSE]) else 0
    s$bg <- pmax(s$z %*% s$alpha + 1e-8, 1e-8) + s$fixed
    cc_means(s, a$p)
    s
  })
  a
}

cc_refine_region <- function(z, ids, selective) {
  a <- cc_region_state(z, ids)
  score <- function() cc_dimension(a$data, a$p, cc_D(a$data), a$price)
  input <- score()
  support <- function() {
    present <- cc_mass(a$data) > 0
    colnames(present) <- rownames(a$p)
    present
  }
  previous <- support()
  rounds <- 0L
  for (i in seq_len(4L)) {
    if (!nrow(a$p)) break
    cc_end_refine(a, list(refinement_sweeps = if (i == 1L) 10L else 0L,
                          penalty_search = TRUE, selective = selective))
    rounds <- i
    current <- score()
    present <- support()
    if (!nrow(a$p) || identical(present, previous)) break
    previous <- present
  }
  if (!rounds || !is.finite(current$objective) || current$objective >= input$objective)
    return(rounds)
  live <- as.integer(rownames(a$p))
  z$p[ids, ] <- 0
  z$p[live, ] <- a$p
  z$bands[live, ] <- a$bands
  for (s in seq_along(z$data)) {
    x <- z$data[[s]]
    y <- a$data[[s]]
    x$h[, ids] <- 0
    x$h[, live] <- y$h
    x$alpha <- y$alpha
    x$bg <- pmax(x$z %*% x$alpha + 1e-8, 1e-8)
    cc_means(x, z$p)
  }
  rounds
}

# A split identity can occur in different samples, where deleting either
# occurrence would throw away real signal. Test one shared spectrum while
# allowing overlapping neighbours to account for contaminating ions.
cc_share_nearby_spectra <- function(z) {
  rank <- nrow(z$p)
  result <- c(tested = 0L, accepted = 0L)
  if (length(z$data) < 2L || rank < 2L) return(result)
  bands <- cbind(rep(Inf, rank), rep(-Inf, rank))
  ranges <- vector("list", length(z$data))
  for (s in seq_along(z$data)) {
    x <- z$data[[s]]
    h <- as(x$h, "dgCMatrix")
    for (k in which(diff(h@p) > 0L)) {
      at <- seq.int(h@p[k] + 1L, h@p[k + 1L])
      height <- h@x[at]
      ix <- h@i[at[height >= .5 * max(height)]] + 1L
      bands[k, 1L] <- min(bands[k, 1L], x$warp[ix[1L]])
      bands[k, 2L] <- max(bands[k, 2L], x$warp[ix[length(ix)]])
    }
    ranges[[s]] <- cc_profile_ranges(h)
  }
  unit <- cc_unit(z$p)
  for (k in seq_len(rank)) {
    if (!any(z$p[k, ] > 0)) next
    near <- which(seq_len(rank) < k & bands[, 1L] <= bands[k, 2L] & bands[, 2L] >= bands[k, 1L])
    if (!length(near)) next
    # The existing local spectral exchange uses cosine > 0.9 for a
    # single-spectral-direction proposal. Only the GIC can accept one.
    partners <- near[drop(unit[near, , drop = FALSE] %*% unit[k, ]) > .9]
    for (j in partners) {
      if (!any(z$p[k, ] > 0) || !any(z$p[j, ] > 0)) next
      pair <- c(k, j)
      overlap <- unique(unlist(lapply(ranges, function(b)
        unlist(lapply(pair[is.finite(b[pair, 1L])], function(i)
          which(b[, 1L] <= b[i, 2L] & b[, 2L] >= b[i, 1L]))))))
      others <- setdiff(overlap, pair)
      result["tested"] <- result["tested"] + 1L
      best <- NULL
      for (remove in pair) {
        ids <- c(remove, setdiff(pair, remove), others)
        q <- z$p[ids, , drop = FALSE]
        weights <- cc_cone(t(q[-1L, , drop = FALSE]), q[1L, ])
        blocks <- cc_blocks(z, ids)
        fit <- cc_mixture_trial(blocks, q, weights, z$price, FALSE)
        if (!is.null(fit) && (is.null(best) || fit$cost < best$fit$cost))
          best <- list(fit = fit, blocks = blocks, ids = ids)
      }
      if (!is.null(best)) {
        cc_commit_block(z, best$blocks, best$ids, best$fit)
        for (b in best$blocks) {
          s <- b$si
          ranges[[s]][best$ids, ] <- cc_profile_ranges(z$data[[s]]$h[, best$ids, drop = FALSE])
        }
        unit[best$ids, ] <- cc_unit(z$p[best$ids, , drop = FALSE])
        result["accepted"] <- result["accepted"] + 1L
      }
    }
  }
  result
}

cc_local_search <- function(z, before, selective) {
  screened <- removed <- region <- rounds <- 0L
  pairs <- c(tested = 0L, accepted = 0L)
  pass <- function(limit) {
    for (i in seq_len(limit)) {
      scan <- cc_screen_support(z)
      screened <<- screened + sum(scan$flagged)
      count <- cc_try_transfers(z, scan)
      if (!count) return(scan)
      removed <<- removed + count
    }
    cc_screen_support(z)
  }
  if (z$price > 0 && nrow(z$p)) {
    scan <- pass(4L)
    if (any(scan$shared)) {
      ids <- cc_support_region(scan$ranges, which(scan$shared))
      region <- length(ids)
      rounds <- cc_refine_region(z, ids, selective)
      pass(2L)
    }
    cc_trim(z)
    pairs <- cc_share_nearby_spectra(z)
  }
  cc_trim(z)
  stat <- cc_dimension(z$data, z$p, cc_D(z$data), z$price)
  if (stat$objective > before + 1e-5 + 1e-8 * abs(before))
    stop("Support search increased the original-ion objective.")
  if (z$trace) cc_history(z, "local support search", 1L, "original price", stat$D)
  list(before = before, after = stat$objective, retained = "local support search",
       screened = screened, removed_by_transfer = removed, regional_factors = region,
       regional_rounds = rounds, identity_pairs = pairs)
}
