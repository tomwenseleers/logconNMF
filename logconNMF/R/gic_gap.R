# Match reference columns across samples without using their areas or spectra.
cc_truth_profiles <- function(truth, times, samples) {
  if (inherits(truth,"GCMS_simulated_samples"))
    truth <- stats::setNames(truth$samples,vapply(truth$samples,"[[","","sample_id"))
  single <- is.matrix(truth) || inherits(truth, "Matrix") ||
    (length(times)==1L && is.list(truth) && any(vapply(truth[c("H","W_true","H_true")],
                                 function(h)!is.null(dim(h)), TRUE)))
  if (single) truth <- list(truth)
  if (!is.list(truth) || length(truth) != length(times))
    stop("Supply one truth.profiles matrix or profile record per sample.", call.=FALSE)
  if (!is.null(names(truth)) && length(times) > 1L) {
    if (anyDuplicated(names(truth)) || !setequal(names(truth), samples))
      stop("Names of truth.profiles must match the sample IDs.", call.=FALSE)
    truth <- truth[samples]
  }
  pieces <- Map(function(record, time) {
    if (is.matrix(record) || inherits(record, "Matrix")) record <- list(H=record)
    if (is.list(record) && is.null(record[["H"]]) && !is.null(record[["W_true"]]))
      record <- list(H=record$W_true)
    if (is.list(record) && is.null(record[["H"]]) && !is.null(record[["H_true"]])) {
      h <- record$H_true
      present <- record$present
      if (!is.logical(present) || length(present)!=ncol(h) || anyNA(present) ||
          length(record$exposure)!=nrow(h) || any(!is.finite(record$exposure)) || any(record$exposure<=0))
        stop("Invalid simulated occurrence or exposure in truth.profiles.",call.=FALSE)
      record <- list(H=h[,present,drop=FALSE]*record$exposure,id=which(present),time=record$time)
    }
    if (!is.list(record) || !(is.matrix(record$H) || inherits(record$H, "Matrix")) ||
        !(is.numeric(record$H) || inherits(record$H, "dMatrix")))
      stop("Each truth.profiles record needs a numeric scan-by-compound matrix H.", call.=FALSE)
    H <- methods::as(record$H, "dgCMatrix")
    if (any(!is.finite(H@x)) || any(H@x < 0))
      stop("truth.profiles must contain finite nonnegative heights.", call.=FALSE)
    if (!is.null(record$time)) {
      t <- record$time
      if (!is.numeric(t) || length(t) != nrow(H) || !length(t) || any(!is.finite(t)) ||
          is.unsorted(t, strictly=TRUE)) stop("Invalid scan times in truth.profiles.", call.=FALSE)
      i <- pmax(1L, findInterval(time, t))
      j <- pmin(i+1L, length(t))
      i <- ifelse(abs(t[j]-time) < abs(t[i]-time), j, i)
      if (anyDuplicated(i) || any(abs(t[i]-time) > 1e-7*pmax(1,abs(time))))
        stop("truth.profiles must cover the fitted scan times.", call.=FALSE)
      if (length(i) != nrow(H) || any(i != seq_along(i))) H <- H[i,,drop=FALSE]
    }
    if (nrow(H) != length(time)) stop("truth.profiles must have one row per fitted scan.",call.=FALSE)
    id <- if (is.null(record$id)) colnames(H) else record$id
    if (!is.null(id) && (!is.atomic(id) || length(id) != ncol(H) || anyNA(id) ||
                         any(!nzchar(as.character(id))) || anyDuplicated(as.character(id))))
      stop("Use unique compound IDs for the columns of each truth.profiles matrix.",call.=FALSE)
    list(H=H, id=id)
  }, unname(truth), times)
  labelled <- vapply(pieces,function(x)!is.null(x$id) || !ncol(x$H),TRUE)
  if (any(!labelled) && (any(labelled & vapply(pieces,function(x)ncol(x$H)>0L,TRUE)) ||
                          length(unique(vapply(pieces,function(x)ncol(x$H),0L)))>1L))
    stop("Supply shared compound IDs or use the same columns in every reference matrix.",call.=FALSE)
  lapply(pieces,function(x) {
    id <- if (is.null(x$id)) seq_len(ncol(x$H)) else x$id
    size <- as.numeric(Matrix::colSums(x$H))
    if (any(!is.finite(size))) stop("Non-finite total in truth.profiles.",call.=FALSE)
    live <- which(size>0)
    h <- x$H[,live,drop=FALSE]
    if (length(live)) h <- h %*% Matrix::Diagonal(x=1/size[live])
    list(H=methods::as(h,"dgCMatrix"),id=id[live])
  })
}

# Bin integration and a nonlinear time map can make the supplied counts
# nonconcave on the observed grid. The upper hull of (time, log height)
# gives the least concave majorant on the same support interval.
cc_feasible_reference <- function(H, time) {
  if (!ncol(H)) return(H)
  rows <- values <- vector("list",ncol(H))
  for (k in seq_len(ncol(H))) {
    at <- H@p[k]+seq_len(H@p[k+1L]-H@p[k])
    at <- at[H@x[at] >= .Machine$double.xmin]
    i <- H@i[at]+1L
    v <- H@x[at]
    if (length(i)>1L) {
      z <- log(v)
      hull <- grDevices::chull(time[i],z)
      # chull returns clockwise vertices. From left to right this is the upper hull.
      upper <- rep(hull,2L)[match(1L,hull)+seq_along(hull)-1L]
      upper <- upper[seq_len(match(length(i),upper))]
      if (length(upper)<length(i) || any(diff(i)>1L)) {
        grid <- seq.int(i[1L],i[length(i)])
        v <- exp(stats::approx(time[i[upper]],z[upper],time[grid],ties="ordered")$y)
        i <- grid
      }
    }
    rows[[k]] <- i
    values[[k]] <- v/sum(v)
  }
  Matrix::sparseMatrix(i=unlist(rows,use.names=FALSE),p=c(0L,cumsum(lengths(rows))),
                       x=unlist(values,use.names=FALSE),dims=dim(H),dimnames=dimnames(H))
}

# Ion coefficients and background are estimated from observed counts. An
# arbitrary positive rescaling of one input column therefore carries no signal.
cc_reference_sample <- function(part, ref, group, price, threads) {
  h <- cc_feasible_reference(ref$H,part$time)
  y <- if (is.character(part$observed))
    sinmf_read_matrix_file(part$observed,nrow(h),length(group)) else as.matrix(part$observed)
  B <- methods::as(unclass(as.matrix(part$background)), "dgCMatrix")
  coef <- cc_weighted_nnls(cbind(h,B),y,part$spec$dispersion)
  part$H <- h
  part$S <- coef[seq_len(ncol(h)),,drop=FALSE]
  part$alpha <- coef[ncol(h)+seq_len(ncol(B)),,drop=FALSE]
  part$observed <- y
  part$windows <- conditional_windows(h)
  part$grouped_dictionary <- FALSE
  x <- cc_sample(part,group,price,threads)
  x$id <- ref$id
  x$norms <- rowSums(x$initial_p)
  x
}

cc_reference_dictionary <- function(data) {
  labels <- sort(unique(unlist(lapply(data,function(x)x$id),use.names=FALSE)))
  K <- length(labels)
  p <- matrix(0,K,ncol(data[[1L]]$y))
  bands <- cbind(rep(Inf,K),rep(-Inf,K))
  for (x in data) {
    k <- match(x$id,labels)
    x$h <- matrix(0,length(x$t),K)
    if (!length(k)) next
    x$h[,k] <- as.matrix(x$init %*% Matrix::Diagonal(x=x$norms))
    p[k,] <- p[k,,drop=FALSE] + x$initial_p
    bands[k,1] <- pmin(bands[k,1],x$warp[x$grid[,1]])
    bands[k,2] <- pmax(bands[k,2],x$warp[x$grid[,2]])
  }
  live <- which(rowSums(p)>0)
  if (length(live)!=K) for (x in data) x$h <- x$h[,live,drop=FALSE]
  p <- p[live,,drop=FALSE]
  list(p=p/pmax(rowSums(p),.Machine$double.xmin),bands=bands[live,,drop=FALSE])
}

cc_gap_traces <- function(result, price) {
  trace <- rbind(data.frame(start="standard",result$trace),
                 data.frame(start="true profiles",result$reference_trace))
  trace$sweep <- stats::ave(as.numeric(trace$block=="spectra and background"),
                            trace$start,FUN=cumsum) + 0.5*(trace$block=="profiles")
  trace$GIC <- trace$D + price*(trace$spectral_df + trace$shape_df)
  structure(list(trace=trace,reference_seconds=result$reference_seconds),class="logconNMF_gap")
}

#' Empirical GIC gap from matched reference profiles
#'
#' The second fit starts from unit-area profiles on the observed scan grid.
#' A least concave majorant of log height gives a feasible starting shape
#' when irregular acquisition or bin integration alters log-concavity.
#' Both fits estimate spectra, amplitudes, shapes and background from the
#' same data. The block traces include the public GIC at each checkpoint.
#' With a different selection penalty, its decrease need not imply a GIC decrease.
#'
#' @md
#' @param object A logconNMF or GCMSprime fit computed with `truth.profiles`.
#' @return A plottable object with the standard and reference GIC, their signed
#'   `difference`, `percent = 100 * difference / abs(reference_gic)`, both
#'   block traces and elapsed seconds for reference initialization and fitting.
#'   A zero reference GIC gives an undefined percentage. A negative difference
#'   means the standard fit reached a lower GIC. The reference is a feasible
#'   fitted solution. A global optimality bound would require a separate proof.
#'   An independently fitted GCMSprime directory returns a named list of these
#'   diagnostics, one per sample.
#' @example examples/gicGap.R
#' @seealso [logconNMF()]
#' @export
gicGap <- function(object) {
  x <- object$gic_gap
  if (!inherits(x, "logconNMF_gap")) {
    if (is.list(x) && length(x) && all(vapply(x, inherits, TRUE, "logconNMF_gap")))
      return(x)
    stop("Fit with truth.profiles to obtain a GIC diagnostic.", call.=FALSE)
  }
  last <- function(which) utils::tail(x$trace$GIC[x$trace$start==which],1L)
  x$fitted_gic <- as.numeric(GIC(object))
  x$trace$GIC <- x$trace$GIC + (x$fitted_gic-last("standard"))
  x$reference_gic <- last("true profiles")
  x$difference <- x$fitted_gic-x$reference_gic
  x$percent <- if(x$reference_gic==0) NA_real_ else 100*x$difference/abs(x$reference_gic)
  x
}

#' @rdname gicGap
#' @param x The returned diagnostic.
#' @param type Show both panels, the complete GIC trace, or the final half of
#'   each trajectory on a percentage scale relative to the reference endpoint.
#' @param ... Other arguments passed to the base plot.
#' @export
plot.logconNMF_gap <- function(x,type=c("both","GIC","gap"),...) {
  type <- match.arg(type)
  t <- x$trace
  arms <- list(standard=t[t$start=="standard",],reference=t[t$start=="true profiles",])
  colors <- c("#17649a","#aa4c13")
  if (type=="both") {
    old <- graphics::par(mfrow=c(1,2))
    on.exit(graphics::par(old))
  }
  for (panel in if (type=="both") c("GIC","gap") else type) {
    series <- arms
    label <- "GIC"
    if (panel=="gap") {
      percent <- x$reference_gic!=0
      label <- if (percent) "GIC difference from reference (%)" else "GIC difference from reference"
      series <- lapply(series,function(z) {
        z <- z[z$sweep>=max(z$sweep)/2,,drop=FALSE]
        z$GIC <- (z$GIC-x$reference_gic)*if (percent) 100/abs(x$reference_gic) else 1
        z
      })
    }
    both <- rbind(series[[1L]],series[[2L]])
    title <- if (panel=="GIC") "GIC across updates" else if (is.finite(x$percent))
      sprintf("Final gap %.3g%%",x$percent) else "Final sweeps"
    args <- utils::modifyList(list(x=range(both$sweep),y=range(both$GIC),type="n",
                                   main=title,las=1,xlab="Profile and spectrum sweeps",ylab=label),list(...))
    do.call(graphics::plot,args)
    if (panel=="gap") graphics::abline(h=0,lty=3,col="grey60")
    for(i in seq_along(series)) {
      z <- series[[i]]
      graphics::lines(z$sweep,z$GIC,col=colors[i],lwd=2,lty=i)
      graphics::points(utils::tail(z$sweep,1),utils::tail(z$GIC,1),pch=19,col=colors[i])
    }
    if (panel=="GIC") graphics::legend("topright",c("Standard fit","True profiles and matches"),
                                      col=colors,lty=1:2,lwd=2,bty="n",cex=.8)
  }
  invisible(x)
}

#' @rdname gicGap
#' @method print logconNMF_gap
#' @export
print.logconNMF_gap <- function(x,...) {
  cat(sprintf("Standard GIC %.8g, reference GIC %.8g, empirical gap %.4g%%\n",
              x$fitted_gic,x$reference_gic,x$percent))
  invisible(x)
}
