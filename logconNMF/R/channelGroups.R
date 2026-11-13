#' Learn and reuse a common channel aggregation
#'
#' Cluster outcome channels using size-adjusted Poisson distances and Ward
#' linkage. The resulting map can be reused across samples or learned from
#' reference spectra supplied as rows with `normalize=TRUE`.
#' @param y Numeric matrix with channels in columns.
#' @param k Maximum number of channel groups.
#' @param channels Unique channel identifiers, normally m/z values.
#' @param normalize Give every input row unit total, useful for spectral libraries.
#' @param ncores Native thread count, defaulting to physical cores.
#' @param groups A map returned by `channelGroups()`.
#' @return `channelGroups()` returns a compact reusable map accepted as
#'   `precluster` by [logconNMF()]. `aggregateChannels()`
#'   returns a matrix whose columns sum the mapped channels. Unseen channels
#'   are collected in the map's reserved group. Save maps with `saveRDS()`.
#' @examples
#' y <- matrix(rpois(100*12, 3), 100, 12)
#' groups <- channelGroups(y, k=4L, ncores=1L)
#' z <- aggregateChannels(y, groups)
#' stopifnot(sum(y) == sum(z))
#' @md
#' @export
channelGroups <- function(y,k=32L,channels=colnames(y),normalize=FALSE,ncores=NULL) {
  Y<-as.matrix(y);m<-ncol(Y)
  if(is.null(channels))channels<-as.character(seq_len(m))
  if(!is.numeric(Y)||!m||!nrow(Y)||any(!is.finite(Y)|Y<0)||
      length(channels)!=m||anyNA(channels)||anyDuplicated(channels))stop("Use finite nonnegative data and unique channel identifiers.")
  if(length(k)!=1L||!is.finite(k)||k<1||k!=floor(k))stop("k must be a positive integer.")
  if(!is.logical(normalize)||length(normalize)!=1L||is.na(normalize))stop("normalize must be TRUE or FALSE.")
  if(normalize)Y<-Y/pmax(rowSums(Y),1e-30)
  k<-min(as.integer(k),m);live<-colSums(Y)>0
  group<-rep.int(k,m)
  if(k==1L)group[]<-1L else if(any(live)) {
    count<-min(k-1L,sum(live))
    group[live]<-if(count==1L)1L else glmCluster(Y[,live,drop=FALSE],
      family=poisson("identity"),method="poissondist",size_fac=colSums(Y[,live,drop=FALSE]),
      k=count,nthreads=sinmf_core_count(ncores))$clusterids
  }
  structure(list(channels=as.character(channels),group=as.integer(group),k=k,
    unknown_group=k,normalize=normalize,method="size-adjusted Poisson Ward",version=1L),
    class="nmf_channel_groups")
}

sinmf_channel_index <- function(groups,channels) {
  if(!inherits(groups,"nmf_channel_groups")||length(groups$group)!=length(groups$channels)||
      length(groups$k)!=1L||anyNA(groups$group)||any(groups$group<1L|groups$group>groups$k))
    stop("Use a channelGroups() map.")
  z<-groups$group[match(as.character(channels),groups$channels)]
  z[is.na(z)]<-groups$unknown_group
  z
}

#' @rdname channelGroups
#' @export
aggregateChannels <- function(y,groups,channels=colnames(y)) {
  if(is.null(channels))channels<-as.character(seq_len(ncol(y)))
  if(length(channels)!=ncol(y))stop("channels must identify every column.")
  group<-sinmf_channel_index(groups,channels)
  out<-as.matrix(y%*%sparseMatrix(i=seq_along(group),j=group,x=1,dims=c(length(group),groups$k)))
  colnames(out)<-paste0("group",seq_len(groups$k));out
}

# Learn from summed time-bin counts without retaining original matrices.
sinmf_common_groups <- function(samples,precluster,ncores) {
  if(inherits(precluster,"nmf_channel_groups"))return(precluster)
  if(is.null(precluster)||identical(precluster,FALSE))return(NULL)
  k<-if(identical(precluster,TRUE))32L else precluster
  if(length(k)!=1L||!is.finite(k)||k<1||k!=floor(k))stop("Invalid precluster.")
  first<-samples[[1L]]$y;if(is.character(first))first<-readRDS(first)
  if(k>=ncol(first))return(NULL)
  domain<-range(unlist(lapply(samples,function(x)range(x$time))))
  pool<-NULL;axis<-NULL;bins<-2048L
  for(x in samples) {
    Y<-if(is.character(x$y))readRDS(x$y)else as.matrix(x$y)
    mass<-colnames(Y);if(is.null(mass))mass<-as.character(seq_len(ncol(Y)))
    if(is.null(axis)){axis<-mass;pool<-matrix(0,bins,ncol(Y))}
    if(!identical(axis,mass))stop("Sample channel names and ordering must agree.")
    at<-pmax(1L,pmin(bins,1L+floor((x$time-domain[1L])/max(diff(domain),1e-20)*bins)))
    z<-rowsum(pmax(Y,0),at,reorder=FALSE);ix<-as.integer(rownames(z));pool[ix,]<-pool[ix,,drop=FALSE]+z
  }
  channelGroups(pool,k,axis,ncores=ncores)
}
