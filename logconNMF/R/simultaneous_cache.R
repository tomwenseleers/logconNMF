# One entry belongs to an active coupled runtime, never to a saved model.
.sinmf_cache <- new.env(parent=emptyenv())
sinmf_cache_reset <- function(enabled=FALSE) {
  rm(list=ls(.sinmf_cache,all.names=TRUE),envir=.sinmf_cache)
  .sinmf_cache$enabled<-isTRUE(enabled);.sinmf_cache$eligible<-FALSE
  .sinmf_cache$key<-NULL;.sinmf_cache$Y<-NULL
  invisible(NULL)
}
sinmf_cache_reset()
sinmf_cache_restore <- function(previous) {
  sinmf_cache_reset();list2env(previous,envir=.sinmf_cache)
  invisible(NULL)
}
sinmf_cache_sample <- function(kind) {
  .sinmf_cache$eligible<-isTRUE(.sinmf_cache$enabled)&&kind%in%1:2
  if(!.sinmf_cache$eligible){.sinmf_cache$key<-NULL;.sinmf_cache$Y<-NULL}
  invisible(NULL)
}
sinmf_cached_read <- function(path,n,m,channels,rows) {
  # Original observations are immutable during a runtime. Larger samples and
  # other families retain the streaming reader; the normalizer doubles this cap.
  if(!isTRUE(.sinmf_cache$eligible)||basename(path)!="observed.bin"||
     8*as.double(n)*m>64*1024^2)return(NULL)
  key<-paste(path,n,m,sep=";")
  if(!identical(.sinmf_cache$key,key)) {
    Y<-slcn_read_block(path,as.integer(n),seq_len(m),seq_len(n))
    attr(Y,"sinmf_log_factorial")<-slcn_log_factorial(Y,1L)
    .sinmf_cache$Y<-Y;.sinmf_cache$key<-key
  }
  if(identical(channels,seq_len(m))&&identical(rows,seq_len(n)))return(.sinmf_cache$Y)
  Y<-.sinmf_cache$Y[rows,channels,drop=FALSE]
  attr(Y,"sinmf_log_factorial")<-attr(.sinmf_cache$Y,"sinmf_log_factorial")[rows,channels,drop=FALSE]
  Y
}
