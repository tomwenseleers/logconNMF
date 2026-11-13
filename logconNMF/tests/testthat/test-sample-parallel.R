test_that("channel maps conserve counts and support saved reference maps", {
  set.seed(715);Y<-matrix(rpois(80*12,3),80,12,dimnames=list(NULL,as.character(40:51)))
  groups<-channelGroups(Y,4L,ncores=1L)
  f<-tempfile();saveRDS(groups,f);on.exit(unlink(f))
  Z<-aggregateChannels(Y,readRDS(f))
  expect_identical(unname(rowSums(Z)),unname(rowSums(Y)))
  expect_equal(ncol(Z),4L)
  expect_identical(aggregateChannels(Y[,12:1],groups),Z)
  unseen<-cbind(Y,`999`=rep(7,nrow(Y)))
  added<-aggregateChannels(unseen,groups)-Z
  expect_identical(added[,groups$unknown_group],rep(7,nrow(Y)))
  expect_equal(sum(added),7*nrow(Y))
})

test_that("native sample slices retain original rows and channel order", {
  Y<-matrix(as.double(seq_len(300)),30,10);path<-paste0(tempfile(),"_\u00e9.bin")
  on.exit(unlink(path));sinmf_write_matrix_file(Y,path)
  for(rows in list(4:12,c(19L,3L,5L)))for(channels in list(1:10,c(9L,2L)))
    expect_identical(sinmf_read_matrix_file(path,30,10,channels,rows),Y[rows,channels,drop=FALSE])
  expect_error(sinmf_read_matrix_file(path,30,10,11L),"Invalid matrix slice")
  con<-file(path,"wb");writeBin(1,con);close(con)
  expect_error(sinmf_read_matrix_file(path,30,10),"Truncated")
})

test_that("sample scheduling preserves fits for every supported noise family", {
  for(family in c("poisson","quasipoisson","gaussian","nb","nb1","Gamma","Gamma_scale")) {
    set.seed(422);n<-70L;m<-8L;H<-convMatrix(dnorm(-5:5,sd=1.5),n)[,c(18L,48L)]
    H<-H%*%Matrix::Diagonal(x=1/colSums(H));S<-matrix(runif(2*m,120,800),2,m)
    phi<-if(family=="nb1")2 else if(family=="gaussian")3 else .5
    if(family=="quasipoisson")phi<-2
    if(family=="poisson")phi<-1
    mu<-as.matrix(H%*%S)+2
    Y<-switch(family,poisson=matrix(rpois(n*m,mu),n,m),
      quasipoisson=matrix(rpois(n*m,mu/phi)*phi,n,m),gaussian=mu+matrix(rnorm(n*m,sd=sqrt(phi)),n,m),
      nb=matrix(rnbinom(n*m,mu=mu,size=20),n,m),nb1=matrix(rnbinom(n*m,mu=mu,size=mu/(phi-1)),n,m),
      Gamma=matrix(rgamma(n*m,shape=1/phi,scale=mu*phi),n,m),
      Gamma_scale=matrix(rgamma(n*m,shape=mu/phi,scale=phi),n,m))
    if(family%in%c("Gamma","Gamma_scale"))Y[Y<.5]<-0
    storage.mode(Y)<-"double"
    spec<-slcn_noise(Y,slcn_family(family,if(family=="nb")NULL else phi,
      if(family=="nb")20 else NULL),1,1L,"none",.5)$spec
    f<-list(H=H,S=S,background=matrix(1,n,1),alpha=matrix(2,1,m),spec=spec)
    desc<-lapply(1:2,function(s)list(y=Y,time=seq_len(n)+s*.2,fit=f,
      weights=seq(.8,1.2,length.out=n),offset=rep(0,n),sample_id=paste0("s",s),small_width=9))
    paths<-c(tempfile(),tempfile());on.exit(unlink(paths,recursive=TRUE),add=TRUE)
    ctl<-list(initial_cycles=1L,cycles=1L,rounds=1L,q=12L)
    one<-simultaneousNMF(desc,family=family,path=paths[1],precluster=3L,ncores=1L,control=ctl)
    two<-simultaneousNMF(desc,family=family,path=paths[2],precluster=one$channel_groups,ncores=2L,control=ctl)
    expect_equal(two$P,one$P,tolerance=1e-10,info=family)
    expect_equal(two$areas,one$areas,tolerance=1e-10,info=family)
    expect_equal(GIC(two),GIC(one),tolerance=1e-10,info=family)
    expect_equal(fitted(two),fitted(one),tolerance=1e-10,info=family)
    expect_null(two$runtime);expect_null(two$store$cache)
    expect_equal(two$parallel,list(ncores=2L,workers=2L,threads=1L))
    expect_true(all(two$history$delta[two$history$stage=="simultaneous_cycle"]<=1e-4))
    expect_equal(sinmf_sample_data(two$store,1L,aggregate=TRUE),aggregateChannels(Y,two$channel_groups),ignore_attr=TRUE)
  }
})

test_that("parallel initialization retains each sample's sparse PSF", {
  set.seed(234);n<-90L;m<-6L;S<-matrix(rexp(2*m,1/4000),2,m)
  desc<-lapply(c(1.4,2.2),function(sd) {
    X<-convMatrix(dnorm(-9:9,sd=sd),n)
    list(y=matrix(rpois(n*m,as.matrix(X[,c(25,65)]%*%S)+2),n,m),
      X=X,time=seq_len(n),small_width=10,background=matrix(1,n,1))
  })
  paths<-c(tempfile(),tempfile());on.exit(unlink(paths,recursive=TRUE))
  control<-list(cycles=0L,initial_cycles=0L,rounds=0L,
    initial=list(profile_cycles=0L,original_cleanup=FALSE))
  one<-simultaneousNMF(desc,path=paths[1],precluster=FALSE,ncores=1L,control=control)
  two<-simultaneousNMF(desc,path=paths[2],precluster=FALSE,ncores=2L,control=control)
  expect_equal(one$P,two$P,tolerance=1e-10)
  expect_equal(one$areas,two$areas,tolerance=1e-10)
  expect_equal(one$L,4L)
  expect_false(identical(readRDS(two$paths[1])$H,readRDS(two$paths[2])$H))
  expect_s4_class(readRDS(two$paths[1])$H,"dgCMatrix")
})

test_that("numeric AIC weight agrees with the named criterion during initialization", {
  set.seed(912);n<-60L;m<-4L;X<-convMatrix(dnorm(-5:5,sd=1.3),n)
  S<-matrix(runif(2*m,500,1500),2,m)
  desc<-lapply(1:2,function(s)list(y=matrix(rpois(n*m,as.matrix(X[,c(18,43)]%*%S)+1),n,m),
    X=X,time=seq_len(n),background=matrix(1,n,1)))
  paths<-c(tempfile(),tempfile());on.exit(unlink(paths,recursive=TRUE))
  control<-list(cycles=0L,initial_cycles=0L,rounds=0L,
    initial=list(profile_cycles=0L,original_cleanup=FALSE))
  named<-simultaneousNMF(desc,lambda='aic',path=paths[1],precluster=FALSE,ncores=1L,control=control)
  numeric<-simultaneousNMF(desc,lambda=2,path=paths[2],precluster=FALSE,ncores=1L,control=control)
  expect_equal(named$P,numeric$P,tolerance=1e-12)
  expect_equal(named$areas,numeric$areas,tolerance=1e-12)
  expect_equal(named$objective,numeric$objective,tolerance=1e-12)
})

test_that("incremental fusion candidates equal a complete proposal rebuild", {
  set.seed(731);m<-10L;K<-24L
  P<-matrix(runif(6*m,.1,1),6,m)[rep(1:6,4),]+matrix(runif(K*m,0,.01),K,m)
  P<-P/rowSums(P)
  state<-list(P=P,map=seq_len(K),scale=rep(1,K),channel_groups=NULL,
    control=list(rt_tolerance=1.1))
  occ<-data.frame(root=seq_len(K),sample=rep(1:4,each=6),
    anchor=rep(seq_len(6),4)+runif(K,-.04,.04),area=runif(K,1,5))
  for(method in c("cosine","hellinger")) {
    z<-state;o<-occ;queue<-sinmf_proposals(z,method,occ=o)
    for(step in 1:15) {
      if(!nrow(queue))break
      edge<-queue[1,];q<-queue[-1,,drop=FALSE]
      weight<-vapply(c(edge$a,edge$b),function(g)sum(o$area[o$root==g]),0)
      z$P[edge$a,]<-as.numeric(crossprod(weight,z$P[c(edge$a,edge$b),,drop=FALSE]))/sum(weight)
      z$map[z$map==edge$b]<-edge$a;o$root[o$root==edge$b]<-edge$a
      queue<-sinmf_refresh_proposals(z,method,o,z$P,q,c(edge$a,edge$b))
      expected<-sinmf_proposals(z,method,occ=o)
      rownames(queue)<-rownames(expected)<-NULL
      if(nrow(expected))expect_equal(queue,expected,tolerance=0)else expect_equal(nrow(queue),0L)
    }
  }
})
