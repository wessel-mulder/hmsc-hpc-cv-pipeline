# Numerical and sampler checks; no full analysis or simulation files are written.
root<-normalizePath(if(dir.exists('data-new')) '.' else '../../..')
out<-file.path(root,'scripts/1-atlas-analyses/homogenization')
source(file.path(out,'null-functions.R'))
Rcpp::sourceCpp(file.path(out,'homogenization-null.cpp'))
ref_beta<-function(x) {
  a<-tcrossprod(x); u<-outer(rowSums(x),rowSums(x),'+')-a
  d<-1-a/u;mean(d[lower.tri(d)])
}
as_int<-function(x) {storage.mode(x)<-'integer';x}
for(x in list(matrix(1L,4,3),diag(1L,4),as_int(matrix(c(1,0,1,1,1,0,0,1,1),3,byrow=TRUE))))
  stopifnot(abs(null_spatial_beta(x)-ref_beta(x))<1e-12)
# Anchor species prevents empty cells while a second species has 4 possible allocations.
x<-cbind(rep(1L,4),c(1L,1L,0L,0L));y<-cbind(rep(1L,4),c(1L,0L,1L,0L))
patterns<-character(4000)
for(i in seq_along(patterns)) {
  z<-simulate_spatial_null(x,y,1L,i,TRUE)
  m<-z$debug_last
  stopifnot(all(m %in% 0:1),identical(colSums(m),colSums(y)),
    all(colSums(x==1 & m==0)==colSums(x==1 & y==0)),
    all(colSums(x==0 & m==1)==colSums(x==0 & y==1)),
    abs(z$beta_later-ref_beta(m))<1e-12)
  patterns[i]<-paste(m[,2],collapse='')
}
counts<-table(patterns);stopifnot(length(counts)==4,all(abs(counts-1000)<150))
# Reproducibility and extreme feasible transition counts.
a<-simulate_spatial_null(x,y,7L,58L); b<-simulate_spatial_null(x,y,7L,58L)
stopifnot(identical(a,b))
for(endpoint in list(x,cbind(rep(1L,4),rep(0L,4)),cbind(rep(1L,4),rep(1L,4)),
                    cbind(rep(1L,4),1L-x[,2]))) {
  z<-simulate_spatial_null(x,endpoint,10L,82L,TRUE)
  stopifnot(identical(colSums(z$debug_last),colSums(endpoint)))
}
x_empty<-matrix(c(1L,0L,1L,0L,0L,1L),3,byrow=TRUE)
y_empty<-matrix(c(1L,0L,0L,1L,0L,1L),3,byrow=TRUE)
errors<-vapply(1:30,function(s) tryCatch({simulate_spatial_null(x_empty,y_empty,1L,s);FALSE},
  error=function(e) grepl('Empty simulated community',conditionMessage(e))),logical(1))
stopifnot(any(errors),!all(errors))
stopifnot(inherits(try(null_spatial_beta(matrix(0L,2,2)),silent=TRUE),'try-error'))
# Known inclusive tails, ties, add-one correction and capped two-sided comparison.
z<-null_inference(c(-1,0,0,1),0)
stopifnot(z$below==1,z$equal==2,z$above==1,z$p_lower==.8,z$p_upper==.8,z$p_two_sided==1)
z<-null_inference(c(0,1,2),-1)
stopifnot(z$p_lower==.25,z$p_upper==1,z$p_two_sided==.5,
  z$basic_low==2*z$observed_delta-z$null_q975)
stopifnot(isTRUE(all.equal(p.adjust(c(.01,.04,.2),'holm'),c(.03,.08,.2))))
cat('PASS: exact Jaccard, allocation frequencies, transition constraints, extreme counts, RNG, empty failures and tails.\n')
print(counts)
