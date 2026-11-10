null_hashes <- function(root,out) {
  paths<-c(file.path(root,c('data-new/0-pre/raw/design/studyDesign.csv',
    'data-new/0-pre/derived/environmental/X_environment.csv',
    'data-new/0-pre/raw/occurrences/Y_occurrences.csv',
    'data-new/0-pre/raw/effort/dækning atlas1+2.xlsx')),
    file.path(out,c('homogenization-data.R','null-functions.R','homogenization-null.cpp','run-homogenization-null.R')))
  setNames(unname(tools::md5sum(paths)),sub(paste0(root,'/'),'',paths,fixed=TRUE))
}
null_inference <- function(delta,observed) {
  B<-length(delta); stopifnot(B>0,all(is.finite(delta)),is.finite(observed))
  below<-sum(delta<observed); equal<-sum(delta==observed); above<-sum(delta>observed)
  lower<-(1+below+equal)/(B+1); upper<-(1+above+equal)/(B+1)
  bin_ci<-function(k) unname(binom.test(k,B)$conf.int)
  lc<-bin_ci(below+equal); uc<-bin_ci(above+equal)
  data.frame(simulations=B,observed_delta=observed,null_mean=mean(delta),null_sd=sd(delta),
    null_q025=unname(quantile(delta,.025)),null_q975=unname(quantile(delta,.975)),
    observed_minus_null=observed-mean(delta),below=below,equal=equal,above=above,
    p_lower=lower,p_upper=upper,p_two_sided=min(1,2*min(lower,upper)),
    lower_tail_mc_low=lc[1],lower_tail_mc_high=lc[2],
    upper_tail_mc_low=uc[1],upper_tail_mc_high=uc[2],minimum_one_sided_p=1/(B+1))
}
atomic_rds <- function(value,path) {
  tmp<-paste0(path,'.tmp'); saveRDS(value,tmp,compress=FALSE)
  if(!file.rename(tmp,path)) stop('Failed atomic checkpoint write: ',path)
}
