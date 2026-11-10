# Pure functions shared by the runner and read-only report compatibility checks.
bootstrap_hashes <- function(root, out) {
  paths <- c('data-new/0-pre/raw/design/studyDesign.csv',
    'data-new/0-pre/derived/environmental/X_environment.csv',
    'data-new/0-pre/raw/occurrences/Y_occurrences.csv',
    'data-new/0-pre/raw/effort/dækning atlas1+2.xlsx')
  paths <- c(file.path(root,paths),file.path(out,c('homogenization-data.R',
    'bootstrap-functions.R','run-spatial-bootstrap.R')))
  setNames(unname(tools::md5sum(paths)),sub(paste0(root,'/'),'',paths,fixed=TRUE))
}
jaccard_matrix <- function(x) {
  a <- tcrossprod(x); u <- outer(rowSums(x),rowSums(x),'+')-a
  stopifnot(all(u>0)); d <- 1-a/u; diag(d)<-0; d
}
block_sums <- function(d,g) {
  # Rowsum twice gives ordered block-to-block sums, including within-block pairs.
  rowsum(t(rowsum(d,g,reorder=FALSE)),g,reorder=FALSE)
}
weighted_beta <- function(A,M,sizes) {
  n <- drop(crossprod(sizes,M)); stopifnot(all(n>=2))
  colSums(M*(A %*% M))/(n*(n-1))
}
bootstrap_summary <- function(draws,observed) {
  do.call(rbind,lapply(split(draws,interaction(draws$config,draws$batch,drop=TRUE)),function(z) {
    do.call(rbind,lapply(names(observed),function(pair) {
      x <- z[[paste0('delta_',pair)]]; q <- quantile(x,c(.025,.975),names=FALSE)
      data.frame(config=z$config[1],width_km=z$width_km[1],origin=z$origin[1],
        batch=z$batch[1],pair=pair,observed=unname(observed[pair]),
        bootstrap_mean=mean(x),displacement=mean(x)-observed[pair],se=sd(x),
        percentile_low=q[1],percentile_high=q[2],
        basic_low=2*observed[pair]-q[2],basic_high=2*observed[pair]-q[1],row.names=NULL)
    }))
  }))
}
