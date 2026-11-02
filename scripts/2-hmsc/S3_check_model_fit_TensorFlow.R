source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc);library(coda)})
# Diagnose all retained draws; sampler transient was already discarded.
diagnostics <- function(x,block) {
  if (length(x)<2L) stop('At least two chains are needed for convergence checks')
  mat <- as.matrix(x)
  rhat <- tryCatch(gelman.diag(x,autoburnin=FALSE,multivariate=FALSE)$psrf[,1],
                   error=function(e) {warning(block,': ',conditionMessage(e));rep(NA_real_,ncol(mat))})
  ess <- tryCatch(effectiveSize(x),error=function(e) rep(NA_real_,ncol(mat)))
  data.frame(block=block,parameter=colnames(mat),Rhat=unname(rhat),ESS=unname(ess),
             mean=colMeans(mat),lower=apply(mat,2,quantile,.025),upper=apply(mat,2,quantile,.975),row.names=NULL)
}
for (name in model_names()) {
  p <- paths(name); m <- read_fit(p)
  destination <- file.path(p$results,'MCMC_diagnostics.csv');check_destination(destination)
  co <- convertToCodaObject(m)
  blocks <- list(Beta=co$Beta,Gamma=co$Gamma,Rho=co$Rho)
  for (r in seq_len(m$nr)) {
    blocks[[paste0('Alpha_',r)]] <- co$Alpha[[r]]
    blocks[[paste0('Omega_',r)]] <- co$Omega[[r]]
  }
  blocks <- blocks[!vapply(blocks,is.null,logical(1))]
  result <- do.call(rbind,lapply(names(blocks),function(b) diagnostics(blocks[[b]],b)))
  write.csv(result,destination,row.names=FALSE)
  summary <- do.call(rbind,lapply(split(result,result$block),function(d) {
    finite <- is.finite(d$Rhat)
    data.frame(block=d$block[1],parameters=nrow(d),undefined_Rhat=sum(!finite),
               Rhat_above_1.1=sum(d$Rhat>1.1,na.rm=TRUE),
               max_Rhat=if(any(finite)) max(d$Rhat[finite]) else NA_real_,
               min_ESS=if(any(is.finite(d$ESS))) min(d$ESS[is.finite(d$ESS)]) else NA_real_)
  }))
  write.csv(summary,file.path(p$results,'MCMC_diagnostic_summary.csv'),row.names=FALSE)
  open_png(file.path(p$results,'MCMC_convergence'))
  par(mfrow=c(2,1),mar=c(9,4,3,1))
  for (metric in c('Rhat','ESS')) {
    values <- lapply(split(result[[metric]],result$block),function(z) z[is.finite(z)])
    values <- values[lengths(values)>0]
    if (length(values)) {boxplot(values,las=2,ylab=metric,main=name,outline=FALSE);if(metric=='Rhat') abline(h=1.1,col='red',lty=2)}
    else {plot.new();title(paste('No finite',metric))}
  }
  dev.off()
  provenance(p$results,'S3_check_model_fit_TensorFlow',p$fitted,
             'All retained draws; autoburnin=FALSE. No automatic convergence approval.')
}
