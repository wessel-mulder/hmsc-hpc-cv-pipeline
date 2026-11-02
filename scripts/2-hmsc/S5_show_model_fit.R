source('scripts/2-hmsc/helpers/pipeline.R')
for (name in model_names()) {
  p <- paths(name); result <- read_objects(p$fit,c('MF','MFCV'))
  metrics <- intersect(names(result$MF),names(result$MFCV))
  destination <- file.path(p$results,'model_fit_by_species.csv');check_destination(destination)
  m <- read_objects(p$unfitted,'models')$models[[1]]
  output <- list()
  open_png(file.path(p$results,'model_fit'))
  for (metric in metrics) {
    x <- result$MF[[metric]]; y <- result$MFCV[[metric]]
    if (!is.numeric(x) || !is.numeric(y) || length(x)!=m$ns || length(y)!=m$ns) next
    output[[metric]] <- data.frame(species=colnames(m$Y),metric=metric,training=x,CV=y)
    keep <- is.finite(x)&is.finite(y)
    if (!any(keep)) {plot.new();title(paste(metric,': no finite estimates'));next}
    limits <- range(c(x[keep],y[keep]))
    if (diff(limits)==0) limits <- limits+c(-.05,.05)
    plot(x[keep],y[keep],xlim=limits,ylim=limits,pch=19,cex=.6,
         xlab='Training fit',ylab=paste(cfg$nfolds,'fold CV fit'),main=paste(name,metric,sep='\n'))
    abline(0,1,lty=2)
    if (metric %in% c('AUC','O.AUC')) {abline(h=.5,v=.5,col='grey')}
  }
  dev.off()
  write.csv(do.call(rbind,output),destination,row.names=FALSE)
  provenance(p$results,'S5_show_model_fit',c(p$fit,p$unfitted))
}
