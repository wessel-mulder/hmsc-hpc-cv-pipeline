source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc);library(corrplot);library(writexl)})
for (name in model_names()) {
  p <- paths(name); m <- read_fit(p)
  destination <- file.path(p$results,'parameter_estimates.rds');check_destination(destination)
  # Parameter summaries always use every retained draw, without further thinning.
  co <- convertToCodaObject(m)
  Beta <- getPostEstimate(m,parName='Beta'); Gamma <- getPostEstimate(m,parName='Gamma')
  Omega <- computeAssociations(m)
  VP <- computeVariancePartitioning(m)
  MF <- if (file.exists(p$fit)) read_objects(p$fit,'MF')$MF else {
    preds <- computePredictedValues(m,thin=1,expected=TRUE,nParallel=cfg$nparallel)
    fit <- evaluateModelFit(m,predY=preds);rm(preds);fit
  }
  summaries <- list(Beta=summarise_draws(as.matrix(co$Beta)),Gamma=summarise_draws(as.matrix(co$Gamma)))
  if (!is.null(co$Rho)) summaries$Rho <- summarise_draws(as.matrix(co$Rho))
  alpha_summaries <- list()
  for (r in seq_len(m$nr)) {
    # Each column is one latent factor; keep all factors and all chains.
    alpha <- summarise_draws(as.matrix(co$Alpha[[r]]))
    alpha$random_level <- names(m$ranLevels)[r]
    alpha_summaries[[r]] <- alpha
    summaries[[paste0('Omega_',r)]] <- summarise_draws(as.matrix(co$Omega[[r]]))
  }
  Alpha <- do.call(rbind,alpha_summaries)
  saveRDS(list(Beta=Beta,Gamma=Gamma,Omega=Omega,Alpha=Alpha,VP=VP,MF=MF,
               summaries=summaries,draws=sum(lengths(m$postList)),XScalePar=m$XScalePar),destination)
  for (block in names(summaries)) write.csv(summaries[[block]],file.path(p$results,paste0('parameter_',block,'.csv')),row.names=FALSE)
  write.csv(Alpha,file.path(p$results,'parameter_Alpha_all_factors.csv'),row.names=FALSE)
  saveRDS(VP,file.path(p$results,paste0(name,'VP.rds')))
  saveRDS(MF,file.path(p$results,paste0(name,'MF.rds')))
  write.csv(VP$vals,file.path(p$results,paste0(name,'parameter_estimates_VP_.csv')))
  for (block in c('Beta','Y')) if (!is.null(VP$R2T[[block]]))
    write.csv(VP$R2T[[block]],file.path(p$results,paste0(name,'parameter_estimates_VP_R2T_',block,'.csv')))
  export_estimate <- function(estimate,labels,label_name,file) {
    tables <- lapply(c('mean','support','supportNeg'),function(field) {
      x <- as.data.frame(t(estimate[[field]]));colnames(x) <- m$covNames
      out <- cbind(data.frame(label=labels),x);names(out)[1] <- label_name;out
    })
    names(tables) <- c('Posterior mean','Pr(x>0)','Pr(x<0)')
    writexl::write_xlsx(tables,file)
  }
  export_estimate(Beta,m$spNames,'Species',file.path(p$results,paste0(name,'parameter_estimates_Beta_.xlsx')))
  export_estimate(Gamma,m$trNames,'Traits',file.path(p$results,paste0(name,'parameter_estimates_Gamma_.xlsx')))
  for (r in seq_len(m$nr)) {
    tabs <- lapply(list(Omega[[r]]$mean,Omega[[r]]$support,1-Omega[[r]]$support),function(x) cbind(Species=m$spNames,as.data.frame(x)))
    names(tabs) <- c('Posterior mean','Pr(x>0)','Pr(x<0)')
    writexl::write_xlsx(tabs,file.path(p$results,paste0('parameter_estimates_Omega_',names(m$ranLevels)[r],'.xlsx')))
  }
  open_png(file.path(p$results,'parameter_estimates'),width=3000,height=2200)
  plotVariancePartitioning(m,VP,main=paste(name,'proportion of explained variance'))
  # Scale each species column by its own fit, without changing species order.
  R2 <- if (!is.null(MF$TjurR2)) MF$TjurR2 else MF$R2
  if (!is.null(R2) && all(is.finite(R2))) {
    VPr <- VP;VPr$vals <- sweep(VP$vals,2,R2,'*')
    plotVariancePartitioning(m,VPr,main=paste(name,'fit-weighted variance partition'))
  }
  for (param in c('Sign','Mean')) {
    plotBeta(m,post=Beta,supportLevel=.95,param=param,plotTree=TRUE,
             covNamesNumbers=c(TRUE,FALSE),spNamesNumbers=c(FALSE,FALSE),cex=c(.6,.6,.8))
    title(paste(name,'Beta',param))
    plotGamma(m,post=Gamma,supportLevel=.95,param=param,
              covNamesNumbers=c(TRUE,FALSE),trNamesNumbers=c(TRUE,FALSE),cex=c(.6,.6,.8))
    title(paste(name,'Gamma',param))
  }
  for (r in seq_len(m$nr)) {
    z <- Omega[[r]]
    signs <- (z$support>.95 | z$support<.05)*sign(z$mean)
    rownames(signs) <- colnames(signs) <- rep('',m$ns)
    corrplot(signs,method='color',type='lower',col=c('blue','white','red'),
             mar=c(0,0,3,0),main=paste(name,'residual associations',names(m$ranLevels)[r]))
  }
  dev.off()
  provenance(p$results,'S7_Show_parameter_estimates',p$fitted,
             c('All retained posterior draws; support threshold 0.95.',
               'Alpha summaries retain all latent factors; factor identities are not comparable across independent atlas models.'))
}
