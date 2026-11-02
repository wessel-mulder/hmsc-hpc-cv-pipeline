# Environmental response curves only. Cross-atlas tests are retired.
source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages(library(Hmsc))
for (name in model_names()) {
  p <- paths(name); m <- read_fit(p)
  directory <- file.path(p$results,'Preds');dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  post <- posterior_draws(m)
  message('Using ',length(post),' retained draws for ',name)
  # Never standardize new data manually: predict() applies the fitted transformation.
  for (variable in all.vars(m$XFormula)) {
    destination <- file.path(directory,paste0('response_',variable,'.rds'))
    check_destination(destination)
    curves <- list()
    for (effect in c('total','marginal')) {
      gradient <- if (effect=='total') constructGradient(m,focalVariable=variable,ngrid=cfg$ngrid) else
        constructGradient(m,focalVariable=variable,non.focalVariables=1,ngrid=cfg$ngrid)
      predictions <- predict(m,Gradient=gradient,post=post,expected=TRUE,
                             nParallel=cfg$nparallel,useSocket=FALSE,mcmcStep=1)
      prediction_array(predictions,cfg$ngrid,m$ns,length(post))
      curves[[effect]] <- list(Gradient=gradient,predY=predictions)
    }
    saveRDS(list(variable=variable,draws=length(post),posterior_thin=cfg$posterior_thin,
                 XScalePar=m$XScalePar,curves=curves),destination)
    open_png(file.path(directory,paste0('response_',variable)))
    for (effect in names(curves)) {
      z <- curves[[effect]]
      pl <- plotGradient(m,z$Gradient,pred=z$predY,yshow=0,measure='S',showData=TRUE,
                          main=paste(name,variable,effect,'summed occurrence probability',sep='\n'))
      if (inherits(pl,'ggplot')) print(pl)
      for (trait in seq.int(2,m$nt)) {
        pl <- plotGradient(m,z$Gradient,pred=z$predY,index=trait,measure='T',yshow=0,showData=TRUE,
                            main=paste(name,variable,effect,m$trNames[trait],sep='\n'))
        if (inherits(pl,'ggplot')) print(pl)
      }
    }
    dev.off()
    rm(curves);gc()
  }
  provenance(directory,'S6_Making_Predictions',p$fitted,paste('Retained draws:',length(post)))
}
