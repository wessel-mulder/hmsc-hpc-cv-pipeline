# Validate saved model definitions without fitting or generating predictions.
# Rscript scripts/2-hmsc/tests/check-S1-model-definitions.R /path/to/S1/output [old-output-root]
suppressPackageStartupMessages(library(Hmsc))
args<-commandArgs(trailingOnly=TRUE)
if(!length(args))stop('Supply the output root of a completed S1 run')
raw<-new.env();load('data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData',envir=raw)
files<-list.files(args[1],pattern='unfitted_models[.]RData$',recursive=TRUE,full.names=TRUE)
stopifnot(length(files)==3)
for(a in 1:3){
  f<-files[grepl(paste0('_Atlas',a,'_'),files)];stopifnot(length(f)==1)
  e<-new.env();load(f,envir=e);stopifnot(length(e$models)==1)
  m<-e$models[[1]]
  idx<-as.character(raw$Design$atlas)==as.character(a)
  expected_x<-raw$X[idx,,drop=FALSE];expected_y<-raw$Y[idx,,drop=FALSE]
  stopifnot(identical(m$XData,expected_x),identical(m$Y,as.matrix(expected_y)),
            identical(rownames(m$XData),rownames(m$Y)),ncol(m$XData)==10,
            ncol(m$Y)==155,length(m$phyloTree$tip.label)==155,
            setequal(m$phyloTree$tip.label,colnames(m$Y)),
            is.factor(m$TrData$Migration_a3_DOF),is.factor(m$TrData$foraging_guild_consensus),
            m$rL[[1]]$nfMin==5,m$rL[[1]]$nfMax==5)
  stopifnot(max(abs(m$XScaled[,-1,drop=FALSE]-scale(as.matrix(expected_x))))<1e-10,
            max(abs(m$XScalePar[1,-1]-colMeans(expected_x)))<1e-10,
            max(abs(m$XScalePar[2,-1]-vapply(expected_x,sd,numeric(1))))<1e-10,
            all(m$XScaled[,1]==1),qr(m$XScaled)$rank==ncol(m$XScaled))
  if(length(args)>1){
    oldfiles<-list.files(args[2],pattern='unfitted_models[.]RData$',recursive=TRUE,full.names=TRUE)
    oldfile<-oldfiles[grepl(paste0('_Atlas',a,'_'),oldfiles)];stopifnot(length(oldfile)==1)
    old<-new.env();load(oldfile,envir=old);before<-old$models[[1]]
    # Constructor call text changes with local variable names; compare spatial settings.
    current_levels<-m$rL; previous_levels<-before$rL
    for (level in names(current_levels)) {
      current_levels[[level]]$call<-NULL
      previous_levels[[level]]$call<-NULL
    }
    stopifnot(identical(m$Y,before$Y),identical(m$TrData,before$TrData),
              identical(m$studyDesign,before$studyDesign),
              isTRUE(all.equal(m$phyloTree,before$phyloTree)),
              max(abs(m$XScaled-before$XScaled))<1e-10,
              isTRUE(all.equal(current_levels,previous_levels)))
  }
  stopifnot(!file.exists(file.path(dirname(dirname(dirname(f))),'Tests','test_atlases.RData')))
  message('PASS Atlas ',a,': raw predictors, scaling, traits, tree, spatial priors and save interface')
}
