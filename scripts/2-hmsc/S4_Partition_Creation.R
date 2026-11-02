# Ten random folds; four independently seeded initialization files per fold.
source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc);library(jsonify)})
for (name in model_names()) {
  p <- paths(name); hM <- read_fit(p)
  input_md5 <- unname(tools::md5sum(p$fitted))
  threads <- cfg$nfolds*cfg$chains
  destinations <- file.path(p$temp,sprintf('temp_samples_%04d_thin_%02d_thread_%d.rds',cfg$samples,cfg$thin,seq_len(threads)))
  for (destination in destinations) check_destination(destination)
  if (file.exists(p$folds)) {
    previous <- read_objects(p$folds,c('partition','seeds','input_md5','row_ids','idfold'))
    stopifnot(identical(previous$input_md5,input_md5),identical(previous$row_ids,rownames(hM$Y)),
              length(previous$seeds)==threads,identical(sort(unique(previous$partition)),seq_len(cfg$nfolds)))
    partition <- previous$partition; seeds <- previous$seeds
    message('Reusing saved fold assignments and seeds')
  } else {
    partition <- createPartition(hM,nfolds=cfg$nfolds)
    seeds <- sample.int(.Machine$integer.max,threads)
  }
  idfold <- rep(seq_len(cfg$nfolds),each=cfg$chains)
  row_ids <- rownames(hM$Y)
  hM1 <- lapply(seq_len(cfg$nfolds),function(k) make_fold_model(hM,partition!=k))
  save(partition,seeds,input_md5,row_ids,idfold,threads,hM1,hM,file=p$folds)
  write.csv(data.frame(survey=row_ids,fold=partition),file.path(p$temp,'CV_fold_assignments.csv'),row.names=FALSE)
  write.csv(data.frame(thread=seq_len(threads),fold=idfold,seed=seeds),file.path(p$temp,'CV_thread_seeds.csv'),row.names=FALSE)
  for (i in seq_len(threads)) {
    set.seed(seeds[i])
    init <- sampleMcmc(hM1[[idfold[i]]],samples=cfg$samples,thin=cfg$thin,
                       transient=cfg$transient,adaptNf=hM$adaptNf,nChains=1,nParallel=1,
                       updater=hM$updater,verbose=1,engine='HPC')
    saveRDS(to_json(init),destinations[i])
    message('Exported CV thread ',i,'/',threads,' for ',name)
  }
  provenance(p$temp,'S4_Partition_Creation',p$fitted,
             'Random-fold interpolation; full-run predictor scaling retained as in Hmsc CV.')
}
