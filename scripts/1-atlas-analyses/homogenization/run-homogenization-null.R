# Explicit simulation runner; report rendering reads saved results only.
root<-normalizePath(if(dir.exists('data-new')) '.' else '../../..')
out<-file.path(root,'scripts/1-atlas-analyses/homogenization')
source(file.path(out,'homogenization-data.R'))
source(file.path(out,'null-functions.R'))
args<-commandArgs(trailingOnly=TRUE)
arg<-function(key,default) {
  x<-grep(paste0('^--',key,'='),args,value=TRUE)
  if(length(x)) sub(paste0('^--',key,'='),'',x[1]) else default
}
resume<-arg('resume','')
slurm<-suppressWarnings(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','')))
workers<-as.integer(arg('workers',if(is.finite(slurm)&&slurm>0) slurm else 6L))
stopifnot(!is.na(workers),workers>0)
results_root<-file.path(out,'null-results/independent-change')
if(nzchar(resume)) {
  run_dir<-normalizePath(resume,mustWork=TRUE)
  saved<-readRDS(file.path(run_dir,'settings.rds'))
  defaults<-saved$settings
} else defaults<-list(simulations=9999L,seed=20261105L,chunk_size=250L)
settings<-list(simulations=as.integer(arg('simulations',defaults$simulations)),
  seed=as.integer(arg('seed',defaults$seed)),chunk_size=as.integer(arg('chunk-size',defaults$chunk_size)))
stopifnot(all(vapply(settings,function(x) length(x)==1 && !is.na(x) && x>0,logical(1))),settings$simulations>=2)
hashes<-null_hashes(root,out)
Rcpp::sourceCpp(file.path(out,'homogenization-null.cpp'))
matrices<-lapply(matrices,function(x) {storage.mode(x)<-'integer';x})
beta<-vapply(matrices,null_spatial_beta,numeric(1))
# Independent direct check of observed statistics.
for(a in 1:3) {
  x<-matrices[[a]]; sh<-tcrossprod(x); u<-outer(rowSums(x),rowSums(x),'+')-sh
  d<-1-sh/u
  stopifnot(abs(beta[a]-mean(d[lower.tri(d)]))<1e-11)
}
pairs<-list('1-2'=c(1,2),'2-3'=c(2,3),'1-3'=c(1,3))
observed<-vapply(pairs,function(ab) beta[ab[2]]-beta[ab[1]],numeric(1))
constraints<-lapply(pairs,function(ab) {
  x<-matrices[[ab[1]]]; y<-matrices[[ab[2]]]
  data.frame(species=colnames(x),earlier_occupancy=colSums(x),later_occupancy=colSums(y),
    losses=colSums(x==1 & y==0),gains=colSums(x==0 & y==1),row.names=NULL)
})
# Audit filtering without changing the shared data helper (and bootstrap hashes).
all_sites<-sort(unique(as.character(design$site)))
land_common<-sort(Reduce(intersect,eligible_sites))
filter_audit<-data.frame(stage=c('Any atlas design cell','Present in all three atlases',
  'Land >=25% in all three atlases','Also Good/Average Atlas 1','Also Good/Average Atlas 2'),
  cells=c(length(all_sites),length(Reduce(intersect,lapply(1:3,function(a) as.character(design$site[design$atlas==a])))),
    length(land_common),length(intersect(land_common,good1)),length(sites)))
manifest<-design[match(paste0(sites,'_1'),design$survey),c('site','lon','lat')]
stopifnot(identical(as.character(manifest$site),sites))
# Check exact observed-analysis sample using its shared helper in an isolated environment.
reference<-new.env(); sys.source(file.path(out,'homogenization-data.R'),reference)
stopifnot(identical(reference$sites,sites),all(vapply(1:3,function(a)
  isTRUE(all.equal(matrices[[a]],reference$matrices[[a]],check.attributes=TRUE)),logical(1))))
B<-settings$simulations; n_chunks<-ceiling(B/settings$chunk_size)
stopifnot(as.double(settings$seed)+3*1000000+n_chunks<.Machine$integer.max)
tasks<-do.call(c,lapply(seq_along(pairs),function(p) lapply(seq_len(n_chunks),function(ch) {
  first<-(ch-1L)*settings$chunk_size+1L
  list(pair=names(pairs)[p],chunk=ch,first=first,draws=min(settings$chunk_size,B-first+1L),
    seed=as.integer(settings$seed+p*1000000L+ch))
})))
if(nzchar(resume)) {
  if(!identical(saved$settings,settings)||!identical(saved$hashes,hashes)||!identical(saved$tasks,tasks))
    stop('Incompatible resume: settings, input/code hashes or seed schedule changed.')
  if(file.exists(file.path(run_dir,'run-status.txt')) &&
     identical(readLines(file.path(run_dir,'run-status.txt')),'complete')) stop('Run already complete; use its saved results.')
} else {
  dir.create(results_root,recursive=TRUE,showWarnings=FALSE)
  run_dir<-tempfile(paste0(format(Sys.time(),'%Y%m%dT%H%M%S'),'_nsim-',B,'_seed-',settings$seed,'_'),tmpdir=results_root)
  dir.create(run_dir)
  atomic_rds(list(settings=settings,hashes=hashes,tasks=tasks),file.path(run_dir,'settings.rds'))
}
checkpoint_dir<-file.path(run_dir,'checkpoints');dir.create(checkpoint_dir,showWarnings=FALSE)
writeLines('running',file.path(run_dir,'run-status.txt'))
log_path<-file.path(run_dir,'simulation-progress.log'); started<-Sys.time()
log<-function(line) { message(line); cat(format(Sys.time()),line,'\n',file=log_path,append=TRUE) }
checkpoint_path<-function(task) file.path(checkpoint_dir,sprintf('%s_%05d.rds',task$pair,task$chunk))
valid_checkpoint<-function(task) {
  path<-checkpoint_path(task); if(!file.exists(path)) return(FALSE)
  z<-tryCatch(readRDS(path),error=function(e) NULL)
  if(is.null(z)||!identical(z$task,task)||length(z$beta_later)!=task$draws||
     any(!is.finite(z$beta_later))||any(z$beta_later<0|z$beta_later>1))
    stop('Invalid checkpoint: ',path)
  TRUE
}
complete<-vapply(tasks,valid_checkpoint,logical(1)); pending<-which(!complete)
log(paste('Run:',run_dir,'; validated existing chunks:',sum(complete),'/',length(tasks)))
run_task<-function(task) {
  tryCatch({
    ab<-pairs[[task$pair]]
    z<-simulate_spatial_null(matrices[[ab[1]]],matrices[[ab[2]]],task$draws,task$seed)
    list(ok=TRUE,task=task,beta_later=z$beta_later)
  },error=function(e) list(ok=FALSE,task=task,error=conditionMessage(e)))
}
tryCatch({
  batches<-split(pending,ceiling(seq_along(pending)/workers))
  for(ids in batches) {
    result<-if(workers>1 && .Platform$OS.type!='windows')
      parallel::mclapply(tasks[ids],run_task,mc.cores=min(workers,length(ids)),mc.set.seed=FALSE) else lapply(tasks[ids],run_task)
    failed<-list()
    for(k in seq_along(ids)) {
      z<-result[[k]]; task<-tasks[[ids[k]]]
      if(!is.list(z)||!isTRUE(z$ok)) { failed[[length(failed)+1L]]<-z; next }
      atomic_rds(list(task=task,beta_later=z$beta_later),checkpoint_path(task));complete[ids[k]]<-TRUE
    }
    if(length(failed)) {
      saveRDS(list(failures=failed,manifest=manifest,constraints=constraints),file.path(run_dir,'failure-details.rds'))
      stop('Simulation chunk failed; see failure-details.rds. ',paste(vapply(failed,function(z)
        if(is.list(z)) z$error else as.character(z),character(1)),collapse='; '))
    }
    elapsed<-as.numeric(difftime(Sys.time(),started,units='secs'))
    done_now<-sum(complete[pending]); completed_draws<-vapply(names(pairs),function(p)
      sum(vapply(tasks[complete],function(t) if(t$pair==p)t$draws else 0L,integer(1))),integer(1))
    eta<-if(done_now) elapsed*(length(pending)-done_now)/done_now else NA_real_
    log(paste(paste(names(completed_draws),completed_draws,sep='=',collapse=' '),
      '; elapsed',round(elapsed,1),'s; remaining estimate',round(eta,1),'s'))
  }
  draws<-do.call(rbind,lapply(tasks,function(task) {
    z<-readRDS(checkpoint_path(task)); base<-beta[pairs[[task$pair]][1]]
    data.frame(pair=task$pair,draw=seq.int(task$first,length.out=task$draws),
      beta_later=z$beta_later,delta=z$beta_later-base)
  }))
  stopifnot(nrow(draws)==3*B,all(table(draws$pair)==B),!anyDuplicated(draws[c('pair','draw')]))
  inference<-do.call(rbind,lapply(names(pairs),function(p) cbind(pair=p,
    null_inference(draws$delta[draws$pair==p],observed[p]))))
  inference$p_holm<-p.adjust(inference$p_two_sided,'holm')
  metadata<-list(settings=settings,input_code_hashes=hashes,seed_schedule=tasks,
    workers=workers,completed=Sys.time(),session=sessionInfo(),compiler=system2(file.path(R.home('bin'),'R'),c('CMD','config','CXX'),stdout=TRUE),
    compiler_version=tryCatch(system2(strsplit(system2(file.path(R.home('bin'),'R'),c('CMD','config','CXX'),stdout=TRUE),' ')[[1]][1],'--version',stdout=TRUE),error=function(e) conditionMessage(e)),
    model='Baseline-fixed independent allocation of observed per-species gains/losses',
    empty_community_policy='Fail chunk; retain checkpoints; never discard or redraw',
    validation='Binary states, unique transition selections, endpoint occupancy, richness total and observed R agreement checked')
  final_path<-file.path(run_dir,'null-results.RData'); tmp<-paste0(final_path,'.tmp')
  save(draws,inference,observed,beta,constraints,manifest,filter_audit,sites,metadata,file=tmp,compress='xz')
  check<-new.env();load(tmp,envir=check)
  stopifnot(identical(check$draws,draws),identical(check$metadata$input_code_hashes,hashes))
  stopifnot(file.rename(tmp,final_path))
  # Delete only this run's redundant chunk checkpoints after read-back validation.
  unlink(vapply(tasks,checkpoint_path,character(1)))
  if(length(list.files(checkpoint_dir,all.files=TRUE,no..=TRUE))) stop('Unexpected files remain in checkpoint folder')
  unlink(checkpoint_dir,recursive=TRUE)
  writeLines('complete',file.path(run_dir,'run-status.txt'))
  pointer_tmp<-file.path(results_root,'latest-run.txt.tmp');writeLines(basename(run_dir),pointer_tmp)
  stopifnot(file.rename(pointer_tmp,file.path(results_root,'latest-run.txt')))
  log('Complete. Final results written and verified; redundant checkpoints removed.')
  print(inference,row.names=FALSE)
},error=function(e) {
  writeLines('failed',file.path(run_dir,'run-status.txt'));log(paste('FAILED:',conditionMessage(e)));stop(e)
})
