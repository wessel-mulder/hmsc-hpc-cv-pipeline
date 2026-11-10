# Run from repository root or this directory. Each run gets its own result folder.
root <- normalizePath(if (dir.exists('data-new')) '.' else '../../..')
analysis_dir <- file.path(root,'scripts/1-atlas-analyses/unusual-temporal-turnover')
results_root <- file.path(analysis_dir,'simulation-results')
null_model <- 'independent-change'
args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name,default) {
  z <- grep(paste0('^--',name,'='),args,value=TRUE)
  if (!length(z)) default else sub(paste0('^--',name,'='),'',z[1])
}
n_sim <- as.integer(arg_value('simulations','999999'))
chunk_size <- as.integer(arg_value('chunk-size','10000'))
seed <- as.integer(arg_value('seed','20261104'))
slurm_workers <- suppressWarnings(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','')))
workers <- as.integer(arg_value('workers',if (is.finite(slurm_workers)&&slurm_workers>0) slurm_workers else 6L))
stopifnot(length(n_sim)==1L,!is.na(n_sim),n_sim>0L,
          length(chunk_size)==1L,!is.na(chunk_size),chunk_size>0L,
          length(seed)==1L,!is.na(seed),seed>0L,
          length(workers)==1L,!is.na(workers),workers>0L)

input_paths <- file.path(root,c('data-new/0-pre/raw/design/studyDesign.csv',
  'data-new/0-pre/derived/environmental/X_environment.csv',
  'data-new/0-pre/raw/occurrences/Y_occurrences.csv',
  'data-new/0-pre/raw/effort/dækning atlas1+2.xlsx'))
stopifnot(all(file.exists(input_paths)))
input_mtime <- as.numeric(file.info(input_paths)$mtime)
script_path <- file.path(analysis_dir,'unusual-temporal-turnover.R')
script_mtime <- as.numeric(file.info(script_path)$mtime)
timestamp <- format(Sys.time(),'%Y%m%dT%H%M%S')
run_stub <- sprintf('%s_nsim-%d_seed-%d',timestamp,n_sim,seed)
model_dir <- file.path(results_root,null_model)
dir.create(model_dir,recursive=TRUE,showWarnings=FALSE)
run_id <- run_stub; suffix <- 1L
while (dir.exists(file.path(model_dir,run_id))) {
  suffix <- suffix+1L; run_id <- paste0(run_stub,'-',suffix)
}
run_dir <- file.path(model_dir,run_id)
dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
writeLines('running',file.path(run_dir,'run-status.txt'))
progress_log <- file.path(run_dir,'simulation-progress.log')
progress_lock <- file.path(run_dir,'.progress-lock')
progress_start <- Sys.time()
writeLines(c('timestamp\tpair\tcompleted\ttotal\telapsed_seconds\trate_per_second\teta_seconds',
  vapply(c('1-2','2-3','1-3'),function(p) paste(format(progress_start,'%Y-%m-%dT%H:%M:%S'),p,0,n_sim,0,0,NA,sep='\t'),character(1))),progress_log)
log_progress <- function(pair,increment) {
  while (!dir.create(progress_lock,showWarnings=FALSE)) Sys.sleep(0.05)
  on.exit(unlink(progress_lock,recursive=TRUE),add=TRUE)
  lines <- readLines(progress_log,warn=FALSE)
  pair_lines <- lines[grepl(paste0('\t',pair,'\t'),lines,fixed=TRUE)]
  last <- strsplit(tail(pair_lines,1),'\t',fixed=TRUE)[[1]]
  completed <- as.integer(last[3])+increment
  elapsed <- as.numeric(difftime(Sys.time(),progress_start,units='secs'))
  rate <- completed/max(elapsed,1)
  eta <- (n_sim-completed)/rate
  line <- paste(format(Sys.time(),'%Y-%m-%dT%H:%M:%S'),pair,completed,n_sim,
                round(elapsed,1),signif(rate,5),round(eta),sep='\t')
  cat(line,'\n',file=progress_log,append=TRUE,sep='')
  flush.console()
}

design <- read.csv(input_paths[1],check.names=FALSE)
env <- read.csv(input_paths[2],check.names=FALSE)
y <- read.csv(input_paths[3],row.names=1,check.names=FALSE)
stopifnot(all(c('site','atlas','survey','lat','lon') %in% names(design)),
          all(c('survey','kvadratkod','atlas','pct_class_0','pct_class_77') %in% names(env)),
          !anyDuplicated(design$survey),!anyDuplicated(env[c('kvadratkod','atlas')]),
          setequal(design$survey,env$survey),
          all(env$survey==paste0(env$kvadratkod,'_',env$atlas)),
          !anyDuplicated(rownames(y)),!anyDuplicated(names(y)),ncol(y)==201L)
env$water_percent <- env$pct_class_0+env$pct_class_77
stopifnot(!anyNA(env$water_percent),all(env$water_percent>=0 & env$water_percent<=100+1e-6))
by_atlas <- lapply(1:3,function(a) intersect(
  as.character(env$kvadratkod[env$atlas==a & env$water_percent<=75]),
  as.character(design$site[design$atlas==a])))
sites <- sort(Reduce(intersect,by_atlas))
e1 <- readxl::read_excel(input_paths[4],sheet='kvadrater atlas1')
e2 <- readxl::read_excel(input_paths[4],sheet='kvadrater atlas2')
cov1 <- as.character(e1$KVADRATKOD[e1[['UNDS_KODE_']] %in% c(1,2)])
cov2 <- as.character(e2$KVADRATKOD[e2[['UNDS_KOD_1']] %in% c(1,2)])
sites <- sort(Reduce(intersect,list(sites,cov1,cov2)))
stopifnot(length(sites)>0L)
d <- design[design$site %in% sites,,drop=FALSE]
stopifnot(all(vapply(1:3,function(a) all(sites %in% d$site[d$atlas==a]),logical(1))))
matrices <- lapply(1:3,function(a) {
  ids <- paste0(sites,'_',a); stopifnot(all(ids %in% rownames(y)))
  z <- as.matrix(y[ids,,drop=FALSE]); rownames(z)<-sites; storage.mode(z)<-'integer'; z
})
stopifnot(all(vapply(matrices,function(z) !anyNA(z)&&all(z %in% 0:1),logical(1))))
coords <- unique(d[d$atlas==1,c('site','lat','lon')]); coords <- coords[match(sites,coords$site),]
stopifnot(identical(as.character(coords$site),sites),!anyNA(coords[c('lat','lon')]))
pairs <- list(c(1,2),c(2,3),c(1,3)); names(pairs)<-c('1-2','2-3','1-3')
observed <- lapply(names(pairs),function(pair_name) {
  ab <- pairs[[pair_name]]; x <- matrices[[ab[1]]]; z <- matrices[[ab[2]]]
  union_n <- rowSums((x+z)>0); changed_n <- rowSums(x!=z)
  tbi <- changed_n/pmax(1,union_n)
  list(pair=pair_name,cell_tbi=tbi,community_mean=mean(tbi))
})
names(observed)<-names(pairs)

Rcpp::cppFunction(code='
  #include <random>
  #include <unordered_set>
  #include <vector>
  Rcpp::List independent_null(Rcpp::LogicalMatrix x,Rcpp::LogicalMatrix y,
                              Rcpp::NumericVector obs,int nsim,unsigned int seed) {
    int n=x.nrow(),s=x.ncol(); std::mt19937 rng(seed);
    std::vector<std::vector<int> > occ(s),unocc(s);
    std::vector<int> losses(s,0),gains(s,0),richness(n,0);
    for(int j=0;j<s;++j) for(int i=0;i<n;++i) {
      if(x(i,j)){occ[j].push_back(i);richness[i]++;} else unocc[j].push_back(i);
      if(x(i,j)&&!y(i,j))losses[j]++;
      if(!x(i,j)&&y(i,j))gains[j]++;
    }
    Rcpp::NumericVector sum(n),sumsq(n),global_mean(nsim);
    Rcpp::IntegerVector below(n),equal(n),above(n);
    for(int b=0;b<nsim;++b) {
      std::vector<int> changed(n,0),united(richness);
      for(int j=0;j<s;++j) for(int kind=0;kind<2;++kind) {
        const std::vector<int>& pool=kind==0?occ[j]:unocc[j];
        int k=kind==0?losses[j]:gains[j],m=pool.size();
        std::unordered_set<int> chosen; chosen.reserve(k*2+1);
        for(int q=m-k;q<m;++q) {
          std::uniform_int_distribution<int> draw(0,q);
          int t=draw(rng),pick=chosen.insert(t).second?t:q;
          if(pick==q)chosen.insert(q);
          int cell=pool[pick];changed[cell]++;
          if(kind==1)united[cell]++;
        }
        if((int)chosen.size()!=k)Rcpp::stop("A species transition was not sampled uniquely.");
      }
      double total=0.0;
      for(int i=0;i<n;++i) {
        double v=united[i]?(double)changed[i]/united[i]:0.0;
        sum[i]+=v;sumsq[i]+=v*v;total+=v;
        if(v<obs[i])below[i]++;else if(v>obs[i])above[i]++;else equal[i]++;
      }
      global_mean[b]=total/n;
    }
    return Rcpp::List::create(Rcpp::Named("sum")=sum,Rcpp::Named("sumsq")=sumsq,
      Rcpp::Named("below")=below,Rcpp::Named("equal")=equal,
      Rcpp::Named("above")=above,Rcpp::Named("global_mean")=global_mean);
  }')

tasks <- list(); task_i <- 0L; n_chunks <- ceiling(n_sim/chunk_size)
for (chunk_i in seq_len(n_chunks)) {
  draws_this_chunk <- min(chunk_size,n_sim-(chunk_i-1L)*chunk_size)
  for (pair_i in seq_along(pairs)) {
    task_i <- task_i+1L
    tasks[[task_i]] <- list(pair_index=pair_i,pair=names(pairs)[pair_i],
      chunk=chunk_i,nsim=draws_this_chunk,seed=as.integer(seed+pair_i*100000L+chunk_i))
  }
}
task_fun <- function(task) {
  ab <- pairs[[task$pair_index]]; pair <- task$pair
  value <- independent_null(matrices[[ab[1]]]==1L,matrices[[ab[2]]]==1L,
    observed[[pair]]$cell_tbi,task$nsim,task$seed)
  log_progress(pair,task$nsim)
  value
}
message('Running ',length(tasks),' chunks across ',workers,' workers; ',
        n_sim,' simulations per atlas pair. Progress: ',progress_log)
chunk_results <- if(workers>1L&&.Platform$OS.type!='windows') {
  parallel::mclapply(tasks,task_fun,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
} else {
  lapply(tasks,task_fun)
}
failed <- vapply(chunk_results,function(z)inherits(z,'try-error')||is.null(z),logical(1))
if(any(failed))stop('Simulation chunk failed: ',which(failed)[1],'; see ',progress_log)

result <- vector('list',length(pairs)); names(result)<-names(pairs)
for (pair_name in names(pairs)) {
  ix <- which(vapply(tasks,function(z)z$pair==pair_name,logical(1)))
  chunks <- chunk_results[ix]
  draws <- sum(vapply(tasks[ix],function(z)z$nsim,integer(1)))
  sum_site <- Reduce(`+`,lapply(chunks,`[[`,'sum'))
  sumsq_site <- Reduce(`+`,lapply(chunks,`[[`,'sumsq'))
  null_mean <- sum_site/draws
  null_variance <- pmax(0,(sumsq_site-sum_site^2/draws)/(draws-1))
  below <- Reduce(`+`,lapply(chunks,`[[`,'below'))
  equal <- Reduce(`+`,lapply(chunks,`[[`,'equal'))
  above <- Reduce(`+`,lapply(chunks,`[[`,'above'))
  global_mean_draws <- unlist(lapply(chunks,`[[`,'global_mean'),use.names=FALSE)
  obs <- observed[[pair_name]]$cell_tbi
  obs_global <- observed[[pair_name]]$community_mean
  stopifnot(draws==n_sim,length(global_mean_draws)==n_sim,
    all(below+equal+above==n_sim),all(is.finite(null_mean)),
    all(is.finite(null_variance)),all(is.finite(global_mean_draws)))
  result[[pair_name]] <- list(pair=pair_name,observed=obs,null_mean=null_mean,
    null_variance=null_variance,below=below,equal=equal,above=above,
    global_mean_draws=global_mean_draws,global_null_mean=mean(global_mean_draws),
    global_null_variance=var(global_mean_draws),observed_global_mean=obs_global,
    global_upper_p=(1+sum(global_mean_draws>=obs_global))/(n_sim+1),n_sim=n_sim)
  message(pair_name,': ',draws,' draws; observed grid mean=',round(obs_global,5),
    '; null grid mean=',round(result[[pair_name]]$global_null_mean,5),
    '; per-site tail counts validated.')
}
site_summary <- do.call(rbind,lapply(result,function(z)data.frame(
  pair=z$pair,site=sites,observed_tbi=z$observed,null_mean=z$null_mean,
  null_variance=z$null_variance,below=z$below,equal=z$equal,above=z$above)))
metadata <- list(null_model=null_model,run_id=run_id,run_timestamp=timestamp,
  simulations_per_pair=n_sim,chunk_size=chunk_size,workers=workers,base_seed=seed,
  input_paths=input_paths,input_mtime=input_mtime,script_path=script_path,
  script_mtime=script_mtime,n_cells=length(sites),n_species=ncol(y))
save(root,analysis_dir,results_root,run_dir,run_id,null_model,sites,coords,matrices,
  pairs,observed,result,site_summary,metadata,n_sim,chunk_size,workers,seed,
  file=file.path(run_dir,'simulation-results.RData'),compress='xz')

writeLines('complete',file.path(run_dir,'run-status.txt'))
writeLines(file.path(null_model,run_id),file.path(results_root,'latest-run.txt'))
source(file.path(analysis_dir,'analyse-site-turnover.R'))
analyse_site_turnover_run(run_dir)
message('Completed run ',run_id,'; results: ',run_dir)
