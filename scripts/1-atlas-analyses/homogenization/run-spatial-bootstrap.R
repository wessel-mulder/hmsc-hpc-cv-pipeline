# Explicit runner: report rendering never launches this simulation.
root <- normalizePath(if(dir.exists('data-new')) '.' else '../../..')
out <- file.path(root,'scripts/1-atlas-analyses/homogenization')
source(file.path(out,'homogenization-data.R'))
source(file.path(out,'bootstrap-functions.R'))
args <- commandArgs(trailingOnly=TRUE)
arg <- function(key,default) {
  x <- grep(paste0('^--',key,'='),args,value=TRUE)
  if(length(x)) sub(paste0('^--',key,'='),'',x[1]) else default
}
B <- as.integer(arg('replicates',2000)); seed <- as.integer(arg('seed',20261105))
stopifnot(!is.na(B),B>=2,!is.na(seed),seed>0,seed<.Machine$integer.max-20000)
D <- lapply(matrices,jaccard_matrix)
beta <- vapply(D,function(d) mean(d[lower.tri(d)]),numeric(1))
observed <- setNames(c(beta[2]-beta[1],beta[3]-beta[2],beta[3]-beta[1]),c('1-2','2-3','1-3'))
coords <- design[match(paste0(sites,'_1'),design$survey),c('site','lon','lat')]
stopifnot(identical(as.character(coords$site),sites),!anyNA(coords))
xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(coords,coords=c('lon','lat'),crs=4326),25832))/1000
quadrant <- paste(ifelse(xy[,1]<=median(xy[,1]),'West','East'),
  ifelse(xy[,2]<=median(xy[,2]),'South','North'),sep='-')
quadrant_levels <- sort(unique(quadrant))
# Validation against explicitly expanded occurrences, with unequal block sizes.
validate_weighting <- function(xs,g,m) {
  sizes <- tabulate(g); M <- matrix(m,ncol=1)
  ids <- rep(seq_len(nrow(xs[[1]])),times=m[g])
  for(x in xs) {
    d <- jaccard_matrix(x); A <- block_sums(d,g)
    expanded <- jaccard_matrix(x[ids,,drop=FALSE])
    stopifnot(abs(weighted_beta(A,M,sizes)-mean(expanded[lower.tri(expanded)]))<1e-12)
  }
}
toy <- matrix(c(1,0,1, 1,1,0, 0,1,1, 1,0,0, 0,1,0),ncol=3,byrow=TRUE)
validate_weighting(list(toy,toy[,3:1]),c(1,1,2,3,3),c(2,0,3))
validate_weighting(list(toy),c(1,1,2,3,3),c(1,1,1))
results_root <- file.path(out,'bootstrap-results'); dir.create(results_root,showWarnings=FALSE)
run_dir <- tempfile(pattern=paste0(format(Sys.time(),'%Y%m%dT%H%M%S'),'_B-',B,'_seed-',seed,'_'),tmpdir=results_root)
dir.create(run_dir); writeLines('running',file.path(run_dir,'run-status.txt'))
log_path <- file.path(run_dir,'progress.log')
log <- function(x) { message(x); cat(format(Sys.time()),x,'\n',file=log_path,append=TRUE) }
layouts <- list(); all_draws <- list(); influences <- list(); block_audit <- list(); j<-0L
for(width in c(10,20,40,80)) for(origin in c('unshifted','east','north','both')) {
  j<-j+1L; config <- paste0(width,'km-',origin)
  offset <- c(if(origin %in% c('east','both')) width/2 else 0,
              if(origin %in% c('north','both')) width/2 else 0)
  block_xy <- floor(sweep(xy,2,offset,'-')/width)
  key <- paste(block_xy[,1],block_xy[,2],sep='_')
  g <- match(key,unique(key)); K <- max(g); sizes <- tabulate(g)
  stopifnot(K>=2,length(g)==N,sum(sizes)==N)
  A <- lapply(D,block_sums,g=g)
  stopifnot(max(abs(vapply(A,function(a) weighted_beta(a,matrix(1,K,1),sizes),numeric(1))-beta))<1e-12)
  validate_weighting(matrices,g,rep(c(0,2,1),length.out=K))
  regional <- sapply(quadrant_levels,function(q) tabulate(g[quadrant==q],nbins=K))
  layout <- data.frame(site=sites,lon=coords$lon,lat=coords$lat,config=config,
    width_km=width,origin=origin,block=g,bx=block_xy[,1],by=block_xy[,2],quadrant=quadrant)
  layouts[[config]]<-layout
  block_audit[[config]]<-data.frame(config=config,width_km=width,origin=origin,
    block=1:K,cells=sizes)
  # Leave each complete block out of the observed dataset once.
  Mloo <- matrix(1,K,K); diag(Mloo)<-0
  loo <- sapply(A,weighted_beta,M=Mloo,sizes=sizes)
  influences[[config]]<-data.frame(config=config,block=1:K,cells_removed=sizes,
    delta_12=loo[,2]-loo[,1],delta_23=loo[,3]-loo[,2],delta_13=loo[,3]-loo[,1])
  for(batch in 1:2) {
    draw_seed <- seed+(batch-1L)*10000L+j
    set.seed(draw_seed)
    # Generate once, making RNG independent of numerical calculation batch size.
    M <- rmultinom(B,K,rep(1/K,K))
    stopifnot(all(colSums(M)==K))
    calculate <- function(cols) sapply(A,weighted_beta,M=M[,cols,drop=FALSE],sizes=sizes)
    means <- do.call(rbind,lapply(split(seq_len(B),ceiling(seq_len(B)/200)),calculate))
    # Repeat a subset under a different numerical batch split.
    ii<-seq_len(min(B,7)); stopifnot(max(abs(means[ii,,drop=FALSE]-calculate(ii)))<1e-12)
    set.seed(draw_seed); stopifnot(identical(M,rmultinom(B,K,rep(1/K,K))))
    n <- drop(crossprod(sizes,M)); unique_cells <- drop(crossprod(sizes,M>0))
    shares <- t(regional) %*% M; shares <- sweep(shares,2,n,'/')
    z <- data.frame(config=config,width_km=width,origin=origin,batch=batch,
      draw=seq_len(B),seed=draw_seed,atlas1=means[,1],atlas2=means[,2],atlas3=means[,3],
      'delta_1-2'=means[,2]-means[,1],'delta_2-3'=means[,3]-means[,2],
      'delta_1-3'=means[,3]-means[,1],total_cells=n,unique_cells=unique_cells,
      unique_blocks=colSums(M>0),duplicate_pair_fraction=drop(crossprod(sizes,M*(M-1)))/(n*(n-1)),
      check.names=FALSE)
    for(q in seq_along(quadrant_levels)) z[[paste0('share_',quadrant_levels[q])]]<-shares[q,]
    stopifnot(all(is.finite(as.matrix(z[c('atlas1','atlas2','atlas3','delta_1-2','delta_2-3','delta_1-3')]))),
      all(means>=0 & means<=1),all(n>=2),max(abs(z[['delta_1-2']]+z[['delta_2-3']]-z[['delta_1-3']]))<1e-12)
    all_draws[[paste(config,batch)]]<-z
    log(paste(config,'batch',batch,'completed',B,'replicates;',K,'blocks'))
  }
}
draws<-do.call(rbind,all_draws); rownames(draws)<-NULL
block_audit<-do.call(rbind,block_audit); influences<-do.call(rbind,influences)
summary_table<-bootstrap_summary(draws,observed)
stopifnot(all(summary_table$percentile_low<=summary_table$percentile_high),
  max(abs(summary_table$basic_low-(2*summary_table$observed-summary_table$percentile_high)))<1e-12)
metadata<-list(replicates_per_batch=B,seed=seed,batches=2,projection='EPSG:25832',
  input_code_hashes=bootstrap_hashes(root,out),generated=Sys.time(),session=sessionInfo(),
  scheme='K occupied non-overlapping blocks sampled uniformly with replacement; raw paired histories',
  quadrant_medians_km=apply(xy,2,median),quadrant_levels=quadrant_levels,
  original_quadrant_share=prop.table(table(quadrant)),validation='expanded-matrix, observed, pairing, RNG, batching and interval checks passed')
save(draws,layouts,block_audit,influences,summary_table,metadata,observed,beta,sites,
  file=file.path(run_dir,'bootstrap-results.RData'),compress='xz')
writeLines('complete',file.path(run_dir,'run-status.txt'))
writeLines(basename(run_dir),file.path(results_root,'latest-run.txt'))
log(paste('Completed',nrow(draws),'replicates. Saved:',run_dir))
