# Posterior uncertainty in Tjur R2-scaled VP, November atlas models.
# Run from the repository root:
#   Rscript scripts/4-post/compare-scaled-vp-changes.R
#   Rscript scripts/4-post/compare-scaled-vp-changes.R --spatial
#   Rscript scripts/4-post/compare-scaled-vp-changes.R --self-test
# Outputs: three atlas RDS files, comparison.rds, summary.md. No figures/CSVs.
# All scaled VP levels are percentages; temporal changes are percentage points.

vp_diagnostic_thresholds <- list(rhat = 1.1, ess_per_chain = 100, mcse_ratio = .05,
  rhat_strict = TRUE)

vp_options <- function(args) {
  opt <- list(model_root=Sys.getenv("CH1_MODEL_OUTPUT_DIR","data-new/2-hmsc/derived/models"),
    model_prefix="2026-11-06_05-11-30",
    model_suffix="MinOccs5_CoverageGoodAverage_EnvironmentOnly",
    samples=250L,thin=100L,chains=4L,seed=369L,contrast_draws=10000L,
    output_dir="scripts/4-post/compare-scaled-vp-changes",self_test=FALSE,spatial=FALSE)
  if ("--spatial" %in% args) {
    opt$spatial <- TRUE
    opt$model_prefix <- "2026-11-02_06-55-21"
    opt$model_suffix <- "MinOccs5_CoverageGoodAverage"
    opt$output_dir <- file.path(opt$output_dir,"spatial")
    args <- args[args!="--spatial"]
  }
  i <- 1L
  while(i<=length(args)) {
    if(args[i]=="--self-test") {opt$self_test<-TRUE;i<-i+1L;next}
    key <- gsub("-","_",sub("^--","",args[i]))
    if(!startsWith(args[i],"--") || !key %in% setdiff(names(opt),c("self_test","spatial")) || i==length(args))
      stop("Unknown option or missing value: ",args[i])
    opt[[key]] <- args[i+1L];i<-i+2L
  }
  for(k in c("samples","thin","chains","seed","contrast_draws")) {
    v<-suppressWarnings(as.numeric(opt[[k]]))
    if(length(v)!=1L || !is.finite(v) || v<1 || v!=floor(v) || v>.Machine$integer.max)
      stop(k," must be a positive integer")
    opt[[k]]<-as.integer(v)
  }
  if(opt$samples<4L || opt$chains<2L) stop("At least four draws per chain and two chains required")
  if(opt$contrast_draws %% (opt$samples*opt$chains)!=0L)
    stop("contrast-draws must be a multiple of samples * chains for balanced reuse")
  opt
}

vp_require <- function(ok,message) {if(!isTRUE(ok)) stop(message,call.=FALSE)}
vp_summaries <- function(x,directional=TRUE) {
  x<-as.matrix(x)
  q<-t(apply(x,2,quantile,probs=c(.025,.5,.975),names=FALSE))
  ans<-data.frame(mean=colMeans(x),median=q[,2],lower=q[,1],upper=q[,3],
    p_increase=colMeans(x>0),p_decrease=colMeans(x<0),p_zero=colMeans(x==0),row.names=NULL)
  if(!directional) ans<-ans[,c("mean","median","lower","upper")]
  ans
}
vp_flat <- function(a) matrix(a,nrow=prod(dim(a)[1:2]))
vp_metrics <- function(x) {
  x<-as.matrix(x)
  q<-t(apply(x,1,quantile,probs=c(.1,.25,.5,.75,.9),names=FALSE))
  colnames(q)<-c("q10","q25","median","q75","q90")
  cbind(q,iqr=q[,4]-q[,2])
}
vp_cdf <- function(sorted,grid) {
  t(vapply(seq_len(nrow(sorted)),function(i)
    findInterval(grid,sorted[i,])/ncol(sorted),numeric(length(grid))))
}
vp_balanced_indices <- function(n,B,seed) {
  vp_require(B%%n==0L,"Balanced contrast draw count must be a multiple of source draws")
  set.seed(seed)
  ans<-vapply(1:3,function(a) unlist(lapply(seq_len(B/n),function(j)sample.int(n)),
    use.names=FALSE),integer(B))
  colnames(ans)<-paste0("atlas",1:3)
  vp_require(all(apply(ans,2,function(x)all(tabulate(x,nbins=n)==B/n))),"Unbalanced draw indices")
  ans
}
vp_joint_interval <- function(x) {
  x<-as.matrix(x);centre<-colMeans(x);s<-apply(x,2,sd)
  scale<-ifelse(s==0,1,s)
  z<-sweep(sweep(x,2,centre,"-"),2,scale,"/")
  max_deviation<-apply(abs(z),1,max)
  critical<-unname(quantile(max_deviation,.95,type=1))
  coverage<-mean(max_deviation<=critical)
  vp_require(coverage>=.95-1e-10,"Joint envelope covers fewer than 95% of sampled vectors")
  list(lower=centre-critical*s,upper=centre+critical*s,critical=critical,
    empirical_coverage=coverage,centre=centre,sd=s)
}
vp_one_draw <- function(x,y,beta) {
  vp_require(ncol(x)==nrow(beta) && ncol(y)==ncol(beta) && nrow(x)==nrow(y),
    "Incompatible X, Y and Beta dimensions")
  vp_require(all(is.finite(x)) && all(is.finite(beta)),"Nonfinite X or Beta")
  vp_require(all(is.na(y) | y %in% c(0,1)),"Y must be binary or missing")
  vp_require(identical(colnames(y),colnames(beta)),"Species order mismatch")
  vp_require(identical(colnames(x),rownames(beta)),"Coefficient order mismatch")
  vp_require(colnames(x)[1]=="(Intercept)" && all(x[,1]==1),"Constant intercept must be first")
  present<-y==1;absent<-y==0;present[is.na(present)]<-FALSE;absent[is.na(absent)]<-FALSE
  n1<-colSums(present);n0<-colSums(absent)
  vp_require(all(n1>0 & n0>0),"Every species needs recorded presences and absences")
  probability<-pnorm(x%*%beta)
  tjur<-colSums(probability*present)/n1-colSums(probability*absent)/n0
  parts<-sweep(beta[-1,,drop=FALSE]^2,1,diag(cov(x))[-1],"*")
  denominator<-colSums(parts)
  vp_require(all(is.finite(denominator) & denominator>0),"Undefined environmental VP denominator")
  relative<-sweep(parts,2,denominator,"/")
  scaled<-sweep(relative,2,tjur,"*")*100
  vp_require(max(abs(colSums(relative)-1))<1e-10 &&
    max(abs(colSums(scaled)-100*tjur))<1e-10,"Draw-level partition sum failure")
  list(relative=relative,scaled=scaled,tjur=tjur,probability=probability)
}

vp_diagnostics <- function(a,key,helpers) {
  a<-array(a,c(dim(a)[1:2],nrow(key)))
  do.call(rbind,lapply(seq_len(nrow(key)),function(j)
    cbind(key[j,,drop=FALSE],mean=mean(a[,,j]),helpers$diagnose_matrix(a[,,j]))))
}
vp_source_status <- function(a,b) ifelse(a=="pass" & b=="pass","pass","source diagnostics flagged/unassessable")

vp_hmsc_expected <- function(m,d) {
  prediction<-getS3method("predict","Hmsc")
  # Some restricted macOS sessions return NA from detectCores(). Keep Hmsc's
  # prediction implementation unchanged except for its serial CPU discovery.
  if(is.na(parallel::detectCores())) environment(prediction)<-list2env(
    list(detectCores=function(...)1L),parent=environment(prediction))
  # Repeating the draw preserves the draw-by-level array for one random level.
  post <- if (m$nr > 0L) list(d,d) else list(d)
  prediction(m,post=post,expected=TRUE,nParallel=1)[[1]]
}

vp_read_atlas <- function(file,atlas,opt,helpers,reference=NULL) {
  started<-proc.time()[[3]]
  e<-new.env(parent=emptyenv());load(file,e)
  vp_require(exists("fitted_model",e,inherits=FALSE),paste("Missing fitted_model:",file))
  m<-e$fitted_model$posteriors
  vp_require(is.list(m) && !is.null(m$postList),"Expected a merged Hmsc list posterior")
  vp_require(m$nr==as.integer(opt$spatial) && m$ncRRR==0 && is.matrix(m$X) && all(m$distr[,1]==2),
    "Expected a matrix-X binary-probit model with the requested spatial structure")
  if (opt$spatial) {
    vp_require(identical(m$rLNames,"site") && m$ranLevels[[1]]$xDim==0L,
      "Spatial mode requires a single site level with no covariate-dependent loadings")
    site_levels <- levels(m$dfPi$site)
    site_index <- match(as.character(m$studyDesign$site),site_levels)
    vp_require(!anyNA(site_index) && identical(site_index,as.integer(m$Pi[,1])),
      "Spatial sites and fitted design are misaligned")
  }
  vp_require(is.null(m$Loff) || all(m$Loff==0),"Nonzero prediction offsets are unsupported")
  vp_require(length(m$postList)==opt$chains && all(lengths(m$postList)==opt$samples),
    "Saved chain/draw counts do not match requested sampling settings")
  vp_require(isTRUE(all.equal(as.numeric(m$thin),as.numeric(opt$thin))),"Saved thinning label mismatch")
  species<-m$spNames;variables<-setdiff(m$covNames,"(Intercept)")
  vp_require(length(species)==155L && length(variables)==10L &&
    !anyDuplicated(species) && !anyDuplicated(variables),"Expected 155 unique species and ten variables")
  vp_require(identical(colnames(m$X),m$covNames) && identical(colnames(m$Y),species) &&
    identical(rownames(m$X),rownames(m$Y)),"Saved design/response names do not align")
  assignment<-attr(m$X,"assign")
  vp_require(length(assignment)==11L && assignment[1]==0 &&
    identical(as.integer(assignment[-1]),seq_len(10L)),"Expected one simple predictor per VP group")
  if(!is.null(reference)) {
    vp_require(setequal(species,reference$species) && setequal(variables,reference$variables),
      "Species or predictor sets differ among atlases")
    species<-reference$species;variables<-reference$variables
  }
  sites<-as.character(m$studyDesign$site)
  vp_require(length(sites)==nrow(m$X) && !anyNA(sites) && !anyDuplicated(sites),"Invalid fitted-cell identifiers")
  ni<-opt$samples;nc<-opt$chains;ns<-length(species);nv<-length(variables)
  ids<-expand.grid(iteration=seq_len(ni),chain=seq_len(nc))
  scaled<-array(NA_real_,c(ni,nc,ns,nv),dimnames=list(iteration=seq_len(ni),chain=seq_len(nc),
    species=species,variable=variables))
  tjur<-array(NA_real_,c(ni,nc,ns),dimnames=dimnames(scaled)[1:3])
  relative_sum<-matrix(0,nv,ns,dimnames=list(variables,species))
  spatial_community <- matrix(0,ni,nc)
  selected<-rbind(c(1L,1L),c(ceiling(ni/2),min(2L,nc)),c(ni,nc))
  validation<-list(prediction_max_error=0,draw_vp_max_error=0)
  # Hoist design covariance and presence masks outside the full-draw loop.
  design_covariance <- cov(m$X)
  variance<-diag(design_covariance)[-1]
  present<-m$Y==1;absent<-m$Y==0;present[is.na(present)]<-FALSE;absent[is.na(absent)]<-FALSE
  vp_require(all(is.na(m$Y) | m$Y %in% c(0,1)),"Nonbinary responses")
  n1<-colSums(present);n0<-colSums(absent)
  vp_require(all(n1>0 & n0>0),"Every species needs recorded presences and absences")
  for(k in seq_len(nc)) {
    message("Atlas ",atlas,": calculating chain ",k," / ",nc)
    for(i in seq_len(ni)) {
      d<-m$postList[[k]][[i]]
      vp_require(identical(dim(d$Beta),c(11L,155L)) && all(is.finite(d$Beta)) &&
        all(is.finite(d$sigma)) && all(abs(d$sigma-1)<1e-10),"Invalid Beta or non-unit probit residual variance")
      linear_predictor <- m$X %*% d$Beta
      fixed_share <- rep(1,ns)
      if (opt$spatial) {
        eta <- d$Eta[[1]]
        lambda <- d$Lambda[[1]]
        vp_require(nrow(eta)==length(site_levels) && ncol(eta)==nrow(lambda) && ncol(lambda)==ns,
          "Spatial factors and loadings have incompatible dimensions")
        if (!is.null(rownames(eta))) eta <- eta[site_levels,,drop=FALSE]
        linear_predictor <- linear_predictor + eta[site_index,,drop=FALSE] %*% lambda
        fixed <- colSums(d$Beta * (design_covariance %*% d$Beta))
        random <- colSums(lambda^2)
        vp_require(all(is.finite(fixed+random) & fixed+random>0),"Invalid fixed/spatial VP denominator")
        fixed_share <- fixed/(fixed+random)
      }
      p<-pnorm(linear_predictor)
      t<-colSums(p*present)/n1-colSums(p*absent)/n0
      spatial_community[i,k] <- mean(t * (1-fixed_share)) * 100
      parts<-sweep(d$Beta[-1,,drop=FALSE]^2,1,variance,"*")
      denominator<-colSums(parts)
      vp_require(all(is.finite(denominator) & denominator>0),"Undefined environmental VP denominator")
      rel<-sweep(parts,2,denominator,"/")
      if (opt$spatial) rel <- sweep(rel,2,fixed_share,"*")
      vp_require(max(abs(colSums(rel)-fixed_share))<1e-10,"Environmental VP sum failure")
      dimnames(rel)<-list(setdiff(m$covNames,"(Intercept)"),m$spNames)
      rel<-rel[variables,species,drop=FALSE];t<-t[match(species,m$spNames)]
      scaled[i,k,,]<-t(sweep(rel,2,t,"*"))*100
      tjur[i,k,]<-t;relative_sum<-relative_sum+rel
      if(any(selected[,1]==i & selected[,2]==k)) {
        expected<-vp_hmsc_expected(m,d)
        err<-max(abs(expected-p));vp_require(err<1e-10,"Expected probabilities differ from Hmsc")
        validation$prediction_max_error<-max(validation$prediction_max_error,err)
        one<-m;one$postList<-list(list(d))
        official<-suppressWarnings(Hmsc::computeVariancePartitioning(one))$vals[variables,species,drop=FALSE]
        err<-max(abs(official-rel));vp_require(err<1e-10,"Draw VP differs from Hmsc")
        validation$draw_vp_max_error<-max(validation$draw_vp_max_error,err)
      }
    }
  }
  relative_mean<-relative_sum/(ni*nc)
  csv<-file.path(dirname(dirname(dirname(file))),"Results",paste0(basename(dirname(dirname(dirname(file)))),"parameter_estimates_VP_.csv"))
  vp_require(file.exists(csv),paste("Missing legacy VP export:",csv))
  raw<-read.csv(csv,row.names=1,check.names=FALSE)
  vp_require(setequal(colnames(raw),species) && all(c(variables,"TjurR2") %in% rownames(raw)),"Legacy export mismatch")
  export_difference <- max(abs(relative_mean-as.matrix(raw[variables,species,drop=FALSE])))
  if (!opt$spatial) {
    validation$mean_relative_export_max_error <- export_difference
    vp_require(export_difference<1e-8,"Mean draw-level relative VP does not reproduce export")
  }
  flat<-array(scaled,c(ni*nc,ns,nv))
  community<-array(NA_real_,c(ni,nc,nv),dimnames=dimnames(scaled)[c(1,2,4)])
  metrics<-c("q10","q25","median","q75","q90","iqr")
  distributions<-array(NA_real_,c(ni,nc,length(metrics),nv),
    dimnames=list(iteration=seq_len(ni),chain=seq_len(nc),metric=metrics,variable=variables))
  for(v in seq_len(nv)) {
    community[,,v]<-matrix(rowMeans(flat[,,v]),ni,nc)
    distributions[,,,v]<-array(vp_metrics(flat[,,v]),c(ni,nc,length(metrics)))
  }
  total_error<-max(abs(rowSums(vp_flat(community))+as.vector(spatial_community)-rowMeans(vp_flat(tjur))*100))
  vp_require(total_error<1e-10,"Scaled community totals do not equal mean Tjur")
  validation$community_total_max_error<-total_error
  mean_product<-t(matrix(colMeans(vp_flat(scaled)),ns,nv))
  product_means<-sweep(relative_mean,2,colMeans(vp_flat(tjur)),"*")*100
  legacy_product<-sweep(as.matrix(raw[variables,species,drop=FALSE]),2,as.numeric(raw["TjurR2",species]),"*")*100
  scaling_audit<-data.frame(variable=variables,posterior_mean_scaled_percent=rowMeans(mean_product),
    product_posterior_means_percent=rowMeans(product_means),legacy_figure_percent=rowMeans(legacy_product))
  message("Atlas ",atlas,": diagnosing scaled VP, fit, community and distribution summaries")
  scaled_key<-expand.grid(species=species,variable=variables,stringsAsFactors=FALSE)
  distribution_key<-expand.grid(metric=metrics,variable=variables,stringsAsFactors=FALSE)
  ans<-list(atlas=atlas,domain="all fitted cells of this atlas",species=species,variables=variables,
    sites=sites,survey_ids=rownames(m$X),n_cells=nrow(m$X),draw_index=ids,
    scaled_vp_percent=scaled,tjur=tjur,community_percent=community,distribution_percent=distributions,
    scaled_summary=vp_diagnostics(scaled,scaled_key,helpers),
    tjur_summary=vp_diagnostics(tjur,data.frame(species=species),helpers),
    community_summary=vp_diagnostics(community,data.frame(variable=variables),helpers),
    distribution_summary=vp_diagnostics(distributions,distribution_key,helpers),
    negative_tjur_count=sum(tjur<0),negative_tjur_fraction=mean(tjur<0),
    scaling_audit=scaling_audit,validation=validation,
    provenance=list(input=normalizePath(file),input_md5=unname(tools::md5sum(file)),
      export=normalizePath(csv),export_md5=unname(tools::md5sum(csv)),formula=deparse(m$XFormula),
      XScalePar=m$XScalePar,settings=opt,units="scaled VP in percent; Tjur in fraction",
      elapsed_seconds=proc.time()[[3]]-started))
  if (opt$spatial) {
    ans$spatial_community_percent <- spatial_community
    # Hmsc's exported spatial-fit environmental VP multiplies the separately
    # averaged fixed share and within-fixed allocation.
    ans$mean_relative_export_max_difference <- export_difference
  }
  ans
}

vp_compare <- function(atlases,opt) {
  ns<-length(atlases[[1]]$species);nv<-length(atlases[[1]]$variables)
  n<-nrow(atlases[[1]]$draw_index);B<-opt$contrast_draws
  indices<-vp_balanced_indices(n,B,opt$seed)
  pairs<-list(c(1L,2L),c(2L,3L),c(1L,3L));labels<-c("1970s to 1990s","1990s to 2010s","1970s to 2010s")
  variables<-atlases[[1]]$variables;species<-atlases[[1]]$species
  grid<-seq(min(vapply(atlases,function(a)min(a$scaled_vp_percent),numeric(1))),
    max(vapply(atlases,function(a)max(a$scaled_vp_percent),numeric(1))),length.out=201)
  community<-matrix(NA_real_,B,3L*nv)
  community_key<-expand.grid(variable=variables,pair=labels,stringsAsFactors=FALSE)
  community_key$source_diagnostic_status<-NA_character_
  distribution_draws<-matrix(NA_real_,B,3L*nv*6L)
  distribution_rows<-list();species_rows<-list();cdf_rows<-list();w_rows<-list()
  w_draws<-matrix(NA_real_,B,3L*nv)
  for(a in 1:3) atlases[[a]]$cdf_summary<-NULL
  for(v in seq_len(nv)) {
    message("Comparing variable ",v," / ",nv,": ",variables[v])
    native<-lapply(atlases,function(a)array(a$scaled_vp_percent,c(n,ns,nv))[,,v])
    sorted<-lapply(native,function(x)t(apply(x,1,sort)))
    cdfs<-lapply(sorted,vp_cdf,grid=grid)
    metrics<-lapply(native,vp_metrics)
    for(a in 1:3) atlases[[a]]$cdf_summary<-rbind(atlases[[a]]$cdf_summary,
      cbind(data.frame(variable=variables[v],grid_percent=grid),vp_summaries(cdfs[[a]],directional=FALSE)))
    for(pair_id in 1:3) {
      pair<-pairs[[pair_id]];early<-pair[1];late<-pair[2];j<-(pair_id-1L)*nv+v
      ie<-indices[,early];il<-indices[,late]
      community[,j]<-rowMeans(native[[late]][il,,drop=FALSE])-rowMeans(native[[early]][ie,,drop=FALSE])
      community_key$source_diagnostic_status[j]<-vp_source_status(
        atlases[[early]]$community_summary$status[v],atlases[[late]]$community_summary$status[v])
      delta<-native[[late]][il,,drop=FALSE]-native[[early]][ie,,drop=FALSE]
      ss<-cbind(data.frame(pair=labels[pair_id],variable=variables[v],species=species),vp_summaries(delta))
      source_rows<-((v-1L)*ns+1L):(v*ns)
      ss$source_diagnostic_status<-vp_source_status(atlases[[early]]$scaled_summary$status[source_rows],
        atlases[[late]]$scaled_summary$status[source_rows])
      species_rows[[j]]<-ss
      md<-metrics[[late]][il,,drop=FALSE]-metrics[[early]][ie,,drop=FALSE]
      cols<-((j-1L)*6L+1L):(j*6L);distribution_draws[,cols]<-md
      ds<-cbind(data.frame(pair=labels[pair_id],variable=variables[v],metric=colnames(md)),vp_summaries(md))
      source_rows<-((v-1L)*6L+1L):(v*6L)
      ds$source_diagnostic_status<-vp_source_status(atlases[[early]]$distribution_summary$status[source_rows],
        atlases[[late]]$distribution_summary$status[source_rows])
      distribution_rows[[j]]<-ds
      cd<-cdfs[[late]][il,,drop=FALSE]-cdfs[[early]][ie,,drop=FALSE]
      cdf_rows[[j]]<-cbind(data.frame(pair=labels[pair_id],variable=variables[v],grid_percent=grid),vp_summaries(cd))
      # Equal-weight empirical distributions: W1 is mean absolute difference of sorted values.
      wd<-rowMeans(abs(sorted[[late]][il,,drop=FALSE]-sorted[[early]][ie,,drop=FALSE]))
      w_draws[,j]<-wd
      w_rows[[j]]<-cbind(data.frame(pair=labels[pair_id],variable=variables[v]),vp_summaries(wd,directional=FALSE))
    }
  }
  joint<-vp_joint_interval(community)
  community_summary<-cbind(community_key,vp_summaries(community),
    joint_lower=joint$lower,joint_upper=joint$upper)
  vp_require(nrow(community_summary)==30L,"Expected 30 community contrasts")
  vp_require(max(abs(community[,2L*nv+seq_len(nv)]-
    community[,seq_len(nv)]-community[,nv+seq_len(nv)]))<1e-10,"Overlapping temporal contrasts inconsistent")
  list(atlases=atlases,comparison=list(draw_triplets=indices,community_draws_pp=community,
    community_summary=community_summary,joint_envelope=joint,
    species_summary=do.call(rbind,species_rows),distribution_draws_pp=distribution_draws,
    distribution_summary=do.call(rbind,distribution_rows),cdf_grid_percent=grid,
    cdf_difference_summary=do.call(rbind,cdf_rows),wasserstein_draws_pp=w_draws,
    wasserstein_summary=do.call(rbind,w_rows)))
}

vp_write_summary <- function(atlases,result,opt,runtime,file) {
  fmt<-function(x)sprintf("%.3f",x)
  table_md<-function(x) {
    c(paste0("| ",paste(names(x),collapse=" | ")," |"),
      paste0("| ",paste(rep("---",ncol(x)),collapse=" | ")," |"),
      apply(x,1,function(z)paste0("| ",paste(z,collapse=" | ")," |")))
  }
  cs<-result$community_summary
  tab<-data.frame(Pair=cs$pair,Variable=cs$variable,Mean_pp=fmt(cs$mean),
    CI95=paste0("[",fmt(cs$lower),", ",fmt(cs$upper),"]"),
    Joint95=paste0("[",fmt(cs$joint_lower),", ",fmt(cs$joint_upper),"]"),
    Pr_increase=fmt(cs$p_increase),Diagnostics=cs$source_diagnostic_status)
  diag_rows<-do.call(rbind,lapply(atlases,function(a) {
    do.call(rbind,lapply(c("scaled","tjur","community","distribution"),function(k) {
      d<-a[[paste0(k,"_summary")]]
      data.frame(Atlas=a$atlas,Quantity=k,Pass=sum(d$status=="pass"),
        Flagged=sum(d$status=="flagged"),Unassessable=sum(d$status=="cannot assess"))
    }))
  }))
  audit<-do.call(rbind,lapply(atlases,function(a)cbind(Atlas=a$atlas,a$scaling_audit)))
  for(j in 3:ncol(audit)) audit[,j]<-fmt(audit[,j])
  lines<-c("# Posterior changes in Tjur R²-scaled VP","",
    sprintf("Completed %s. Local runtime: %.1f seconds (%.2f minutes).",Sys.time(),runtime,runtime/60),"",
    "## What was calculated","",
    sprintf("Model family: %s (%s).",if (opt$spatial) "spatial" else "environment-only",opt$model_prefix),"",
    "All 155 species, ten predictors, four chains and every saved draw are retained. Each atlas uses its own fitted cells; there is no common-cell restriction or model refitting.","",
    "For each draw: compute Hmsc environmental VP; compute expected occurrence probabilities from that same draw; calculate species Tjur R²; multiply VP by Tjur R². Values are percentages and changes are percentage points. Negative draw-level fit values are retained.","",
    if (opt$spatial) c("Predictions include the fitted site effect. Environmental shares include the site contribution in the partition denominator. Within each draw, the environmental and site shares sum to the full Tjur R² after scaling. The legacy spatial VP export averages the fixed share and its allocation across variables separately before multiplying them.",""),
    paste0("Cells by atlas: ",paste(vapply(atlases,function(a)a$n_cells,integer(1)),collapse=", "),"."),
    paste0("Negative Tjur draw/species entries by atlas: ",paste(vapply(atlases,function(a)a$negative_tjur_count,integer(1)),collapse=", "),"."),"",
    sprintf("Independently combine whole atlas posterior draws using %d balanced draw triplets, seed %d. Reuse triplets across all contrasts; retain dependence among species and variables. These triplets do not create new observations or increase the original MCMC effective sample size.",opt$contrast_draws,opt$seed),"",
    "## Community comparisons","",
    "Differences are later minus earlier. CI95 is a marginal credible interval. Joint95 is a joint posterior envelope across all 30 community contrasts, not a Holm adjustment or frequentist familywise-error guarantee.","",
    sprintf("Joint-envelope empirical coverage: %.4f. %d/30 marginal intervals and %d/30 joint intervals exclude zero.",
      result$joint_envelope$empirical_coverage,sum(cs$lower>0 | cs$upper<0),sum(cs$joint_lower>0 | cs$joint_upper<0)),
    "These counts must be interpreted alongside the diagnostic flags; flagged source posteriors undermine interval reliability.","",table_md(tab),"",
    "## Convergence and numerical uncertainty","",
    "Diagnostics use the original chain arrays. The checks require rank-normalised R-hat < 1.1, bulk and tail ESS of at least 100 per chain, and MCSE of the mean at most 5% of posterior SD. Community and distribution contrasts carry the diagnostic status of the contributing atlas summaries.","",table_md(diag_rows),"",
    "Directional probabilities have finite MCMC and draw-pairing error. They are not p-values; 1 - Pr(increase) is not sent to Holm. Native MCMC MCSE is recorded in atlas summaries. Repeated contrast triplets do not remedy insufficient effective sample sizes.","",
    "## Product-of-means versus mean-of-products","",
    "The new posterior mean of VP × Tjur can differ from mean(VP) × mean(Tjur), because VP and fit depend on the same parameters. The legacy figure value is included separately; differences can also reflect the legacy fit estimate. This script does not overwrite that figure.","",table_md(audit),"",
    "## Species and distribution results","",
    "comparison.rds contains species_summary (4,650 paired species contrasts), distribution_summary (quantile and IQR changes), cdf_difference_summary (pointwise bands for cumulative-distribution differences), and wasserstein_summary (descriptive distributional distance). Atlas RDS files contain cumulative-distribution summaries and full scaled VP draws.","",
    "Quantile changes show whether importance shifted broadly or mainly in a tail; IQR changes describe spread. Cumulative curves show the proportion of species at or below a given VP value. CDF bands are pointwise, not simultaneous over the grid or variables, and CDF-specific convergence has not been assessed. A positive Wasserstein distance is not evidence by itself against a no-change null; continuous uncertain posteriors can produce positive distances even for indistinguishable distributions.","",
    "## What this does and does not establish","",
    "This quantifies uncertainty in changes in explanatory importance conditional on the fitted models and observed data, including dependence represented by the Hmsc posterior. It does not prove causality, isolate changing response slopes, correct spatial confounding, account for every source of survey error, or establish performance at new sites. Environmental variation, fitted cell coverage, species distributions, response estimates and overall fit can all affect scaled VP.","",
    "## Reproducibility","","```r",
    sprintf("x <- readRDS(%s)",encodeString(file.path(opt$output_dir,"comparison.rds"),quote='"')),
    "x$community_summary","x$species_summary","x$provenance","```","",
    "Run `Rscript scripts/4-post/compare-scaled-vp-changes.R --self-test` for synthetic checks. Use `--spatial` for the 2 November spatial fits or run without arguments for the 6 November environment-only fits. Options include --model-root, --model-prefix, --model-suffix, --samples, --thin, --chains, --seed, --contrast-draws and --output-dir. Model-root defaults to CH1_MODEL_OUTPUT_DIR when set. Contrast draws must be a multiple of the retained draw count.")
  writeLines(lines,file,useBytes=TRUE)
}

vp_self_test <- function() {
  set.seed(5)
  x<-cbind("(Intercept)"=1,a=seq(-2,2,length.out=30),b=rnorm(30))
  b<-matrix(c(0,1,.2,0,-1,.4),3,2,dimnames=list(colnames(x),c("sp1","sp2")))
  y<-cbind(sp1=as.numeric(x[,"a"]>0),sp2=as.numeric(x[,"a"]<0))
  z<-vp_one_draw(x,y,b)
  vp_require(max(abs(colSums(z$scaled)-100*z$tjur))<1e-10,"Synthetic sums")
  reverse<-vp_one_draw(x,y[,2:1],b[,2:1]);vp_require(max(abs(z$scaled-reverse$scaled[,2:1]))<1e-10,"Species rematching")
  negative<-vp_one_draw(x,y,-b);vp_require(all(negative$tjur<0) && all(negative$scaled<0),"Negative fit must remain negative")
  expect_error<-function(expr)vp_require(inherits(tryCatch({force(expr);NULL},error=identity),"error"),"Invalid input did not fail")
  expect_error(vp_one_draw(x,y,b*0));expect_error(vp_one_draw(x,y[,2:1],b))
  expect_error(vp_one_draw(x,y*0,b));expect_error(vp_one_draw(x,y+2,b))
  expect_error(vp_options(c("--contrast-draws","10001")))
  ii<-vp_balanced_indices(40,400,369);vp_require(identical(ii,vp_balanced_indices(40,400,369)),"Reproducible triplets")
  vp_require(!identical(ii[,1],ii[,2]),"Independent atlas schedules")
  xx<-matrix(rnorm(1000),100,10);delta<-(xx+2)-xx
  vp_require(max(abs(vp_summaries(delta)$mean-2))<1e-12,"Known location shift")
  vp_require(max(abs(vp_metrics(xx+2)[,1:5]-vp_metrics(xx)[,1:5]-2))<1e-12,"Known quantile shift")
  sx<-t(apply(xx,1,sort));grid<-seq(-5,5,length.out=201)
  vp_require(max(abs(vp_cdf(sx,grid)-vp_cdf(sx,grid)))==0 && all(rowMeans(abs(sx-sx))==0),"Identical distributions")
  vp_require(all(vp_summaries(xx-xx)$lower==0),"Identical paired inputs")
  zero<-vp_joint_interval(matrix(0,100,30));vp_require(all(zero$lower==0) && zero$empirical_coverage==1,"Degenerate joint envelope")
  joint<-vp_joint_interval(matrix(rnorm(30000),1000,30));vp_require(joint$empirical_coverage>=.95,"Joint coverage")
  message("All synthetic self-tests passed.")
  invisible(TRUE)
}

vp_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  started<-proc.time()[[3]];opt<-vp_options(args)
  if(opt$self_test) return(vp_self_test())
  for(package in c("Hmsc","posterior"))
    if(!requireNamespace(package,quietly=TRUE)) stop("Required package: ",package)
  helpers<-new.env();helper_path<-"scripts/2-hmsc/model-diagnostics/helpers.R";source(helper_path,local=helpers)
  helpers$diag_thresholds <- vp_diagnostic_thresholds
  folders<-sprintf("%s_Atlas%d_%s",opt$model_prefix,1:3,opt$model_suffix)
  files<-file.path(opt$model_root,folders,"Models","Fitted",
    sprintf("HPC_samples_%04d_thin_%02d_chains_%d.Rdata",opt$samples,opt$thin,opt$chains))
  vp_require(all(file.exists(files)),paste("Missing model(s):",paste(files[!file.exists(files)],collapse=", ")))
  hashes<-tools::md5sum(files)
  atlases<-list()
  for(a in 1:3) atlases[[a]]<-vp_read_atlas(files[a],a,opt,helpers,if(a>1)atlases[[1]] else NULL)
  results<-vp_compare(atlases,opt);atlases<-results$atlases;comparison<-results$comparison
  vp_require(identical(hashes,tools::md5sum(files)),"Input model files changed during analysis")
  runtime<-proc.time()[[3]]-started
  comparison$provenance<-list(settings=opt,created=as.character(Sys.time()),compute_seconds=runtime,
    models=normalizePath(files),model_md5=hashes,script_md5=tools::md5sum(c(
      "scripts/4-post/compare-scaled-vp-changes.R",helper_path)),
    versions=vapply(c("Hmsc","posterior"),function(p)as.character(utils::packageVersion(p)),character(1)),
    R=R.version.string,session=capture.output(sessionInfo()),diagnostic_thresholds=helpers$diag_thresholds,
    units="VP levels in percent; changes and W1 in percentage points; Tjur in fraction; CDF in fraction",
    inference="independent fitted-atlas posteriors; fixed species set; all fitted cells per atlas",
    community_diagnostics_all_pass=all(vapply(atlases,function(a)all(a$community_summary$status=="pass"),logical(1))))
  dir.create(opt$output_dir,recursive=TRUE,showWarnings=FALSE)
  for(a in 1:3) saveRDS(atlases[[a]],file.path(opt$output_dir,paste0("atlas",a,".rds")))
  comparison$provenance$runtime_seconds<-proc.time()[[3]]-started
  saveRDS(comparison,file.path(opt$output_dir,"comparison.rds"))
  vp_write_summary(atlases,comparison,opt,proc.time()[[3]]-started,file.path(opt$output_dir,"summary.md"))
  message("Completed in ",round(proc.time()[[3]]-started,1)," seconds. Results: ",normalizePath(opt$output_dir))
  print(comparison$community_summary[,c("pair","variable","mean","lower","upper","joint_lower","joint_upper","p_increase","source_diagnostic_status")],row.names=FALSE)
  invisible(comparison)
}

if(sys.nframe()==0L) vp_main()
