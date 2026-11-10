# Shared, read-only posterior diagnostics. Matrices are iteration x chain.
diag_thresholds <- list(rhat=1.01, ess_per_chain=100, mcse_ratio=.05)
safe_stat <- function(f, x, ...) tryCatch(suppressWarnings(as.numeric(f(x,...))), error=function(e) NA_real_)
diagnose_matrix <- function(z) {
  stopifnot(is.matrix(z))
  valid <- all(is.finite(z))
  q <- if(valid) quantile(z,c(.025,.5,.975),names=FALSE) else rep(NA_real_,3)
  d <- c(rhat=safe_stat(posterior::rhat,z),ess_bulk=safe_stat(posterior::ess_bulk,z),
    ess_tail=safe_stat(posterior::ess_tail,z),mcse_mean=safe_stat(posterior::mcse_mean,z))
  s <- if(valid) sd(as.vector(z)) else NA_real_
  end <- safe_stat(posterior::mcse_quantile,z,probs=c(.025,.975))
  if(length(end)!=2) end<-rep(NA_real_,2)
  ratio <- d['mcse_mean']/s
  flags <- c(rhat=d['rhat']>diag_thresholds$rhat,
    bulk=d['ess_bulk']<diag_thresholds$ess_per_chain*ncol(z),
    tail=d['ess_tail']<diag_thresholds$ess_per_chain*ncol(z),mcse=ratio>diag_thresholds$mcse_ratio)
  cannot <- !valid || any(!is.finite(c(d,ratio,end))) || s==0 || ncol(z)<2 || nrow(z)<4
  data.frame(median=q[2],sd=s,lower=q[1],upper=q[3],width=q[3]-q[1],
    rhat=unname(d[1]),ess_bulk=unname(d[2]),ess_tail=unname(d[3]),mcse_mean=unname(d[4]),
    mcse_sd_ratio=unname(ratio),mcse_lower=end[1],mcse_upper=end[2],
    p_positive=if(valid) mean(z>0) else NA_real_,p_negative=if(valid) mean(z<0) else NA_real_,
    wrong_sign=if(valid) 1-max(mean(z>0),mean(z<0)) else NA_real_,
    spans_zero=if(valid) q[1]<=0 && q[3]>=0 else NA,
    flag_rhat=unname(flags[1]),flag_bulk=unname(flags[2]),flag_tail=unname(flags[3]),flag_mcse=unname(flags[4]),
    status=if(cannot) 'cannot assess' else if(any(flags)) 'flagged' else 'pass',row.names=NULL)
}
chain_summary <- function(z) do.call(rbind,lapply(seq_len(ncol(z)),function(k) {
  v<-z[,k]; q<-if(all(is.finite(v))) quantile(v,c(.025,.5,.975),names=FALSE) else rep(NA,3)
  data.frame(chain=k,mean=mean(v),sd=sd(v),lower=q[1],median=q[2],upper=q[3],
    p_positive=mean(v>0),p_negative=mean(v<0))
}))
extract_beta <- function(m) {
  n<-length(m$postList[[1]]); c<-length(m$postList); p<-length(m$covNames); s<-length(m$spNames)
  stopifnot(n>0,c>0,all(lengths(m$postList)==n),!anyDuplicated(m$spNames),!anyDuplicated(m$covNames),
    is.matrix(m$XScaled),ncol(m$XScaled)==p,ncol(m$Y)==s,identical(colnames(m$Y),m$spNames),
    identical(colnames(m$XScaled),m$covNames),all(is.finite(m$XScaled)))
  a<-array(NA_real_,c(n,c,p*s))
  for(k in seq_len(c)) for(i in seq_len(n)) {
    b<-m$postList[[k]][[i]]$Beta
    stopifnot(identical(dim(b),c(p,s)))
    a[i,k,]<-as.vector(b)
  }
  key<-expand.grid(predictor=m$covNames,species=m$spNames,stringsAsFactors=FALSE)
  key$parameter<-seq_len(nrow(key)); key$intercept<-key$predictor=='(Intercept)'
  list(draws=a,key=key)
}
summarise_array <- function(a,key) do.call(rbind,lapply(seq_len(dim(a)[3]),function(j) cbind(key[j,,drop=FALSE],diagnose_matrix(a[,,j]))))
write_table <- function(x,name,out) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='NA')
save_plot <- function(p,name,out,w=12,h=7) ggplot2::ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=140,bg='white',limitsize=FALSE)
plot_beta <- function(d,a,out) {
  g<-ggplot2::ggplot; aes<-ggplot2::aes
  for(metric in c('rhat','ess_bulk','ess_tail')) for(part in c(FALSE,TRUE)) {
    x<-d[d$intercept==part,,drop=FALSE]; if(!nrow(x)) next
    save_plot(g(x,aes(predictor,species,fill=.data[[metric]]))+ggplot2::geom_tile()+
      ggplot2::scale_fill_viridis_c(na.value='grey65')+ggplot2::theme_minimal(base_size=10)+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),axis.text.y=ggplot2::element_text(size=5))+
      ggplot2::labs(title=paste('Beta',metric,if(part) 'intercepts' else 'environmental slopes'),x=NULL,y=NULL),
      paste0('beta-',metric,if(part) '-intercepts' else '-slopes'),out,w=if(part) 6 else 12,h=max(9,nrow(x)/length(unique(x$predictor))*.085))
  }
  dist<-do.call(rbind,lapply(c('rhat','ess_bulk','ess_tail','mcse_sd_ratio'),function(v) data.frame(metric=v,value=d[[v]])))
  save_plot(g(dist,aes(value))+ggplot2::geom_histogram(bins=40,na.rm=TRUE)+ggplot2::facet_wrap(~metric,scales='free')+ggplot2::theme_minimal(), 'diagnostic-distributions',out)
  bad<-unique(c(head(order(d$rhat,decreasing=TRUE,na.last=NA),5),head(order(d$ess_bulk,na.last=NA),5),head(order(d$ess_tail,na.last=NA),5)))
  good<-head(which(d$status=='pass')[order(d$species[d$status=='pass'],d$predictor[d$status=='pass'])],3)
  sel<-unique(c(bad,good));write_table(transform(d[sel,],selection=ifelse(sel %in% bad,'worst diagnostic','passing comparison')),'selected-coefficients',out)
  for(j in sel) {
    z<-a[,,j]; dd<-data.frame(iteration=rep(seq_len(nrow(z)),ncol(z)),chain=factor(rep(seq_len(ncol(z)),each=nrow(z))),value=as.vector(z))
    dd$rank<-rank(dd$value,ties.method='average')
    label<-paste(d$species[j],d$predictor[j]); prefix<-sprintf('coefficient-%04d',j)
    save_plot(g(dd,aes(iteration,value,colour=chain))+ggplot2::geom_line(linewidth=.25)+ggplot2::theme_minimal()+ggplot2::labs(title=label),paste0(prefix,'-trace'),out)
    save_plot(g(dd,aes(rank,fill=chain))+ggplot2::geom_histogram(bins=20)+ggplot2::facet_wrap(~chain,ncol=1)+ggplot2::theme_minimal()+ggplot2::labs(title=label),paste0(prefix,'-ranks'),out)
    save_plot(g(dd,aes(value,colour=chain))+ggplot2::geom_density(na.rm=TRUE)+ggplot2::theme_minimal()+ggplot2::labs(title=label),paste0(prefix,'-density'),out)
    ac<-do.call(rbind,lapply(seq_len(ncol(z)),function(k) {
      v<-z[,k]; aa<-if(all(is.finite(v)) && sd(v)>0) as.numeric(acf(v,plot=FALSE,lag.max=min(50,nrow(z)-1))$acf) else NA_real_
      data.frame(lag=seq_along(aa)-1,acf=aa,chain=factor(k))
    }))
    save_plot(g(ac,aes(lag,acf,colour=chain))+ggplot2::geom_line(na.rm=TRUE)+ggplot2::theme_minimal()+ggplot2::labs(title=label),paste0(prefix,'-acf'),out)
  }
}
add_species <- function(d,m) {
  meta<-data.frame(species=m$spNames,occurrences=colSums(m$Y>0,na.rm=TRUE),n_observed=colSums(!is.na(m$Y)))
  meta$prevalence<-meta$occurrences/meta$n_observed
  tr<-m$TrData
  if(!is.null(tr) && nrow(tr)==nrow(meta)) {
    stopifnot(!is.null(rownames(tr)),setequal(rownames(tr),m$spNames))
    tr<-tr[match(m$spNames,rownames(tr)),,drop=FALSE]
  } else tr<-data.frame(row.names=m$spNames)
  for(v in c('foraging_guild_consensus','Migration_a3_DOF','species_thermal_index')) meta[[v]]<-if(v %in% names(tr)) tr[[v]] else NA
  cbind(d,meta[match(d$species,meta$species),setdiff(names(meta),'species'),drop=FALSE])
}
support_summary <- function(d) do.call(rbind,lapply(c(FALSE,TRUE),function(intercept) do.call(rbind,lapply(c(.85,.95,.975),function(t) {
  x<-d[d$intercept==intercept,,drop=FALSE]; chosen<-x$status=='pass' & pmax(x$p_positive,x$p_negative)>=t;chosen[is.na(chosen)]<-FALSE
  n<-sum(chosen);e<-sum(x$wrong_sign[chosen])
  data.frame(family=if(intercept) 'intercepts' else 'environmental slopes',threshold=t,n_coefficients=nrow(x),n_passing=sum(x$status=='pass'),selected=n,expected_wrong_sign=e,expected_wrong_sign_proportion=if(n) e/n else NA_real_)
}))))
concentration <- function(d,out) {
  d$flagged<-d$status=='flagged';d$cannot_assess<-d$status=='cannot assess';d$weak_support<-pmax(d$p_positive,d$p_negative)<.95
  # Broad uncertainty is descriptive and defined within each predictor.
  d$broad_uncertainty<-ave(d$width,d$predictor,FUN=function(x) x>=quantile(x,.75,na.rm=TRUE))==1
  for(v in c('predictor','species','foraging_guild_consensus','Migration_a3_DOF')) {
    x<-if(v=='predictor') d else d[!d$intercept,,drop=FALSE]
    groups<-split(x,ifelse(is.na(x[[v]]),'Missing',as.character(x[[v]])))
    tab<-do.call(rbind,lapply(names(groups),function(n) {z<-groups[[n]];data.frame(group=n,n_species=length(unique(z$species)),n_coefficients=nrow(z),n_flagged=sum(z$flagged),n_cannot_assess=sum(z$cannot_assess),fraction_flagged=mean(z$flagged),n_rhat=sum(z$flag_rhat,na.rm=TRUE),n_bulk=sum(z$flag_bulk,na.rm=TRUE),n_tail=sum(z$flag_tail,na.rm=TRUE),n_mcse=sum(z$flag_mcse,na.rm=TRUE),median_width=median(z$width,na.rm=TRUE),fraction_spanning_zero=mean(z$spans_zero,na.rm=TRUE))}))
    write_table(tab,paste0('concentration-',v),out)
    if(v!='species') save_plot(ggplot2::ggplot(tab,ggplot2::aes(reorder(group,fraction_flagged),fraction_flagged))+ggplot2::geom_col()+ggplot2::geom_text(ggplot2::aes(label=paste0('s=',n_species,', c=',n_coefficients)),hjust=-.05,size=2.5)+ggplot2::coord_flip(ylim=c(0,1.3))+ggplot2::theme_minimal()+ggplot2::labs(x=v,y='Fraction flagged (all coefficients in group)'),paste0('concentration-',v),out,h=max(6,nrow(tab)*.25))
  }
  for(y in c('rhat','width')) for(x in c('prevalence','species_thermal_index')) save_plot(ggplot2::ggplot(d,ggplot2::aes(.data[[x]],.data[[y]],colour=status))+ggplot2::geom_point(alpha=.45,size=.7,na.rm=TRUE)+ggplot2::facet_wrap(~predictor,scales='free_y')+ggplot2::theme_minimal(),paste0(y,'-vs-',x),out,h=10)
  write_table(as.data.frame(table(d[,c('intercept','status','broad_uncertainty','weak_support')],useNA='ifany')),'uncertainty-support-overlap',out)
  d
}
# Low-rank algebra avoids allocating cell x species x draw prediction arrays.
spatial_cor <- function(x,eta,lambda) {
  xc<-scale(x,center=TRUE,scale=TRUE); ec<-scale(eta,center=TRUE,scale=FALSE)
  variance<-colSums(lambda*(crossprod(ec)%*%lambda))/(nrow(x)-1)
  sweep(crossprod(xc,ec)%*%lambda/(nrow(x)-1),2,sqrt(pmax(variance,0)),'/')
}
spatial_diagnostics <- function(m,b,shared,out) {
  levels<-which(vapply(m$ranLevels,function(r) isTRUE(r$sDim>0),logical(1)))
  if(!length(levels)) return('Inapplicable: no spatial random levels.')
  if(!is.matrix(m$X) || !identical(colnames(m$X),m$covNames)) stop('Spatial projection requires a common matrix design matching saved Beta.')
  sites<-as.character(m$studyDesign$site)
  stopifnot(length(sites)==nrow(m$X),!anyNA(sites),!anyDuplicated(sites),nrow(m$Pi)==length(sites))
  n<-dim(b$draws)[1]; nc<-dim(b$draws)[2]; s<-length(m$spNames); p<-length(m$covNames)
  notes<-character()
  for(r in levels) {
    rn<-names(m$ranLevels)[r]; if(is.null(rn)||is.na(rn)||!nzchar(rn)) rn<-paste0('level',r)
    prefix<-paste0('spatial-',make.names(rn))
    for(scope in c('all-cells','shared-cells')) {
      ix<-if(scope=='all-cells') seq_along(sites) else match(shared,sites)
      if(length(ix)<2) {notes<-c(notes,paste(rn,scope,'Skipped: fewer than two shared cells.'));next}
      stopifnot(!anyNA(ix))
      x<-m$X[ix,,drop=FALSE]; variable<-which(apply(x,2,sd)>0)
      if(!length(variable)) {notes<-c(notes,paste(rn,scope,'No nonconstant predictors.'));next}
      key<-expand.grid(predictor=colnames(x)[variable],species=m$spNames,stringsAsFactors=FALSE)
      corrs<-array(NA_real_,c(n,nc,nrow(key)));projected<-array(NA_real_,dim(b$draws))
      fullrank<-qr(x)$rank==ncol(x)
      fieldrows<-list(); ii<-0L; dimensions<-list()
      for(k in seq_len(nc)) for(i in seq_len(n)) {
        draw<-m$postList[[k]][[i]]; raw<-draw$Eta[[r]];lam<-draw$Lambda[[r]]
        stopifnot(all(m$Pi[ix,r]>=1),all(m$Pi[ix,r]<=nrow(raw)),ncol(raw)==nrow(lam),ncol(lam)==s)
        eta<-raw[m$Pi[ix,r],,drop=FALSE]
        stopifnot(all(is.finite(eta)),all(is.finite(lam)))
        cr<-spatial_cor(x[,variable,drop=FALSE],eta,lam);corrs[i,k,]<-as.vector(cr)
        fr<-suppressWarnings(cor(x[,variable,drop=FALSE],eta))
        ii<-ii+1L
        fieldrows[[ii]]<-data.frame(iteration=i,chain=k,predictor=rep(colnames(x)[variable],ncol(eta)),saved_factor=rep(seq_len(ncol(eta)),each=length(variable)),correlation=as.vector(fr))
        dimensions[[ii]]<-data.frame(chain=k,iteration=i,n_factors=ncol(eta))
        if(fullrank) projected[i,k,]<-as.vector(qr.coef(qr(x),eta)%*%lam)
        if(i %in% c(1,n) && k %in% c(1,nc)) {
          direct<-suppressWarnings(cor(x[,variable,drop=FALSE],eta%*%lam))
          stopifnot(isTRUE(all.equal(cr,direct,tolerance=1e-9,check.attributes=FALSE)))
          rotation<-qr.Q(qr(matrix(sin(seq_len(ncol(eta)^2)),ncol(eta))))
          stopifnot(max(abs(eta%*%lam-(eta%*%rotation)%*%(t(rotation)%*%lam)))<1e-8)
        }
      }
      tag<-paste(prefix,scope,sep='-')
      cs<-summarise_array(corrs,key);cs$n_cells<-length(ix);write_table(cs,paste0(tag,'-combined-correlations'),out)
      ch<-do.call(rbind,lapply(seq_len(nrow(key)),function(j) cbind(key[j,],chain_summary(corrs[,,j]))));ch$n_cells<-length(ix)
      write_table(ch,paste0(tag,'-combined-chains'),out)
      fr<-do.call(rbind,fieldrows);write_table(do.call(rbind,dimensions),paste0(tag,'-factor-dimensions'),out)
      fs<-do.call(rbind,lapply(split(fr,interaction(fr$chain,fr$predictor,fr$saved_factor,drop=TRUE)),function(z) data.frame(chain=z$chain[1],predictor=z$predictor[1],saved_factor=z$saved_factor[1],n_draws=nrow(z),n_cells=length(ix),mean=mean(z$correlation),sd=sd(z$correlation),lower=quantile(z$correlation,.025,na.rm=TRUE),upper=quantile(z$correlation,.975,na.rm=TRUE))))
      write_table(fs,paste0(tag,'-saved-field-chains'),out)
      save_plot(ggplot2::ggplot(fr,ggplot2::aes(correlation,colour=factor(chain)))+ggplot2::geom_density(na.rm=TRUE)+ggplot2::facet_grid(saved_factor~predictor)+ggplot2::theme_minimal(base_size=8)+ggplot2::labs(title=paste(rn,scope,'saved factor basis (unaltered)'),colour='Chain',subtitle='Factor identity and sign are conditional on the saved basis; dimensions may vary.'),paste0(tag,'-saved-fields'),out,w=16,h=10)
      save_plot(ggplot2::ggplot(ch,ggplot2::aes(predictor,median,colour=factor(chain)))+ggplot2::geom_boxplot(outlier.size=.25)+ggplot2::coord_flip()+ggplot2::theme_minimal()+ggplot2::labs(title=paste(rn,scope,'Eta × Lambda correlation'),subtitle=paste(length(ix),'cells; species distributions of chain medians'),colour='Chain'),paste0(tag,'-combined'),out)
      if(scope=='all-cells' && fullrank) {
        ps<-summarise_array(projected,b$key);ts<-summarise_array(projected+b$draws,b$key)
        write_table(ps,paste0(prefix,'-projected'),out);write_table(ts,paste0(prefix,'-beta-plus-projection'),out)
        compensation<-do.call(rbind,lapply(seq_len(nrow(b$key)),function(j) data.frame(b$key[j,],chain=seq_len(nc),correlation=vapply(seq_len(nc),function(k) suppressWarnings(cor(projected[,k,j],b$draws[,k,j])),numeric(1)),row.names=NULL)))
        write_table(compensation,paste0(prefix,'-within-chain-compensation'),out)
        pp<-rbind(data.frame(component='Beta',b$key,rhat=vapply(seq_len(nrow(b$key)),function(j) diagnose_matrix(b$draws[,,j])$rhat,numeric(1))),data.frame(component='Projected space',b$key,rhat=ps$rhat),data.frame(component='Beta + projected space',b$key,rhat=ts$rhat))
        save_plot(ggplot2::ggplot(pp,ggplot2::aes(predictor,rhat,colour=component))+ggplot2::geom_boxplot(outlier.size=.3)+ggplot2::geom_hline(yintercept=1.01,lty=2)+ggplot2::coord_flip()+ggplot2::theme_minimal(),paste0(prefix,'-component-convergence'),out)
      }
      notes<-c(notes,paste(rn,scope,':',length(ix),'cells;',if(fullrank) 'full rank design.' else 'rank deficient: individual projections skipped.'))
    }
  }
  paste(notes,collapse='\n')
}
fit_metrics <- function(m,input,workers) {
  cache<-file.path(dirname(input),sub('HPC_','MF_preCV_',basename(input),fixed=TRUE))
  provenance<-list(input_md5=unname(tools::md5sum(input)),Hmsc=as.character(utils::packageVersion('Hmsc')),mode='in-sample all saved draws expected=TRUE',version=1L)
  use<-FALSE
  if(file.exists(cache)) {e<-new.env(parent=emptyenv());load(cache,e);use<-exists('MF',e,inherits=FALSE)&&exists('MF_provenance',e,inherits=FALSE)&&identical(e$MF_provenance,provenance)}
  if(use) {message('Reusing validated fit cache');return(e$MF)}
  message('Computing in-sample fit with all saved draws')
  preds<-Hmsc::computePredictedValues(m,nParallel=workers)
  MF<-Hmsc::evaluateModelFit(m,predY=preds);rm(preds)
  MF_provenance<-provenance
  save(MF,MF_provenance,file=cache)
  MF
}
read_model <- function(input) {e<-new.env(parent=emptyenv());load(input,e);if(is.null(e$fitted_model$posteriors$postList)) stop('Missing fitted_model$posteriors$postList: ',input);e$fitted_model$posteriors}
print_heatmaps <- function(d,out) {
  species<-sort(unique(d$species));groups<-split(species,ceiling(seq_along(species)/32))
  d$coefficient_type<-ifelse(d$intercept,'Intercept','Environmental slopes')
  d$coefficient_type<-factor(d$coefficient_type,levels=c('Intercept','Environmental slopes'))
  for(metric in c('rhat','ess_bulk','ess_tail')) for(i in seq_along(groups)) {
    z<-d[d$species %in% groups[[i]],,drop=FALSE]
    values<-d[[metric]][is.finite(d[[metric]])]; limits<-if(length(values)) range(values) else c(0,1)
    if(diff(limits)==0) limits<-limits+c(-.001,.001)
    save_plot(ggplot2::ggplot(z,ggplot2::aes(predictor,species,fill=.data[[metric]]))+
      ggplot2::geom_tile()+ggplot2::facet_grid(~coefficient_type,scales='free_x',space='free_x')+
      ggplot2::scale_fill_viridis_c(limits=limits,na.value='grey65')+
      ggplot2::theme_minimal(base_size=10)+ggplot2::theme(axis.text.y=ggplot2::element_text(size=9),axis.text.x=ggplot2::element_text(size=8,angle=45,hjust=1))+
      ggplot2::labs(title=paste('Beta',metric,'- species',min(match(groups[[i]],species)),'to',max(match(groups[[i]],species))),x=NULL,y=NULL),
      sprintf('print-heatmap-%s-%02d',metric,i),out,w=11,h=6.3)
  }
}
render_diagnostics <- function(out,template,title,notes) {
  directory<-normalizePath(out,mustWork=TRUE)
  beta_file<-file.path(directory,'beta-diagnostics.csv')
  if(file.exists(beta_file)) print_heatmaps(read.csv(beta_file),directory)
  parameters<-list(directory=directory,title=title,notes=notes)
  local_template<-file.path(directory,'report-source.Rmd')
  if(!file.copy(template,local_template,overwrite=TRUE)) stop('Could not stage report template.')
  for(fmt in c('html_document','pdf_document')) rmarkdown::render(local_template,
    output_format=fmt,output_file=if(fmt=='html_document') 'report.html' else 'report.pdf',
    output_dir=directory,intermediates_dir=directory,params=parameters,envir=new.env(parent=globalenv()),quiet=TRUE)
}

run_diagnostics <- function(args,helper_dir) {
  for(pkg in c('Hmsc','posterior','ggplot2','rmarkdown','knitr')) if(!requireNamespace(pkg,quietly=TRUE)) stop('Required package: ',pkg)
  if(!length(args)||!nzchar(args[1])) stop('Supply an explicit model-folder pattern: S3_check_model_fit_TensorFlow.R PATTERN [samples=250] [thin=100] [chains=4] [workers=4]')
  if(length(args)>5) stop('At most five positional arguments are supported.')
  settings<-c(250L,100L,4L,4L)
  if(length(args)>1) for(i in 2:length(args)) {v<-suppressWarnings(as.numeric(args[i]));if(length(v)!=1||!is.finite(v)||v<1||v!=floor(v)) stop('Sampling arguments must be positive integers.');settings[i-1]<-v}
  root<-Sys.getenv('CH1_MODEL_OUTPUT_DIR','HmscOutputs')
  folders<-list.dirs(root,recursive=FALSE,full.names=TRUE);folders<-sort(folders[grepl(args[1],basename(folders),fixed=TRUE)])
  if(!length(folders)) stop('No model folders matched under ',root)
  filename<-sprintf('HPC_samples_%04d_thin_%02d_chains_%d.Rdata',settings[1],settings[2],settings[3])
  inputs<-file.path(folders,'Models','Fitted',filename)
  if(any(!file.exists(inputs))) stop('Missing requested fitted inputs: ',paste(inputs[!file.exists(inputs)],collapse=', '))
  hashes<-tools::md5sum(inputs)
  site_sets<-lapply(inputs,function(f) {m<-read_model(f);as.character(m$studyDesign$site)})
  shared<-sort(Reduce(intersect,site_sets));all_beta<-list();overview<-list()
  template<-file.path(helper_dir,'report.Rmd');stamp<-sub('.Rdata$','',filename)
  for(idx in seq_along(inputs)) {
    input<-inputs[idx];modelname<-basename(folders[idx]);message('Diagnosing ',modelname)
    m<-read_model(input);b<-extract_beta(m)
    stopifnot(dim(b$draws)[1]==settings[1],dim(b$draws)[2]==settings[3],m$samples==settings[1],m$thin==settings[2])
    out<-file.path(folders[idx],'Results','model-diagnostics',stamp);dir.create(out,recursive=TRUE,showWarnings=FALSE)
    # A rerun replaces only PNGs produced by this workflow; stale selections must not enter reports.
    old_png<-list.files(out,pattern='\\.png$',full.names=TRUE);if(length(old_png)) unlink(old_png)
    d<-add_species(summarise_array(b$draws,b$key),m)
    ch<-do.call(rbind,lapply(seq_len(nrow(b$key)),function(j) cbind(b$key[j,],chain_summary(b$draws[,,j]))));write_table(ch,'beta-chains',out)
    half<-summarise_array(b$draws[seq.int(floor(dim(b$draws)[1]/2)+1,dim(b$draws)[1]),,,drop=FALSE],b$key)
    half$full_median<-d$median;half$full_rhat<-d$rhat;half$median_shift_in_full_sd<-(half$median-d$median)/d$sd
    write_table(half,'beta-second-half',out)
    d<-concentration(d,out);write_table(d,'beta-diagnostics',out);plot_beta(d,b$draws,out)
    write_table(support_summary(d),'directional-support',out)
    for(v in c('Gamma','rho')) {
      first<-m$postList[[1]][[1]][[v]]
      if(!is.null(first) && length(first) && !all(is.na(first))) {
        a<-array(NA_real_,c(settings[1],settings[3],length(first)))
        for(k in seq_len(settings[3])) for(i in seq_len(settings[1])) {z<-m$postList[[k]][[i]][[v]];stopifnot(length(z)==length(first));a[i,k,]<-as.vector(z)}
        key<-if(v=='Gamma') expand.grid(predictor=m$covNames,trait=m$trNames,stringsAsFactors=FALSE) else data.frame(parameter=seq_along(first))
        write_table(summarise_array(a,key),paste0(tolower(v),'-diagnostics'),out)
      }
    }
    notes<-spatial_diagnostics(m,b,shared,out)
    MF<-fit_metrics(m,input,settings[4]);ft<-data.frame(species=m$spNames)
    for(v in names(MF)) {if(length(MF[[v]])!=nrow(ft)) stop('Unexpected fit metric size: ',v);ft[[v]]<-as.numeric(MF[[v]])}
    sp<-read.csv(file.path(out,'concentration-species.csv'));ft$fraction_beta_flagged<-sp$fraction_flagged[match(ft$species,sp$group)]
    ft$prevalence<-d$prevalence[match(ft$species,d$species)];write_table(ft,'in-sample-fit',out)
    for(v in names(MF)) for(x in c('prevalence','fraction_beta_flagged')) save_plot(ggplot2::ggplot(ft,ggplot2::aes(.data[[x]],.data[[v]]))+ggplot2::geom_point()+ggplot2::theme_minimal()+ggplot2::labs(title=paste('In-sample',v)),paste0('fit-',make.names(v),'-vs-',x),out)
    provenance<-list(input=normalizePath(input),md5=unname(hashes[idx]),settings=setNames(settings,c('samples','thin','chains','workers')),thresholds=diag_thresholds,transient=m$transient,chain_draw_counts=lengths(m$postList),n_cells=nrow(m$Y),shared_cells=length(shared),n_species=m$ns,n_predictors=length(m$covNames),saved_scale='Hmsc saved Beta on X scale; projections use X, not XScaled',XScalePar=m$XScalePar,created=as.character(Sys.time()),versions=sapply(c('Hmsc','posterior','ggplot2','rmarkdown','knitr'),function(p) as.character(utils::packageVersion(p))))
    dput(provenance,file.path(out,'provenance.R'));capture.output(sessionInfo(),file=file.path(out,'session-info.txt'))
    render_diagnostics(out,template,modelname,notes)
    all_beta[[idx]]<-transform(d,model=modelname)
    overview[[idx]]<-data.frame(model=modelname,n_cells=nrow(m$Y),shared_cells=length(shared),n_beta=nrow(d),passing=sum(d$status=='pass'),flagged=sum(d$status=='flagged'),cannot_assess=sum(d$status=='cannot assess'),median_rhat=median(d$rhat,na.rm=TRUE),report=normalizePath(file.path(out,'report.html')),report_pdf=normalizePath(file.path(out,'report.pdf')))
    rm(m,b);gc(verbose=FALSE)
    stopifnot(identical(unname(tools::md5sum(input)),unname(hashes[idx])))
    message('Completed ',modelname)
  }
  comparison<-file.path(root,'model-diagnostics-comparison',paste0(make.names(args[1]),'-',stamp));dir.create(comparison,recursive=TRUE,showWarnings=FALSE)
  overview_df<-do.call(rbind,overview)
  write_table(overview_df,'model-comparison',comparison)
  comparison_plot<-do.call(rbind,lapply(c('passing','flagged','cannot_assess'),function(v) data.frame(model=overview_df$model,status=v,fraction=overview_df[[v]]/overview_df$n_beta)))
  save_plot(ggplot2::ggplot(comparison_plot,ggplot2::aes(model,fraction,fill=status))+ggplot2::geom_col()+ggplot2::coord_flip()+ggplot2::theme_minimal()+ggplot2::labs(title='Beta screening across selected models',x=NULL),'beta-comparison',comparison)
  shared_tables<-lapply(seq_along(folders),function(i) {
    od<-file.path(folders[i],'Results','model-diagnostics',stamp)
    files<-list.files(od,pattern='shared-cells-combined-correlations\\.csv$',full.names=TRUE)
    if(!length(files)) return(NULL)
    do.call(rbind,lapply(files,function(f) transform(read.csv(f),model=basename(folders[i]),random_level=sub('-shared-cells.*','',basename(f)))))
  })
  shared_df<-do.call(rbind,shared_tables)
  if(!is.null(shared_df)) {
    write_table(shared_df,'shared-cell-spatial-comparison',comparison)
    save_plot(ggplot2::ggplot(shared_df,ggplot2::aes(model,median,colour=model))+ggplot2::geom_boxplot(outlier.size=.4)+ggplot2::facet_grid(random_level~predictor)+ggplot2::theme_minimal(base_size=8)+ggplot2::theme(axis.text.x=ggplot2::element_blank(),legend.position='bottom')+ggplot2::labs(title=paste('Combined spatial overlap on',length(shared),'shared cells'),y='Species posterior median correlation',x=NULL),'shared-cell-spatial-comparison',comparison,w=16,h=7)
  }
  combined<-do.call(rbind,all_beta);write_table(support_summary(combined),'directional-support-combined',comparison)
  render_diagnostics(comparison,template,'Selected model comparison',paste('Primary family: all environmental slopes across selected models. Shared spatial comparison:',length(shared),'cells. See each linked report for model-specific diagnostics.'))
  stopifnot(identical(unname(tools::md5sum(inputs)),unname(hashes)))
  message('Comparison: ',normalizePath(file.path(comparison,'report.html')))
  invisible(comparison)
}
