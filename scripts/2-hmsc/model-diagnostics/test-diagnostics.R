# Rscript scripts/2-hmsc/model-diagnostics/test-diagnostics.R
source('scripts/2-hmsc/model-diagnostics/helpers.R')
set.seed(369)
z<-matrix(rnorm(8000),2000,4)
d<-diagnose_matrix(z)
stopifnot(d$status=='pass',identical(d$rhat,as.numeric(posterior::rhat(z))),identical(d$ess_bulk,as.numeric(posterior::ess_bulk(z))),identical(d$ess_tail,as.numeric(posterior::ess_tail(z))))
shift<-z;shift[,1]<-shift[,1]+3
stopifnot(diagnose_matrix(shift)$flag_rhat)
ar<-replicate(4,as.numeric(arima.sim(list(ar=.99),n=2000)))
stopifnot(diagnose_matrix(ar)$flag_bulk)
stopifnot(diagnose_matrix(matrix(1,50,4))$status=='cannot assess')
invalid<-z;invalid[1,1]<-NA
stopifnot(diagnose_matrix(invalid)$status=='cannot assess')
x<-cbind('(Intercept)'=1,temperature=rnorm(100),constant=2)
eta<-matrix(rnorm(300),100,3);lambda<-matrix(rnorm(15),3,5)
r<-spatial_cor(x[,2,drop=FALSE],eta,lambda)
stopifnot(max(abs(r-cor(x[,2,drop=FALSE],eta%*%lambda)))<1e-12)
q<-qr.Q(qr(matrix(rnorm(9),3)))
stopifnot(max(abs(r-spatial_cor(x[,2,drop=FALSE],eta%*%q,t(q)%*%lambda)))<1e-12,qr(x)$rank<ncol(x))
m<-list(covNames=c('(Intercept)','temperature'),spNames=c('a','b'),XScaled=x[,1:2],X=x[,1:2],Y=matrix(c(0,1),100,2,dimnames=list(NULL,c('a','b'))),ranLevels=list(),TrData=NULL)
m$postList<-lapply(1:4,function(k) lapply(1:50,function(i) list(Beta=matrix(rnorm(4),2,2))))
b<-extract_beta(m)
stopifnot(identical(as.numeric(b$draws[3,2,]),as.vector(m$postList[[2]][[3]]$Beta)),b$key$species[3]=='b',b$key$predictor[3]=='(Intercept)')
d<-add_species(summarise_array(b$draws,b$key),m)
stopifnot(all(is.na(d$Migration_a3_DOF)),grepl('Inapplicable',spatial_diagnostics(m,b,character(),tempdir())))
# Support zero draws must not be treated as perfect directional certainty.
stopifnot(diagnose_matrix(matrix(0,50,4))$wrong_sign==1)
# Downstream cache readers still load MF; a second call validates reuse.
folder<-tempfile();dir.create(folder);input<-file.path(folder,'HPC_test.Rdata');writeLines('fingerprint',input)
MF<-list(AUC=c(.8,.9));MF_provenance<-list(input_md5=unname(tools::md5sum(input)),Hmsc=as.character(packageVersion('Hmsc')),mode='in-sample all saved draws expected=TRUE',version=1L)
save(MF,MF_provenance,file=file.path(folder,'MF_preCV_test.Rdata'))
stopifnot(identical(fit_metrics(NULL,input,1),MF))
cat('Simulation, direct posterior, indexing, spatial invariance, environment-only, missing trait, invalid diagnostic and cache checks passed.\n')
# Projection is on the exported X scale, and preserves the fitted total.
xx<-cbind(1,rnorm(100)); ee<-matrix(rnorm(300),100,3);ll<-matrix(rnorm(6),3,2)
proj<-qr.coef(qr(xx),ee)%*%ll
stopifnot(max(abs(proj-qr.coef(qr(xx),ee%*%ll)))<1e-12)
bb<-matrix(rnorm(4),2,2)
stopifnot(max(abs(xx%*%(bb+proj)-(xx%*%bb+xx%*%proj)))<1e-12)
cat('Projection and summed-component algebra checks passed.\n')
if('--integration' %in% commandArgs(trailingOnly=TRUE)) {
  fixture_root<-tempfile('s3-fixture-');folder<-file.path(fixture_root,'synthetic-environment-only','Models','Fitted');dir.create(folder,recursive=TRUE)
  set.seed(369);ns<-40L
  m<-list(covNames=c('(Intercept)','temperature','constant'),spNames=c('a','b'),
    X=cbind('(Intercept)'=1,temperature=rnorm(30),constant=2),
    Y=matrix(rep(c(0,1,1),20),30,2,dimnames=list(NULL,c('a','b'))),ranLevels=list(),TrData=NULL,
    studyDesign=data.frame(site=paste0('cell',1:30)),ns=2L,samples=ns,thin=1L,XScalePar=rbind(c(0,0,0),c(1,1,1)))
  m$XScaled<-m$X
  m$postList<-lapply(1:2,function(k) lapply(1:ns,function(i) {bb<-matrix(rnorm(6),3,2);bb[3,]<-1;list(Beta=bb)}))
  fitted_model<-list(posteriors=m)
  input<-file.path(folder,'HPC_samples_0040_thin_01_chains_2.Rdata');save(fitted_model,file=input)
  MF<-list(AUC=c(.8,.9),TjurR2=c(.2,.3))
  MF_provenance<-list(input_md5=unname(tools::md5sum(input)),Hmsc=as.character(packageVersion('Hmsc')),mode='in-sample all saved draws expected=TRUE',version=1L)
  save(MF,MF_provenance,file=file.path(folder,'MF_preCV_samples_0040_thin_01_chains_2.Rdata'))
  original_root<-Sys.getenv('CH1_MODEL_OUTPUT_DIR',unset=NA_character_);Sys.setenv(CH1_MODEL_OUTPUT_DIR=fixture_root)
  status<-system2(file.path(R.home('bin'),'Rscript'),c('scripts/2-hmsc/S3a_show_model_fit.R','synthetic','40','1','2','1'),stdout=file.path(fixture_root,'run.log'),stderr=file.path(fixture_root,'run.log'))
  if(is.na(original_root)) Sys.unsetenv('CH1_MODEL_OUTPUT_DIR') else Sys.setenv(CH1_MODEL_OUTPUT_DIR=original_root)
  if(status!=0) stop(paste(readLines(file.path(fixture_root,'run.log')),collapse='\n'))
  output<-file.path(dirname(dirname(folder)),'Results','model-diagnostics','HPC_samples_0040_thin_01_chains_2')
  stopifnot(file.exists(file.path(output,'report.pdf')),file.exists(file.path(output,'report.html')))
  dd<-read.csv(file.path(output,'beta-diagnostics.csv'))
  stopifnot(all(dd$status[dd$predictor=='constant']=='cannot assess'),all(is.na(dd$Migration_a3_DOF)))
  html<-paste(readLines(file.path(output,'report.html'),warn=FALSE),collapse='\n')
  stopifnot(grepl('Inapplicable: no spatial random levels',html),grepl('data:image/png',html,fixed=TRUE))
  cat('Full legacy wrapper integration (environment-only, missing traits, constant predictor/draws, cached MF, HTML and PDF) passed.\n')
}
if('--models' %in% commandArgs(trailingOnly=TRUE)) {
  root<-Sys.getenv('CH1_MODEL_OUTPUT_DIR','HmscOutputs')
  files<-Sys.glob(file.path(root,'*','Results','model-diagnostics','*','provenance.R'))
  if(!length(files)) stop('No completed diagnostics with provenance under ',root)
  for(f in files) {
    out<-dirname(f);p<-dget(f);m<-read_model(p$input);d<-read.csv(file.path(out,'beta-diagnostics.csv'))
    stopifnot(unname(tools::md5sum(p$input))==p$md5,nrow(d)==m$ns*length(m$covNames),
      identical(d$species,rep(m$spNames,each=length(m$covNames))),identical(d$predictor,rep(m$covNames,m$ns)))
    for(j in unique(c(1L,ceiling(nrow(d)/2),nrow(d)))) {
      pi<-match(d$predictor[j],m$covNames);si<-match(d$species[j],m$spNames)
      z<-vapply(m$postList,function(ch) vapply(ch,function(dr) dr$Beta[pi,si],numeric(1)),numeric(p$settings['samples']))
      direct<-diagnose_matrix(z)
      for(v in c('median','sd','lower','upper','rhat','ess_bulk','ess_tail','mcse_mean','mcse_lower','mcse_upper')) stopifnot(isTRUE(all.equal(d[[v]][j],direct[[v]],tolerance=1e-10)))
    }
    chains<-read.csv(file.path(out,'beta-chains.csv'));stopifnot(nrow(chains)==nrow(d)*length(m$postList))
    e<-new.env();load(file.path(dirname(p$input),sub('HPC_','MF_preCV_',basename(p$input),fixed=TRUE)),e)
    stopifnot(is.list(e$MF),identical(e$MF_provenance$input_md5,p$md5))
    selected<-read.csv(file.path(out,'selected-coefficients.csv'))
    bad<-unique(c(head(order(d$rhat,decreasing=TRUE,na.last=NA),5),head(order(d$ess_bulk,na.last=NA),5),head(order(d$ess_tail,na.last=NA),5)))
    good<-head(which(d$status=='pass')[order(d$species[d$status=='pass'],d$predictor[d$status=='pass'])],3)
    stopifnot(setequal(selected$parameter,unique(c(bad,good))))
    html<-paste(readLines(file.path(out,'report.html'),warn=FALSE),collapse='\n')
    stopifnot(grepl('data:image/png',html,fixed=TRUE),grepl('<table',html,fixed=TRUE),file.exists(file.path(out,'report.pdf')))
    for(sf in list.files(out,pattern='shared-cells-combined-correlations\\.csv$',full.names=TRUE)) stopifnot(all(read.csv(sf)$n_cells==p$shared_cells))
    message('Saved model extraction, direct posterior diagnostics, unchanged fingerprint, MF compatibility, deterministic selection, HTML/PDF and shared-cell checks passed: ',basename(dirname(dirname(dirname(out)))))
  }
}
