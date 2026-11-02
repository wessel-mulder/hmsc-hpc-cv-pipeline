# Rebuild the current workflow using scripts/ and data-new/ raw inputs only.
# --verify rebuilds into an empty temporary folder and compares saved outputs.
# --offline requires all 132 monthly CHELSAcruts crops already in data-new/.
args <- commandArgs(trailingOnly=TRUE)
verify <- '--verify' %in% args
for(p in c('ape','readxl','terra','sf','exactextractr','ggplot2','cowplot','dplyr'))
  if(!requireNamespace(p,quietly=TRUE))stop('Missing R package: ',p)
root <- if(verify)tempfile('ch1-reproduction-',tmpdir='/tmp') else 'data-new'
pre <- file.path(root,'0-pre/derived');hmsc <- file.path(root,'2-hmsc/derived')
dir.create(pre,recursive=TRUE,showWarnings=FALSE)
report_dir <- 'data-new/reproducibility'
dir.create(report_dir,recursive=TRUE,showWarnings=FALSE)
report <- c('CURRENT PREPROCESSING REPRODUCTION',paste('R version:',R.version.string),paste('Fresh staging directory:',root),
            'Inputs: raw files in data-new/0-pre/raw. No prior derived environmental/trait/tree tables used for rebuilding.',
            'Study design and occurrence CSVs are supplied raw starting inputs; their upstream compilation is outside this workflow.')
record <- function(x){message(x);report <<- c(report,x)}
status <- system2(file.path(R.home('bin'),'Rscript'),'scripts/0-pre/tests/check-environment.R', stdout=file.path(report_dir,'unit-checks.log'),stderr=file.path(report_dir,'unit-checks.log'))
if(status!=0)stop('Scientific unit checks failed; see unit-checks.log')
record('PASS scientific unit checks')
source('scripts/0-pre/preparing-traits.R')
source('scripts/0-pre/preparing-taxonomy.R')
record('1/5: rebuilding species-level traits and explicit Corvus cornix STI proxy')
prepare_trait_table(output_dir=file.path(pre,'Tr'))
record('2/5: rebuilding IOC rank tree and PNG')
prepare_taxonomy(taxa_file=file.path(pre,'Tr/Trait_table.csv'),output_dir=file.path(pre,'taxonomy'))
Sys.setenv(CH1_PREPROCESSING_DERIVED=pre,CH1_HMSC_DERIVED=hmsc)
record('3/5: rebuilding environmental data from original HILDA states and monthly CHELSAcruts crops')
env_log <- file.path(report_dir,'environment-rebuild.log')
options <- c('scripts/0-pre/preparing-environment.R','--refresh-landuse','--workers=3')
if('--offline' %in% args)options<-c(options,'--offline')
status<-system2(file.path(R.home('bin'),'Rscript'),options,stdout=env_log,stderr=env_log)
if(status!=0)stop('Environmental rebuild failed; see ',env_log)
record('4/5: rebuilding and independently validating S0 RData')
sys.source('scripts/2-hmsc/S0_Data_Definitions.R',envir=new.env(parent=globalenv()))
file_name<-'preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
z<-new.env();load(file.path(hmsc,file_name),envir=z)
x<-read.csv(file.path(pre,'environmental/X_environment.csv'))
t<-read.csv(file.path(pre,'Tr/Trait_table.csv'));tree<-ape::read.tree(file.path(pre,'taxonomy/tree_fromPD.tre'))
stopifnot(nrow(x)==6495,ncol(x)==21,nrow(t)==226,length(tree$tip.label)==226,
          setequal(t$species,tree$tip.label),!anyNA(z$X),!anyNA(z$Tr),ncol(z$X)==10,
          identical(rownames(z$X),rownames(z$Y)),identical(rownames(z$Y),rownames(z$Design)),
          identical(colnames(z$Y),rownames(z$Tr)),all(colnames(z$Y)%in%tree$tip.label))
j<-match(rownames(z$X),x$survey)
stopifnot(all(x$pct_class_0[j]+x$pct_class_77[j]<=75),
          identical(z$X$hh,x$hh[j]),
          identical(z$X$perc_cropland,x$pct_class_22[j]),
          identical(z$X$perc_forest_deciduous,x$pct_class_43[j]+x$pct_class_44[j]),
          identical(z$X$perc_forest_evergreen,x$pct_class_42[j]+x$pct_class_41[j]),
          identical(z$X$perc_forest_mixed_other,x$pct_class_40[j]+x$pct_class_45[j]))
for(a in 1:3)stopifnot(all(colSums(z$Y[z$Design$atlas==a,,drop=FALSE])>=5))
# Compare regenerated objects numerically; binary RData compression bytes may differ.
if(verify){
  baseline<-new.env();load(file.path('data-new/2-hmsc/derived',file_name),envir=baseline)
  for(name in c('X','Y','Tr','Design')){
    eq<-all.equal(z[[name]],baseline[[name]],tolerance=1e-7)
    if(!isTRUE(eq))stop('Saved ',name,' differs from regeneration: ',paste(eq,collapse='; '))
    record(paste('PASS regenerated',name,'matches saved RData'))
  }
  for(path in c('environmental/X_environment.csv','Tr/Trait_table.csv')){
    a<-read.csv(file.path(pre,path));b<-read.csv(file.path('data-new/0-pre/derived',path))
    if(!isTRUE(all.equal(a,b,tolerance=1e-7)))stop('CSV differs from saved baseline: ',path)
    record(paste('PASS regenerated CSV matches:',path))
  }
  old_tree<-ape::read.tree('data-new/0-pre/derived/taxonomy/tree_fromPD.tre')
  stopifnot(isTRUE(all.equal(ape::cophenetic.phylo(tree)[t$species,t$species],
                            ape::cophenetic.phylo(old_tree)[t$species,t$species],tolerance=1e-8)))
  record('PASS regenerated tree has identical tips and patristic distances')
}
record('5/5: regenerating current environmental maps and breeding climate point cloud')
source('scripts/0-pre/plot-environment.R')
plot_environment(file.path(pre,'environmental/X_environment.csv'),file.path(pre,'environmental/figures/current'))
record(paste('PASS:',nrow(z$X),'surveys;',ncol(z$Y),'species;',ncol(z$X),'predictors;',length(tree$tip.label),'tree tips'))
record('S1 was not executed: model definitions still require new paths/predictors and scaling review; fitting has not been reproduced.')
record('Historic expanded-variable plots are superseded by figures/current; no claim of reproducing earlier exploratory schemas.')
record(paste('Outputs verified at:',normalizePath(root)))
report_name <- if(verify)'reproduction-report.txt' else 'rebuild-report.txt'
writeLines(report,file.path(report_dir,report_name))
writeLines(capture.output(sessionInfo()),file.path(report_dir,'sessionInfo.txt'))
# Input/source fingerprints for this run, excluding generated outputs and unused transition rasters.
source_files<-c(list.files('scripts/0-pre',recursive=TRUE,full.names=TRUE),
                'scripts/2-hmsc/S0_Data_Definitions.R','scripts/reproduce-preprocessing.R')
input_files<-list.files('data-new/0-pre/raw',recursive=TRUE,full.names=TRUE)
input_files<-input_files[!grepl('/transitions/|\\.partial',input_files)]
manifest<-data.frame(path=c(source_files,input_files))
manifest$md5<-unname(tools::md5sum(manifest$path))
write.csv(manifest,file.path(report_dir,'source-input-checksums.csv'),row.names=FALSE)
# Retain verification staging for inspection; it is outside the supplied data-new baseline.
message('Report: ',file.path(report_dir,report_name))
