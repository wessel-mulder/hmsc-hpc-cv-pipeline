# Run from repository root. Correlate original environmental predictors with
# species-specific spatial effects Eta %*% Lambda, on 1,697 shared cells.
# No ordinary correlation tests/p-values; cell independence is not assumed.
required <- c("ggplot2", "posterior")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Install: ",paste(missing,collapse=", "))
suppressPackageStartupMessages(library(ggplot2))
model_root <- "data-new/2-hmsc/derived/models"
model_stamp <- "2026-11-02_06-55-21"
output_dir <- "scripts/5-vis/spatial-latents"
variable_labels <- c(tmean_breeding="Temperature",prec_breeding="Precipitation",
  hh="Land-use heterogeneity",perc_urban="Urban cover",perc_cropland="Cropland cover",
  perc_pasture="Pasture cover",perc_forest_deciduous="Deciduous forest cover",
  perc_forest_evergreen="Evergreen forest cover",perc_forest_mixed_other="Mixed/other forest cover",
  perc_grass_shrub="Grass/shrubland cover")
atlas_labels <- c("1970s","1990s","2010s")
atlas_colours <- c("1970s"="#0072B2","1990s"="#D55E00","2010s"="#009E73")
read_fit <- function(atlas) {
  folder <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage",model_stamp,atlas)
  e <- new.env(parent=emptyenv())
  load(file.path(model_root,folder,"Models","Fitted","HPC_samples_0250_thin_100_chains_4.Rdata"),e)
  e$fitted_model$posteriors
}
metadata <- lapply(1:3,function(a) {
  m <- read_fit(a)
  stopifnot(length(m$postList)==4L,all(lengths(m$postList)==250L),m$ns==155L,m$nr==1L,
    setequal(colnames(m$XData),names(variable_labels)),
    identical(rownames(m$XData),rownames(m$studyDesign)))
  list(sites=as.character(m$studyDesign$site),species=m$spNames)
})
common_sites <- sort(Reduce(intersect,lapply(metadata,`[[`,"sites")))
species <- metadata[[1]]$species
stopifnot(length(common_sites)==1697L,!anyDuplicated(species),
  all(vapply(metadata,function(x) identical(x$species,species),logical(1))))

effect_species_rows <- list(); effect_community_rows <- list(); effect_draws <- list()
diagnose <- function(z) c(rhat=posterior::rhat(z),ess_bulk=posterior::ess_bulk(z),
                          ess_tail=posterior::ess_tail(z))
for(atlas in 1:3) {
  m <- read_fit(atlas)
  sites <- as.character(m$studyDesign$site)
  ix <- match(common_sites,sites)
  stopifnot(!anyNA(ix),!anyDuplicated(sites),
    identical(as.character(m$dfPi$site),sites),
    identical(as.integer(m$Pi[,1]),match(sites,levels(m$dfPi$site))))
  x <- as.matrix(m$XData[ix,names(variable_labels),drop=FALSE])
  stopifnot(all(is.finite(x)),all(apply(x,2,sd)>0))
  x <- scale(x) # Positive rescaling leaves Pearson correlations unchanged.
  n <- nrow(x)
  correlations <- array(NA_real_,c(250L,4L,10L,155L),
    dimnames=list(NULL,NULL,names(variable_labels),species))
  for(chain in 1:4) for(iteration in 1:250) {
    draw <- m$postList[[chain]][[iteration]]
    eta <- draw$Eta[[1]][m$Pi[ix,1],,drop=FALSE]
    lambda <- draw$Lambda[[1]]
    stopifnot(identical(dim(eta),c(n,5L)),identical(dim(lambda),c(5L,155L)),
      all(is.finite(eta)),all(is.finite(lambda)))
    centred_eta <- scale(eta,center=TRUE,scale=FALSE)
    # Cov(X, Eta Lambda) and Var(Eta Lambda), avoiding a full site x species
    # matrix for every draw. X has sample SD=1 on the evaluation cells.
    covariance <- crossprod(x,centred_eta) %*% lambda / (n-1)
    effect_variance <- colSums(lambda * (crossprod(centred_eta) %*% lambda)) / (n-1)
    if(any(!is.finite(effect_variance)) || any(effect_variance<=0)) stop("Zero/invalid spatial variance.")
    r <- sweep(covariance,2,sqrt(effect_variance),"/")
    if(any(!is.finite(r)) || any(abs(r)>1+1e-9)) stop("Invalid correlations.")
    correlations[iteration,chain,,] <- pmax(-1,pmin(1,r))
    if(iteration %in% c(1L,250L) && chain %in% c(1L,4L)) {
      direct <- cor(x,eta %*% lambda)
      stopifnot(max(abs(r-direct))<1e-10)
      # An orthogonal factor rotation leaves the combined spatial effect unchanged.
      q <- diag(5)[,c(5,3,1,4,2)]; q[,1] <- -q[,1]
      stopifnot(max(abs(cor(x,(eta %*% q) %*% (t(q) %*% lambda))-r))<1e-10)
    }
  }
  effect_draws[[atlas]] <- correlations
  effect_species_rows[[atlas]] <- do.call(rbind,lapply(1:10,function(v) {
    do.call(rbind,lapply(1:155,function(s) {
      z <- correlations[,,v,s]
      q <- quantile(z,c(.025,.5,.975),names=FALSE)
      d <- diagnose(z)
      data.frame(atlas=atlas,variable=names(variable_labels)[v],species=species[s],
        median_r=q[2],lower_r=q[1],upper_r=q[3],rhat=d["rhat"],
        ess_bulk=d["ess_bulk"],ess_tail=d["ess_tail"],row.names=NULL)
    }))
  }))
  effect_community_rows[[atlas]] <- do.call(rbind,lapply(1:10,function(v) {
    # Summarise all species jointly inside each draw; no species resampling.
    median_abs <- apply(abs(correlations[,,v,]),c(1,2),median)
    mean_abs <- apply(abs(correlations[,,v,]),c(1,2),mean)
    do.call(rbind,lapply(c("median_abs","mean_abs"),function(metric) {
      z <- if(metric=="median_abs") median_abs else mean_abs
      q <- quantile(z,c(.025,.5,.975),names=FALSE); d <- diagnose(z)
      data.frame(atlas=atlas,variable=names(variable_labels)[v],metric=metric,
        estimate=q[2],lower=q[1],upper=q[3],rhat=d["rhat"],
        ess_bulk=d["ess_bulk"],ess_tail=d["ess_tail"],row.names=NULL)
    }))
  }))
  message("Atlas ",atlas,": all ten predictors x 155 species x 1,000 draws checked.")
  rm(m); gc(verbose=FALSE)
}
effect_species_summary <- do.call(rbind,effect_species_rows)
effect_community_summary <- do.call(rbind,effect_community_rows)
stopifnot(nrow(effect_species_summary)==4650L,nrow(effect_community_summary)==60L)
stopifnot(all(is.finite(as.matrix(effect_species_summary[,c("rhat","ess_bulk","ess_tail")]))),
  all(is.finite(as.matrix(effect_community_summary[,c("rhat","ess_bulk","ess_tail")]))))
effect_diagnostic_summary <- do.call(rbind,lapply(1:3,function(a) do.call(rbind,lapply(
  names(variable_labels),function(v) {
    d <- subset(effect_species_summary,atlas==a & variable==v)
    c <- subset(effect_community_summary,atlas==a & variable==v & metric=="median_abs")
    data.frame(atlas=a,variable=v,n_species=nrow(d),
      median_species_rhat=median(d$rhat),n_species_rhat_ge_1_1=sum(d$rhat>=1.1),
      median_abs_rhat=c$rhat,median_abs_ess_bulk=c$ess_bulk)
  }))))

dat <- effect_species_summary
dat$atlas <- factor(atlas_labels[dat$atlas],levels=atlas_labels)
dat$variable <- factor(dat$variable,levels=names(variable_labels),labels=unname(variable_labels))
distribution_plot <- ggplot(dat,aes(atlas,median_r,colour=atlas,fill=atlas)) +
  geom_hline(yintercept=0,linetype="dashed",colour="grey65") +
  geom_violin(trim=TRUE,scale="width",alpha=.18,linewidth=.5) +
  geom_point(position=position_jitter(width=.15,height=0,seed=20261106),size=.7,alpha=.45) +
  geom_boxplot(width=.12,outlier.shape=NA,fill="white",linewidth=.45) +
  facet_wrap(~variable,ncol=2) + scale_colour_manual(values=atlas_colours) +
  scale_fill_manual(values=atlas_colours) + scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,.5)) +
  labs(title="Environmental overlap with the combined spatial random effect",
    subtitle="Each point: one species' median draw-level correlation with Eta × Lambda; same 1,697 cells",
    x=NULL,y="Pearson correlation across cells",
    caption="155 matched species per atlas; all 1,000 draws used. Original predictors, without spatial detrending.\nViolin/boxplot describe species differences, not posterior uncertainty. Spatial effects are inferred, not independently measured.\nOverlap alone is not a measure of coefficient bias or proof of spatial confounding; inspect chain diagnostics.") +
  theme_minimal(base_size=12) + theme(legend.position="none",panel.grid.minor=element_blank(),
    strip.text=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(output_dir,"environment-spatial-effect-correlations.png"),distribution_plot,
  width=12,height=16,dpi=300,bg="white")

d <- subset(effect_community_summary,metric=="median_abs")
d$atlas <- factor(atlas_labels[d$atlas],levels=atlas_labels)
d$variable <- factor(d$variable,levels=rev(names(variable_labels)),labels=rev(unname(variable_labels)))
pd <- position_dodge(width=.65)
strength_plot <- ggplot(d,aes(estimate,variable,colour=atlas,group=atlas)) +
  geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",position=pd,width=.2,linewidth=.6) +
  geom_point(position=pd,size=2.2) + scale_colour_manual(values=atlas_colours) +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.1)) +
  labs(title="How pervasive is environmental overlap with space?",
    subtitle="Median absolute correlation across 155 species, calculated jointly within each posterior draw",
    x="Median across species of |cor(predictor, Eta × Lambda)|",y=NULL,colour="Atlas",
    caption="Points: pooled-draw medians; bars: 2.5–97.5% draw intervals. Same 1,697 cells in all atlases.\nAbsolute correlations prevent positive and negative species responses from cancelling.\nIntervals are conditional model summaries; poor chain agreement limits inference. No significance tests.") +
  theme_minimal(base_size=12) + theme(legend.position="top",panel.grid.minor=element_blank(),
    plot.caption=element_text(hjust=0,size=10))
ggsave(file.path(output_dir,"environment-spatial-effect-correlation-strength.png"),strength_plot,
  width=12,height=9,dpi=300,bg="white")
message("Saved two environment–spatial-effect PNGs; no model changes.")
print(subset(effect_community_summary,metric=="median_abs")[,
  c("atlas","variable","estimate","lower","upper","rhat")],row.names=FALSE)
