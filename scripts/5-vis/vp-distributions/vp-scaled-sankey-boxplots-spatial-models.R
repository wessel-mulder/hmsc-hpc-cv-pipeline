# Run from project root:
# Rscript scripts/5-vis/vp-distributions/vp-scaled-sankey-boxplots-spatial-models.R
# Exploratory paired signed-rank tests on exported species-level VP estimates.
# Not a posterior uncertainty analysis or an omnibus distribution-equality test.
model_root <- "data-new/2-hmsc/derived/models"
model_stamp <- "2026-11-02_06-55-21"
atlas_labels <- c("1970s", "1990s", "2010s")
variable_labels <- c(tmean_breeding="Temperature",prec_breeding="Precipitation",
  hh="Land-use heterogeneity",perc_urban="Urban cover",perc_cropland="Cropland cover",
  perc_pasture="Pasture cover",perc_forest_deciduous="Deciduous forest cover",
  perc_forest_evergreen="Evergreen forest cover",perc_forest_mixed_other="Mixed/other forest cover",
  perc_grass_shrub="Grass/shrubland cover")
reference_species <- NULL; vp_tables <- list(); total_rows <- list()
for(a in 1:3) {
  model <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, a)
  raw <- read.csv(file.path(model_root,model,"Results",paste0(model,"parameter_estimates_VP_.csv")),
                   row.names=1,check.names=FALSE)
  stopifnot(setequal(rownames(raw),c(names(variable_labels),"Random: site","TjurR2")),ncol(raw)==155L,
    !anyDuplicated(colnames(raw)),all(vapply(raw,is.numeric,logical(1))),all(is.finite(as.matrix(raw))))
  if(is.null(reference_species)) reference_species <- colnames(raw)
  stopifnot(setequal(colnames(raw),reference_species))
  raw <- raw[,reference_species,drop=FALSE]
  relative <- as.matrix(raw[names(variable_labels),,drop=FALSE])
  stopifnot(all(relative>=0),all(relative<=1),max(abs(colSums(relative)+as.numeric(raw["Random: site",])-1))<1e-8)
  tjur <- as.numeric(raw["TjurR2",])
  stopifnot(all(tjur>=0),all(tjur<=1))
  vp_tables[[a]] <- list(scaled=sweep(relative,2,tjur,"*")*100)
  spatial_scaled <- as.numeric(raw["Random: site",])*tjur*100
  stopifnot(max(abs(colSums(vp_tables[[a]]$scaled)+spatial_scaled-tjur*100))<1e-8)
  total_rows[[a]] <- data.frame(atlas=a,mean_tjur_percent=mean(tjur)*100,
    environmental_percent=sum(rowMeans(vp_tables[[a]]$scaled)),spatial_percent=mean(spatial_scaled))
}
pairs <- list(c(1L,2L),c(2L,3L),c(1L,3L))
rows <- list(); paired_rows <- list(); row_i <- 0L
for(scale_name in "scaled") for(pair in pairs) for(variable in names(variable_labels)) {
  row_i <- row_i+1L
  earlier <- vp_tables[[pair[1]]][[scale_name]][variable,reference_species]
  later <- vp_tables[[pair[2]]][[scale_name]][variable,reference_species]
  stopifnot(identical(names(earlier),names(later)),all(is.finite(earlier)),all(is.finite(later)))
  delta <- later-earlier
  nonzero <- delta[delta!=0]
  abs_rank <- rank(signif(abs(nonzero),7),ties.method="average")
  statistic <- sum(abs_rank[nonzero>0])
  rank_biserial <- if(length(nonzero)) (2*statistic-sum(abs_rank))/sum(abs_rank) else NA_real_
  if(length(nonzero)) {
    test <- wilcox.test(later,earlier,paired=TRUE,alternative="two.sided",mu=0,
      exact=FALSE,correct=TRUE,conf.int=TRUE,conf.level=.95,tol.root=1e-8,digits.rank=7)
    stopifnot(abs(unname(test$statistic)-statistic)<1e-10)
    # Independently verify equivalent one-sample-difference and shuffled-ID tests.
    direct <- wilcox.test(delta,mu=0,alternative="two.sided",exact=FALSE,correct=TRUE,digits.rank=7)
    reversed <- later[rev(seq_along(later))]
    rematched <- reversed[match(names(earlier),names(reversed))]
    shuffled <- wilcox.test(rematched,earlier,paired=TRUE,exact=FALSE,correct=TRUE,digits.rank=7)
    stopifnot(abs(test$p.value-direct$p.value)<1e-12,abs(test$p.value-shuffled$p.value)<1e-12)
    p_raw <- test$p.value; hl <- unname(test$estimate)
    lower <- test$conf.int[1]; upper <- test$conf.int[2]
    note <- "normal approximation; continuity correction; tie adjustment"
  } else {
    p_raw <- hl <- lower <- upper <- NA_real_
    note <- "All paired differences zero; signed-rank inference undefined"
  }
  scale_sd <- sd(delta)
  skewness <- if(scale_sd>0) mean((delta-mean(delta))^3)/scale_sd^3 else NA_real_
  pair_label <- paste0(atlas_labels[pair[1]]," to ",atlas_labels[pair[2]])
  rows[[row_i]] <- data.frame(scale=scale_name,pair=pair_label,variable=variable,
    n_species=length(delta),n_nonzero=length(nonzero),n_positive=sum(delta>0),
    n_negative=sum(delta<0),n_zero=sum(delta==0),
    n_tied_absolute_differences=length(nonzero)-length(unique(signif(abs(nonzero),7))),
    mean_change_pp=mean(delta),median_change_pp=median(delta),
    hl_pseudomedian_change_pp=hl,hl_ci_lower_pp=lower,hl_ci_upper_pp=upper,
    rank_biserial=rank_biserial,V=statistic,p_raw=p_raw,
    difference_skewness=skewness,method_note=note)
  paired_rows[[row_i]] <- data.frame(scale=scale_name,pair=pair_label,variable=variable,
    species=reference_species,earlier_vp_percent=as.numeric(earlier),
    later_vp_percent=as.numeric(later),difference_pp=as.numeric(delta))
}
wilcoxon_results <- do.call(rbind,rows)
paired_species_values <- do.call(rbind,paired_rows)
stopifnot(nrow(wilcoxon_results)==30L,nrow(paired_species_values)==30L*155L,
          all(wilcoxon_results$n_positive+wilcoxon_results$n_negative+wilcoxon_results$n_zero==155L))
# Scaled-only family requested: ten predictors x three atlas pairs = 30 tests.
# Holm controls familywise error under arbitrary test dependence if the marginal
# p-values are valid. It cannot repair dependence among species within a test.
wilcoxon_results$p_holm <- p.adjust(wilcoxon_results$p_raw,method="holm",n=30L)
stopifnot(all(wilcoxon_results$p_holm>=wilcoxon_results$p_raw,na.rm=TRUE))
wilcoxon_summary <- do.call(rbind,lapply("scaled",function(scale_name) {
  do.call(rbind,lapply(vapply(pairs,function(p)paste0(atlas_labels[p[1]]," to ",atlas_labels[p[2]]),character(1)),
    function(pair_label) {
      d <- wilcoxon_results[wilcoxon_results$scale==scale_name & wilcoxon_results$pair==pair_label,]
      data.frame(scale=scale_name,pair=pair_label,n_tests=nrow(d),
        raw_p_below_005=sum(d$p_raw<.05,na.rm=TRUE),holm_p_below_005=sum(d$p_holm<.05,na.rm=TRUE))
  }))
}))

suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr)})
output_dir <- "scripts/5-vis/vp-distributions"
variable_colours <- c("Temperature"="firebrick3","Precipitation"="dodgerblue3",
  "Land-use heterogeneity"="orchid3","Urban cover"="snow3","Cropland cover"="goldenrod1",
  "Pasture cover"="darkorange","Deciduous forest cover"="#A6D96A",
  "Evergreen forest cover"="#006400","Mixed/other forest cover"="#87CEEB",
  "Grass/shrubland cover"="#B39DDB")

species_vp <- do.call(rbind,lapply(1:3,function(a)data.frame(
  atlas=atlas_labels[a],variable=rep(names(variable_labels),times=length(reference_species)),
  species=rep(reference_species,each=length(variable_labels)),scale="scaled",
  vp_percent=as.vector(vp_tables[[a]]$scaled))))
vp_total_check <- do.call(rbind,total_rows)
stopifnot(max(abs(vp_total_check$environmental_percent+vp_total_check$spatial_percent-
  vp_total_check$mean_tjur_percent))<1e-8)
dat <- species_vp; scale_name <- "scaled"
x_max <- ceiling(max(dat$vp_percent)/5)*5
  box_data <- dat
  box_data$variable <- factor(box_data$variable,levels=names(variable_labels),labels=unname(variable_labels))
  box_data$atlas <- factor(box_data$atlas,levels=atlas_labels)
  atlas_offset <- if(scale_name=="scaled") 1.3 else 1.1
  box_data$x_position <- as.integer(box_data$variable)*4 +
    c(-atlas_offset,0,atlas_offset)[as.integer(box_data$atlas)]
  box_data$box_group <- interaction(box_data$variable,box_data$atlas,drop=TRUE)
  brackets <- wilcoxon_results
  brackets <- brackets[brackets$scale==scale_name & brackets$p_holm<.05,]
  brackets$variable <- factor(brackets$variable,levels=names(variable_labels),labels=unname(variable_labels))
  pair_levels <- c("1970s to 1990s","1990s to 2010s","1970s to 2010s")
  bracket_level <- match(brackets$pair,pair_levels)
  centres <- as.integer(brackets$variable)*4
  brackets$x_start <- centres+c(-atlas_offset,0,-atlas_offset)[bracket_level]
  brackets$x_end <- centres+c(0,atlas_offset,atlas_offset)[bracket_level]
  panel_max <- tapply(box_data$vp_percent,box_data$variable,max)
  # Stack only the retained comparisons; the widest comparison goes last.
  retained_level <- ave(bracket_level,brackets$variable,FUN=function(z)rank(z,ties.method="first"))
  brackets$y <- unname(panel_max[as.character(brackets$variable)]) +
    x_max*(.045+.065*(retained_level-1))
  brackets$label <- vapply(brackets$p_holm,function(p) {
    if(p<.001) "p<0.001" else paste0("p=",formatC(p,format="f",digits=3))
  },character(1))
  if(scale_name=="scaled") brackets$label <- sub("0\\.",".",brackets$label)
  brackets$x_mid <- (brackets$x_start+brackets$x_end)/2
  label_half_gap <- if(scale_name=="scaled") .58 else .43
  brackets$left_end <- brackets$x_mid-label_half_gap
  brackets$right_start <- brackets$x_mid+label_half_gap
  # Guarantee visible segments even for adjacent-atlas comparisons.
  brackets$left_start <- pmin(brackets$x_start,brackets$left_end-.22)
  brackets$right_end <- pmax(brackets$x_end,brackets$right_start+.22)
  stopifnot(nrow(brackets)==sum(wilcoxon_results$scale==scale_name &
    wilcoxon_results$p_holm<.05),all(brackets$p_holm<.05),!anyNA(brackets$y))
  box_plot <- ggplot(box_data,aes(x=x_position,y=vp_percent,fill=variable,alpha=atlas,group=box_group)) +
    geom_boxplot(orientation="x",width=.65,colour="grey35",linewidth=.4,
      outlier.size=.7,outlier.alpha=.35) +
    stat_summary(fun=mean,geom="point",orientation="x",shape=23,size=1.6,
      fill="white",colour="grey25",alpha=1) +
    geom_segment(data=brackets,aes(x=left_start,xend=left_end,y=y,yend=y),
      inherit.aes=FALSE,linewidth=.35,colour="grey25") +
    geom_segment(data=brackets,aes(x=right_start,xend=right_end,y=y,yend=y),
      inherit.aes=FALSE,linewidth=.35,colour="grey25") +
    geom_text(data=brackets,aes(x=x_mid,y=y,label=label),
      inherit.aes=FALSE,size=if(scale_name=="scaled") 3.5 else 2.5,
      colour="grey15") +
    scale_fill_manual(values=variable_colours,guide="none") +
    scale_alpha_manual(name="Atlas (left to right)",values=c("1970s"=.45,"1990s"=.7,"2010s"=1)) +
    scale_x_continuous(breaks=seq_along(variable_labels)*4,
      labels=vapply(unname(variable_labels),function(z)paste(strwrap(z,width=17),collapse="\n"),character(1)),
      expand=expansion(add=1.8)) +
    scale_y_continuous(limits=c(0,max(c(dat$vp_percent,brackets$y))+x_max*.045),
      breaks=pretty(c(0,x_max),n=6),
      labels=function(x)paste0(format(x,trim=TRUE,scientific=FALSE),"%"),
      expand=expansion(mult=c(.005,0))) +
    labs(title=if(scale_name=="scaled") "Spatial models: Tjur R²-scaled variance partitions"
               else "Spatial models: relative variance partitions",
      subtitle="155 matched species per box; only Holm-adjusted paired Wilcoxon p < 0.05 shown (30-test family)",
      x=NULL,y=if(scale_name=="scaled") "Tjur R²-scaled VP (%)" else "Relative VP (%)",
      caption=paste0("November 2 spatial models; within each variable: 1970s, 1990s, 2010s from left to right. Lines connect compared atlas boxes.\n",
        "Box: median and IQR; whiskers: 1.5 × IQR; dots: outliers; white diamonds: means. All species estimates retained.\n",
        if(scale_name=="scaled") "Relative VP × species-specific Tjur R². " else "Relative environmental partitions sum to 100% within species. ",
        "Exploratory tests: posterior uncertainty and dependence among species are not accounted for.")) +
    theme_minimal(base_size=12) +
    theme(legend.position="top",legend.justification="left",
      panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
      axis.text.x=element_text(size=10),plot.title=element_text(face="bold"),
      plot.caption=element_text(hjust=0,size=9),plot.margin=margin(12,12,12,12))
  boxes <- ggplot_build(box_plot)$data[[1]]
  expected_medians <- as.numeric(tapply(box_data$vp_percent,list(box_data$variable,box_data$atlas),median))
  stopifnot(nrow(boxes)==30L,max(abs(sort(boxes$middle)-sort(expected_medians)))<1e-10)

means <- aggregate(vp_percent~atlas+variable,dat,mean)
means$atlas <- factor(means$atlas,levels=atlas_labels)
means$variable <- factor(means$variable,levels=names(variable_labels),
  labels=unname(variable_labels))
totals <- aggregate(vp_percent~atlas,means,sum)
stopifnot(nrow(means)==30L,max(abs(totals$vp_percent-
  vp_total_check$environmental_percent))<1e-10)
panel_a <- ggplot(means,aes(x=atlas,node=variable,fill=variable,value=vp_percent)) +
  ggsankey::geom_sankey_bump(space=0,type="alluvial",alpha=.8,colour=NA,linewidth=0) +
  scale_fill_manual(values=variable_colours,drop=FALSE) +
  scale_x_discrete(expand=expansion(mult=c(.04,.04))) +
  scale_y_continuous(limits=c(0,max(12.5,ceiling(max(totals$vp_percent)/2.5)*2.5)),breaks=seq(0,max(12.5,ceiling(max(totals$vp_percent)/2.5)*2.5),2.5),
    labels=function(x)paste0(x,"%"),expand=expansion(mult=c(0,.01))) +
  labs(title="Community mean variance partitions",
    subtitle="Ribbon widths: mean scaled VP across the same 155 species; spatial partition excluded without renormalisation",
    x=NULL,y="Mean Tjur R²-scaled VP (%)",fill="Variable") +
  theme_minimal(base_size=12) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold"),legend.position="right",
    plot.margin=margin(12,12,12,12))
panel_b <- box_plot +
  labs(title="Species-level variance partitions",caption=NULL,
    subtitle="Only Holm-adjusted paired Wilcoxon p < 0.05 shown; adjustment across all 30 scaled-VP comparisons")
combined <- (panel_a / panel_b) + plot_layout(heights=c(.8,1)) +
  plot_annotation(title="Spatial models: scaled environmental variance partitioning through time",
    tag_levels="A",caption=paste0(
      "November 2 spatial models; ten environmental predictors displayed; 155 matched species. Scaled VP = original relative VP × species-specific Tjur R²; spatial partition omitted without renormalisation.\n",
      "Panel B: within each variable, 1970s, 1990s, 2010s from left to right. Boxes: median and IQR; whiskers: 1.5 × IQR; dots: outliers; white diamonds: means.\n",
      "Exploratory tests of paired shifts; posterior uncertainty and dependence among species are not accounted for."),
    theme=theme(plot.title=element_text(face="bold",size=17),
      plot.tag=element_text(face="bold",size=18),
      plot.caption=element_text(hjust=0,size=10)))
output_path <- file.path(output_dir,"vp-scaled-sankey-boxplots-spatial-models.png")
ggsave(output_path,combined,width=18,height=11,dpi=300,bg="white")
message("Saved: ",output_path,"; all 30 ribbon means checked against boxplot inputs.")

print(vp_total_check,row.names=FALSE)
