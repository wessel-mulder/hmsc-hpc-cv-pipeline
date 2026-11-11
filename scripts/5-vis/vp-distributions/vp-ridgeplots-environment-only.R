# Run from project root:
# Rscript scripts/5-vis/vp-distributions/vp-ridgeplots-environment-only.R
# Generates only the scaled two-column ridge PNG; keeps scaled boxplots in
# memory for the combined figure. Species estimates, not VP posterior draws.
required <- c("ggplot2", "ggridges")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Install: ",paste(missing,collapse=", "))
suppressPackageStartupMessages({library(ggplot2); library(ggridges)})
model_root <- "data-new/2-hmsc/derived/models"
output_dir <- "scripts/5-vis/vp-distributions"
model_stamp <- "2026-11-06_05-11-30"
atlas_labels <- c("1970s","1990s","2010s")
variable_labels <- c(tmean_breeding="Temperature",prec_breeding="Precipitation",
  hh="Land-use heterogeneity",perc_urban="Urban cover",perc_cropland="Cropland cover",
  perc_pasture="Pasture cover",perc_forest_deciduous="Deciduous forest cover",
  perc_forest_evergreen="Evergreen forest cover",perc_forest_mixed_other="Mixed/other forest cover",
  perc_grass_shrub="Grass/shrubland cover")
variable_colours <- c("Temperature"="firebrick3","Precipitation"="dodgerblue3",
  "Land-use heterogeneity"="orchid3","Urban cover"="snow3","Cropland cover"="goldenrod1",
  "Pasture cover"="darkorange","Deciduous forest cover"="#A6D96A",
  "Evergreen forest cover"="#006400","Mixed/other forest cover"="#87CEEB",
  "Grass/shrubland cover"="#B39DDB")
species <- NULL; rows <- list(); total_rows <- list()
for(a in 1:3) {
  folder <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage_EnvironmentOnly",model_stamp,a)
  path <- file.path(model_root,folder,"Results",paste0(folder,"parameter_estimates_VP_.csv"))
  raw <- read.csv(path,row.names=1,check.names=FALSE)
  stopifnot(setequal(rownames(raw),c(names(variable_labels),"TjurR2")),
    ncol(raw)==155L,!anyDuplicated(colnames(raw)),
    all(vapply(raw,is.numeric,logical(1))),all(is.finite(as.matrix(raw))))
  if(is.null(species)) species <- colnames(raw)
  stopifnot(setequal(colnames(raw),species))
  raw <- raw[,species,drop=FALSE]
  relative <- as.matrix(raw[names(variable_labels),,drop=FALSE])
  tjur <- as.numeric(raw["TjurR2",])
  stopifnot(all(relative>=0),all(relative<=1),max(abs(colSums(relative)-1))<1e-8,
    all(tjur>=0),all(tjur<=1))
  scaled <- sweep(relative,2,tjur,"*")
  stopifnot(abs(sum(rowMeans(scaled))-mean(tjur))<1e-10)
  rows[[a]] <- do.call(rbind,lapply("scaled",function(scale_name) {
    values <- if(scale_name=="scaled") scaled else relative
    data.frame(atlas=atlas_labels[a],variable=rep(names(variable_labels),times=155L),
      species=rep(species,each=10L),scale=scale_name,vp_percent=as.vector(values)*100)
  }))
  total_rows[[a]] <- data.frame(atlas=a,mean_tjur_percent=mean(tjur)*100,
    summed_scaled_means_percent=sum(rowMeans(scaled))*100)
}
species_vp <- do.call(rbind,rows)
vp_total_check <- do.call(rbind,total_rows)
stopifnot(nrow(species_vp)==4650L,all(table(species_vp$atlas,species_vp$variable,species_vp$scale)==155L))

# Use the scaled-only 30-test family and species pairing from the Wilcoxon helper.
wilcoxon_analysis <- new.env()
invisible(capture.output(source(file.path(output_dir,"vp-wilcoxon-environment-only.R"),
  local=wilcoxon_analysis)))
ridge_results <- list()
for(scale_name in "scaled") {
  dat <- species_vp[species_vp$scale==scale_name,]
  bandwidths <- vapply(names(variable_labels),function(v) bw.nrd0(dat$vp_percent[dat$variable==v]),numeric(1))
  stopifnot(all(is.finite(bandwidths)),all(bandwidths>0))
  # The same horizontal scale across variables/atlases includes all observations
  # and four bandwidths of upper tails where the theoretical upper bound permits.
  x_max <- min(100,ceiling((max(dat$vp_percent)+4*max(bandwidths))/5)*5)
  grid_x <- seq(0,x_max,length.out=2048)
  density_rows <- lapply(names(variable_labels),function(v) {
    panel <- do.call(rbind,lapply(atlas_labels,function(a) {
      x <- dat$vp_percent[dat$variable==v & dat$atlas==a]
      stopifnot(length(x)==155L)
      bw <- bandwidths[[v]]
      # Reflection at zero and 100 percentage points avoids density below zero.
      kernel <- function(centres) rowMeans(dnorm(outer(grid_x,centres,"-"),sd=bw))
      y <- kernel(x)+kernel(-x)+kernel(200-x)
      mass <- sum(diff(grid_x)*(head(y,-1)+tail(y,-1))/2)
      stopifnot(all(is.finite(y)),all(y>=0),abs(mass-1)<.01)
      data.frame(variable=unname(variable_labels[v]),atlas=a,x=grid_x,density=y,bandwidth=bw)
    }))
    panel$height <- panel$density/max(panel$density)
    panel
  })
  ridge_data <- do.call(rbind,density_rows)
  ridge_data$variable <- factor(ridge_data$variable,levels=unname(variable_labels))
  ridge_data$atlas <- factor(ridge_data$atlas,levels=rev(atlas_labels))
  p <- ggplot(ridge_data,aes(x=x,y=atlas,height=height,fill=variable,alpha=atlas,group=atlas)) +
    ggridges::geom_ridgeline(stat="identity",scale=.8,colour=NA,linewidth=0,min_height=0) +
    facet_grid(variable~.,switch="y") +
    scale_fill_manual(values=variable_colours) +
    scale_alpha_manual(values=c("1970s"=.45,"1990s"=.7,"2010s"=1)) +
    scale_x_continuous(limits=c(0,x_max),breaks=pretty(c(0,x_max),n=8),
      labels=function(x)paste0(format(x,trim=TRUE,scientific=FALSE),"%"),
      expand=expansion(mult=c(0,0))) +
    scale_y_discrete(expand=expansion(add=c(.15,1))) +
    labs(title=if(scale_name=="scaled") "Environment-only models: Tjur R²-scaled variance partitions"
               else "Environment-only models: relative variance partitions",
      subtitle="Each ridge: 155 species' VP estimates; three atlas ridges per variable",
      x=if(scale_name=="scaled") "Tjur R²-scaled VP (%)" else "Relative VP (% of partitioned model-explained variation)",
      y=NULL,
      caption=paste0("November 6 models; ten environmental predictors; all 155 matched species, equally weighted. No outlines.\n",
        if(scale_name=="scaled") "Each species' exported relative VP is multiplied by its Tjur R², matching the Sankey calculation.\n"
        else "Environmental partitions sum to 100% within each species; these models have no spatial random effect.\n",
        "Ridges describe variation among species, not posterior uncertainty. Shared bandwidth within each variable; reflected boundaries.\n",
        "Common horizontal scale across all rows; peak heights normalised within variable. VP shifts do not directly test coefficient changes.")) +
    theme_minimal(base_size=12) +
    theme(legend.position="none",panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
      strip.placement="outside",strip.background=element_rect(fill="grey95",colour=NA),
      strip.text.y.left=element_text(angle=0,hjust=1,face="bold",size=10),axis.text.y=element_text(size=9),
      panel.spacing.y=grid::unit(.2,"lines"),plot.title=element_text(face="bold"),
      plot.caption=element_text(hjust=0,size=9),plot.margin=margin(12,12,12,12))
  two_columns <- p +
    facet_wrap(~variable,nrow=5,ncol=2,dir="v",strip.position="left",
      axes="all_y",axis.labels="all_y") +
    theme(panel.spacing.x=grid::unit(1,"lines"),panel.spacing.y=grid::unit(.35,"lines"))
  layout <- ggplot_build(two_columns)$layout$layout
  stopifnot(identical(as.character(layout$variable[layout$COL==1]),unname(variable_labels[1:5])),
    identical(as.character(layout$variable[layout$COL==2]),unname(variable_labels[6:10])))
  two_column_path <- file.path(output_dir,paste0("vp-ridgeplots-environment-only-",scale_name,"-two-columns.png"))
  ggsave(two_column_path,two_columns,width=18,height=10.5,dpi=300,bg="white")
  box_data <- dat
  box_data$variable <- factor(box_data$variable,levels=names(variable_labels),labels=unname(variable_labels))
  box_data$atlas <- factor(box_data$atlas,levels=atlas_labels)
  atlas_offset <- if(scale_name=="scaled") 1.3 else 1.1
  box_data$x_position <- as.integer(box_data$variable)*4 +
    c(-atlas_offset,0,atlas_offset)[as.integer(box_data$atlas)]
  box_data$box_group <- interaction(box_data$variable,box_data$atlas,drop=TRUE)
  brackets <- wilcoxon_analysis$wilcoxon_results
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
  stopifnot(nrow(brackets)==sum(wilcoxon_analysis$wilcoxon_results$scale==scale_name &
    wilcoxon_analysis$wilcoxon_results$p_holm<.05),all(brackets$p_holm<.05),!anyNA(brackets$y))
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
    labs(title=if(scale_name=="scaled") "Environment-only models: Tjur R²-scaled variance partitions"
               else "Environment-only models: relative variance partitions",
      subtitle="155 matched species per box; only Holm-adjusted paired Wilcoxon p < 0.05 shown (30-test family)",
      x=NULL,y=if(scale_name=="scaled") "Tjur R²-scaled VP (%)" else "Relative VP (%)",
      caption=paste0("November 6 models; within each variable: 1970s, 1990s, 2010s from left to right. Lines connect compared atlas boxes.\n",
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
  ridge_results[[scale_name]] <- list(plot=p,two_column_plot=two_columns,
    box_plot=box_plot,data=ridge_data,x_max=x_max,bandwidths=bandwidths)
  message("Saved: ",two_column_path,"; first five variables down left, remaining five down right.")
  message("Kept scaled boxplots in memory for panel B; all 30 medians checked.")
}
message("Validated all 30 ridges per scale; species scaling reproduces the environment-only Sankey totals.")
print(vp_total_check,row.names=FALSE)
