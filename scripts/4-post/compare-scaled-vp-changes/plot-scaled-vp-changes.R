# Run from the repository root:
# Rscript scripts/4-post/compare-scaled-vp-changes/plot-scaled-vp-changes.R
# Reads saved posterior results only; does not rerun models or comparisons.
if(!requireNamespace("ggplot2",quietly=TRUE)) stop("Required package: ggplot2")
suppressPackageStartupMessages(library(ggplot2))
output_dir <- "scripts/4-post/compare-scaled-vp-changes"
input <- file.path(output_dir,"comparison.rds")
results <- readRDS(input)
dat <- results$community_summary
variable_labels <- c(tmean_breeding="Temperature",prec_breeding="Precipitation",
  hh="Land-use heterogeneity",perc_urban="Urban cover",perc_cropland="Cropland cover",
  perc_pasture="Pasture cover",perc_forest_deciduous="Deciduous forest cover",
  perc_forest_evergreen="Evergreen forest cover",perc_forest_mixed_other="Mixed/other forest cover",
  perc_grass_shrub="Grass/shrubland cover")
pair_levels <- c("1970s to 1990s","1990s to 2010s","1970s to 2010s")
stopifnot(nrow(dat)==30L,setequal(dat$variable,names(variable_labels)),
  setequal(dat$pair,pair_levels),!anyDuplicated(dat[c("pair","variable")]),
  all(is.finite(as.matrix(dat[c("mean","joint_lower","joint_upper")]))),
  all(dat$joint_lower<=dat$joint_upper),results$joint_envelope$empirical_coverage>=.95-1e-10)
dat$direction <- ifelse(dat$joint_lower>0,"Increase",
  ifelse(dat$joint_upper<0,"Decrease","Uncertain"))
dat$direction[dat$source_diagnostic_status!="pass"] <- "Check diagnostics"
dat$direction <- factor(dat$direction,levels=c("Increase","Decrease","Uncertain","Check diagnostics"))
dat$variable_label <- factor(unname(variable_labels[dat$variable]),levels=rev(unname(variable_labels)))
dat$pair_label <- factor(dat$pair,levels=pair_levels,
  labels=c("1970s → 1990s","1990s → 2010s","1970s → 2010s"))
change <- ifelse(abs(dat$mean)<.005,"0.00",sprintf("%+.2f",dat$mean))
dat$label <- paste0(as.character(dat$direction),"\n",change," pp")
dat$text_colour <- ifelse(dat$direction %in% c("Increase","Decrease"),"white","#333333")
counts <- table(dat$direction)
subtitle <- sprintf("Joint 95%% posterior intervals across 30 comparisons: %d increases, %d decreases, %d uncertain",
  counts["Increase"],counts["Decrease"],counts["Uncertain"])
if(!isTRUE(results$provenance$community_diagnostics_all_pass)) subtitle <- paste0(subtitle,
  "\nCaution: source convergence flags limit interpretation of the joint envelope")
p <- ggplot(dat,aes(x=pair_label,y=variable_label,fill=direction)) +
  geom_tile(width=.96,height=.94,colour="white",linewidth=.6) +
  geom_text(aes(label=label,colour=text_colour),size=4.1,lineheight=1.1) +
  scale_colour_identity() +
  scale_fill_manual(values=c("Increase"="#20853F","Decrease"="#B2182B",
    "Uncertain"="#EEEEEE","Check diagnostics"="#FFD08A"),drop=TRUE) +
  scale_x_discrete(position="top",expand=expansion(add=.02)) +
  scale_y_discrete(expand=expansion(add=.02)) +
  labs(title="Changes in community mean scaled VP",subtitle=subtitle,x=NULL,y=NULL,fill=NULL,
    caption=paste0("Cell values: posterior mean change in percentage points (pp), later atlas minus earlier atlas.\n",
      "Green/red: the joint 95% interval lies entirely above/below zero. Grey: the interval includes zero; this does not establish no change.\n",
      "November 6 environment-only models; 155 matched species; each atlas's fitted cells. Scaled VP = relative VP × species-specific Tjur R².\n",
      "These are model-based changes in explanatory importance, not causal effects or tests of changing response slopes.")) +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(),axis.ticks=element_blank(),
    axis.text.x=element_text(size=12,face="bold",colour="grey20",margin=margin(b=9)),
    axis.text.y=element_text(size=11,colour="grey20",margin=margin(r=8)),
    legend.position="bottom",legend.text=element_text(size=11),
    plot.title=element_text(face="bold",size=18),plot.subtitle=element_text(size=11,margin=margin(b=15)),
    plot.caption=element_text(hjust=0,size=9,lineheight=1.15,margin=margin(t=10)),
    plot.margin=margin(15,18,15,15))
output <- file.path(output_dir,"scaled-vp-community-change-directions.png")
ggsave(output,p,width=12,height=9.5,dpi=300,bg="white")
stopifnot(nrow(ggplot_build(p)$data[[1]])==30L)
print(table(dat$pair,dat$direction))
message("Saved: ",normalizePath(output))
