# Run from repository root: Rscript scripts/5-vis/env-viz/land-use-sankey-bump.R
# Mean cell cover, original area percentages; not dominant-category cell counts
# or measured transitions between categories. PNG only.
suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
stopifnot(requireNamespace('ggsankey',quietly=TRUE))
out_dir <- 'scripts/5-vis/env-viz'
e <- new.env(parent=globalenv())
sys.source('scripts/1-atlas-analyses/homogenization/homogenization-data.R',envir=e)
sites <- e$sites
ids <- unlist(lapply(1:3,function(a)paste0(sites,'_',a)))
z <- e$environment[match(ids,e$environment$survey),,drop=FALSE]
stopifnot(identical(z$survey,ids),all(table(z$kvadratkod)==3),
          all(table(z$atlas)==length(sites)))
cover <- data.frame(Cropland=z$pct_class_22,Urban=z$pct_class_11,
  Pasture=z$pct_class_33,'Grass / shrub'=z$pct_class_55,
  'Deciduous forest'=z$pct_class_43+z$pct_class_44,
  'Evergreen forest'=z$pct_class_42+z$pct_class_41,
  'Mixed / other forest'=z$pct_class_40+z$pct_class_45,check.names=FALSE)
stopifnot(all(is.finite(as.matrix(cover))),all(as.matrix(cover)>=0),
          all(rowSums(cover)<=100+1e-6))
colours <- c(Cropland='goldenrod1',Urban='snow3',Pasture='darkorange',
  'Grass / shrub'='#B39DDB','Deciduous forest'='#A6D96A',
  'Evergreen forest'='#006400','Mixed / other forest'='#87CEEB')
plot_data <- do.call(rbind,lapply(1:3,function(a)data.frame(
  atlas=a,category=names(cover),cover_percent=colMeans(cover[z$atlas==a,,drop=FALSE]))))
plot_data$atlas <- factor(plot_data$atlas,levels=1:3,
  labels=c('Atlas 1\n1970s','Atlas 2\n1990s','Atlas 3\n2010s'))
plot_data$category <- factor(plot_data$category,levels=names(colours))
# Keep a numerical, reproducible summary in the console without extra output files.
summary <- data.frame(category=names(cover),
  atlas1=plot_data$cover_percent[1:7],atlas2=plot_data$cover_percent[8:14],
  atlas3=plot_data$cover_percent[15:21])
summary$change_percentage_points <- summary$atlas3-summary$atlas1
print(summary,row.names=FALSE,digits=4)
print(tapply(plot_data$cover_percent,plot_data$atlas,sum))
p <- ggplot(plot_data,aes(x=atlas,node=category,fill=category,value=cover_percent)) +
  ggsankey::geom_sankey_bump(space=0,type='alluvial',alpha=.9,colour=NA,linewidth=0) +
  scale_fill_manual(values=colours,breaks=names(colours)) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),
    labels=function(x)paste0(x,'%'),expand=expansion(mult=0)) +
  scale_x_discrete(expand=expansion(mult=c(.025,.025))) +
  labs(title='Land-use composition across the three atlases',
    subtitle=sprintf('%s identical cells · Seven model land-cover categories',
      format(length(sites),big.mark=',')),
    x=NULL,y='Mean percentage of cell area',fill=NULL,
    caption=paste('Each cell has equal weight; percentages retain the original valid sampled cell-area denominator.',
      'Water and other omitted classes account for the gap to 100%; forest types remain separate.',
      'Ribbons follow category shares through time; they do not represent land converted between categories.',sep='\n')) +
  theme_minimal(base_size=12) + theme(legend.position='right',
    panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
    axis.text.x=element_text(face='bold',size=12),
    plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,lineheight=1.15),
    plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(12,12,12,12))
ggsave(file.path(out_dir,'land-use-sankey-bump.png'),p,width=12,height=7,dpi=300,bg='white')
