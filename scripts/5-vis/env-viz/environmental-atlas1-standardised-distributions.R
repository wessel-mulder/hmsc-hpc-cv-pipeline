# Run from root: Rscript scripts/5-vis/env-viz/environmental-atlas1-standardised-distributions.R
# Atlas 1 mean and sample SD applied unchanged to all three periods.
suppressPackageStartupMessages(library(ggplot2))
input <- 'data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
out_dir <- 'scripts/5-vis/env-viz'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
e <- new.env(parent=emptyenv()); load(input,envir=e)
X <- as.data.frame(e$X); d <- e$Design
stopifnot(identical(rownames(X),rownames(d)),
          !anyDuplicated(d[c('site','atlas')]),all(is.finite(as.matrix(X))))
shared <- sort(Reduce(intersect,lapply(1:3,function(a)
  as.character(d$site[as.character(d$atlas)==as.character(a)]))))
keep <- as.character(d$site) %in% shared
X <- X[keep,,drop=FALSE]; d <- droplevels(d[keep,,drop=FALSE])
stopifnot(length(shared)>0,all(table(d$site)==3),all(table(d$atlas)==length(shared)))
labels <- c(tmean_breeding='Breeding-season temperature',
  prec_breeding='Breeding-season precipitation',hh='Habitat heterogeneity',
  perc_cropland='Cropland cover',perc_urban='Urban cover',
  perc_pasture='Pasture cover',perc_grass_shrub='Grass / shrub cover',
  perc_forest_deciduous='Deciduous forest cover',
  perc_forest_evergreen='Evergreen forest cover',
  perc_forest_mixed_other='Mixed / other forest cover')
units <- c(tmean_breeding='°C',prec_breeding='mm / season',hh='Simpson index',
           setNames(rep('%',7),names(labels)[4:10]))
stopifnot(setequal(names(X),names(labels)))
baseline <- as.character(d$atlas)=='1'
pooled_mean <- vapply(X[baseline,,drop=FALSE],mean,numeric(1))
pooled_sd <- vapply(X[baseline,,drop=FALSE],sd,numeric(1))
stopifnot(all(is.finite(pooled_sd)),all(pooled_sd>0))
Z <- sweep(sweep(as.matrix(X),2,pooled_mean,'-'),2,pooled_sd,'/')
stopifnot(max(abs(colMeans(Z[baseline,,drop=FALSE])))<1e-10,
          max(abs(apply(Z[baseline,,drop=FALSE],2,sd)-1))<1e-10)
stats <- do.call(rbind,lapply(names(labels),function(v) {
  do.call(rbind,lapply(1:3,function(a) {
    q <- quantile(Z[as.character(d$atlas)==as.character(a),v],
                  c(0,.25,.5,.75,1),names=FALSE)
    data.frame(variable=v,atlas=a,min=q[1],q1=q[2],median=q[3],q3=q[4],max=q[5])
  }))
}))
fmt <- function(z) formatC(z,format='f',digits=2)
axis_labels <- setNames(vapply(names(labels),function(v)
  sprintf('%s\nMean %s %s; SD %s %s',labels[[v]],fmt(pooled_mean[[v]]),units[[v]],
          fmt(pooled_sd[[v]]),units[[v]]),character(1)),names(labels))
# Reverse levels places climate first at the top after flipping the coordinates.
stats$variable <- factor(stats$variable,levels=rev(names(labels)))
# Reverse the factor order so each horizontal group reads Atlas 1, 2, 3 top-down.
stats$period <- factor(stats$atlas,levels=3:1,
  labels=c('Atlas 3 — 2010s','Atlas 2 — 1990s','Atlas 1 — 1970s'))
colours <- c('Atlas 1 — 1970s'='#74A9CF','Atlas 2 — 1990s'='#FEE08B',
             'Atlas 3 — 2010s'='#EF8A62')
dodge <- position_dodge(width=.75)
cap <- function(x) pmax(-8,pmin(8,x))
p <- ggplot(stats,aes(x=variable,fill=period,group=interaction(variable,period))) +
  geom_hline(yintercept=0,colour='grey45',linetype='dashed',linewidth=.45) +
  geom_boxplot(aes(ymin=cap(min),lower=cap(q1),middle=cap(median),upper=cap(q3),ymax=cap(max)),
    stat='identity',position=dodge,width=.65,linewidth=.4,colour='grey35') +
  geom_point(aes(y=cap(median),alpha=I(as.numeric(median>=-8 & median<=8))),position=dodge,shape=21,fill='white',
             colour='black',size=1.65,stroke=.45,show.legend=FALSE) +
  # Keep every group in these layers so dodge positions match the boxes.
  # Invisible points mark in-range endpoints; coloured dots flag clipped tails.
  geom_point(aes(y=-8,alpha=I(as.numeric(min < -8))),position=dodge,
             shape=21,colour='grey25',size=2.2,stroke=.4,show.legend=FALSE) +
  geom_point(aes(y=8,alpha=I(as.numeric(max > 8))),position=dodge,
             shape=21,colour='grey25',size=2.2,stroke=.4,show.legend=FALSE) +
  scale_x_discrete(labels=axis_labels) +
  scale_fill_manual(values=colours,breaks=names(colours),name=NULL) +
  scale_y_continuous(breaks=seq(-8,8,by=2),
                     expand=expansion(mult=0)) +
  coord_flip(ylim=c(-8,8),clip='off') +
  labs(title='Environmental distributions relative to Atlas 1',
    subtitle=sprintf('%s identical cells per atlas; all periods use the Atlas 1 mean and SD',length(shared)),
    x=NULL,y='Deviation from Atlas 1 mean (Atlas 1 standard deviations)',
    caption=paste('White dot / black line = median; box = middle 50%; whiskers = minimum–maximum, clipped at ±8 SD.',
      'Coloured dots at −8 or +8 flag observed minima or maxima beyond the displayed range.',
      'Labels give Atlas 1 mean and sample SD in original units; medians outside ±8 SD are not drawn.',
      'Zero marks the Atlas 1 mean. Ranges describe variation among cells, not uncertainty.',sep='\n')) +
  theme_minimal(base_size=10.5) +
  theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
    axis.text.y=element_text(size=9.5,lineheight=1.15,colour='grey15'),
    axis.title.x=element_text(margin=margin(t=10)),
    legend.position='top',legend.justification='left',
    legend.margin=margin(0,0,5,0),
    plot.title=element_text(face='bold',size=17),
    plot.subtitle=element_text(margin=margin(b=6)),
    plot.caption=element_text(hjust=0,lineheight=1.15,margin=margin(t=12)),
    plot.background=element_rect(fill='white',colour=NA),
    plot.margin=margin(14,16,12,14))
ggsave(file.path(out_dir,'environmental-atlas1-standardised-distributions.png'),p,
       width=12,height=9,dpi=240,bg='white')
if(interactive()) print(p)
