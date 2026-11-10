# Run from repository root: Rscript scripts/5-vis/env-viz/environmental-distributions.R
# Descriptive distributions of S0 predictors on the same cells in every atlas.
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
labels <- c(tmean_breeding='Breeding-season temperature (°C)',
  prec_breeding='Breeding-season precipitation (mm / season)',
  hh='Habitat heterogeneity (Simpson index)',
  perc_cropland='Cropland cover (%)',perc_urban='Urban cover (%)',
  perc_pasture='Pasture cover (%)',perc_grass_shrub='Grass / shrub cover (%)',
  perc_forest_deciduous='Deciduous forest cover (%)',
  perc_forest_evergreen='Evergreen forest cover (%)',
  perc_forest_mixed_other='Mixed / other forest cover (%)')
stopifnot(setequal(names(X),names(labels)))
period_labels <- c('Atlas 1\n1970s','Atlas 2\n1990s','Atlas 3\n2010s')
summary <- do.call(rbind,lapply(names(labels),function(v) {
  do.call(rbind,lapply(1:3,function(a) {
    z <- X[[v]][as.character(d$atlas)==as.character(a)]
    q <- quantile(z,c(0,.25,.5,.75,1),names=FALSE)
    data.frame(variable=v,atlas=a,n=length(z),min=q[1],q1=q[2],median=q[3],
               q3=q[4],max=q[5],zero_percent=100*mean(z==0))
  }))
}))
stopifnot(nrow(summary)==3*length(labels),all(summary$n==length(shared)),
          all(summary$min<=summary$q1),all(summary$q1<=summary$median),
          all(summary$median<=summary$q3),all(summary$q3<=summary$max))
summary$predictor <- factor(labels[summary$variable],levels=unname(labels))
summary$period <- factor(summary$atlas,levels=1:3,labels=period_labels)
p <- ggplot(summary,aes(x=period,fill=period)) +
  # Explicit full-range whiskers, rather than the default 1.5-IQR whiskers.
  geom_boxplot(aes(ymin=min,lower=q1,middle=median,upper=q3,ymax=max),
               stat='identity',width=.52,linewidth=.5,colour='grey25') +
  geom_point(aes(y=median),shape=21,fill='white',colour='black',size=2.1,stroke=.55) +
  facet_wrap(~predictor,ncol=2,scales='free_y') +
  scale_fill_manual(values=setNames(c('#74A9CF','#FEE08B','#EF8A62'),period_labels)) +
  scale_y_continuous(expand=expansion(mult=c(.04,.08))) +
  labs(title='Environmental distributions across atlas periods',
    subtitle=sprintf('%s identical training cells per period; predictors in original units',length(shared)),
    x=NULL,y=NULL,
    caption=paste('White dot / black line = median; box = 25th–75th percentiles; whiskers = minimum–maximum.',
      'Panels have different vertical scales. Collapsed boxes indicate tied quartiles (often zeros).',
      'These show variation among cells, not uncertainty in a period estimate.',sep='\n')) +
  theme_minimal(base_size=11) +
  theme(legend.position='none',panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),strip.text=element_text(face='bold',size=10),
    strip.background=element_rect(fill='grey92',colour=NA),
    panel.spacing=grid::unit(1,'lines'),
    plot.title=element_text(face='bold',size=17),
    plot.subtitle=element_text(margin=margin(b=12)),
    plot.caption=element_text(hjust=0,lineheight=1.15,margin=margin(t=12)),
    plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(14,18,14,14))
ggsave(file.path(out_dir,'environmental-distributions.png'),p,
       width=11,height=14,dpi=240,bg='white')
if(interactive()) print(p)
# Keep exact summary values inspectable alongside the PNG without extra data exports.
fmt <- function(z) format(signif(z,5),trim=TRUE,scientific=FALSE)
rows <- vapply(seq_len(nrow(summary)),function(i) {
  z <- summary[i,]
  paste0('| ',labels[[z$variable]],' | ',z$atlas,' | ',z$n,' | ',
    paste(vapply(z[c('min','q1','median','q3','max','zero_percent')],fmt,character(1)),collapse=' | '),' |')
},character(1))
writeLines(c('# Environmental distribution summaries','',
  paste('Source:',input),'',
  sprintf('Same %s cells in all three periods. Original units; sample quantiles use R type 7.',length(shared)),
  'Ranges describe variation among cells, not confidence intervals. Values below use five significant digits.','',
  '| Predictor | Atlas | Cells | Min | Q1 | Median | Q3 | Max | Zero (%) |',
  '|---|---:|---:|---:|---:|---:|---:|---:|---:|',rows),
  file.path(out_dir,'environmental-distribution-summary.md'))
message('Saved environmental-distributions.png and environmental-distribution-summary.md')
