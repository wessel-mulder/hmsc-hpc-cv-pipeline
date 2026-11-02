# Simplified view of the verified final-model-cell area diagnostic.
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(readxl);library(cowplot)})
out <- 'manuscript-outputs/figures/exploratory'
dat <- read.csv(file.path(out,'fig1b-habitat-area-diagnostic.csv'))
totals <- dat |> group_by(period,variable) |> summarise(
  n_associations=sum(supported_species),
  mean_cover_percent=first(mean_cover_percent),
  n_cover_values=n_distinct(mean_cover_percent),.groups='drop') |>
  mutate(period=factor(period,levels=c('1970s','1990s','2010s'))) |>
  arrange(variable,period)
stopifnot(all(totals$n_cover_values==1),all(totals$n_associations<=157),nrow(totals)==15)
cols <- c(Urban='grey55',Cropland='goldenrod',Pasture='darkorange',Forest='#008b45','Grass/shrubland'='#23c776')
p <- ggplot(totals,aes(mean_cover_percent,n_associations,colour=variable,group=variable))+
  geom_path(linewidth=.65,alpha=.6)+geom_point(aes(alpha=period),size=4)+
  scale_colour_manual(values=cols)+scale_alpha_manual(values=c('1970s'=.4,'1990s'=.7,'2010s'=1))+
  scale_x_continuous(labels=scales::label_number(suffix="%"),limits=c(0,NA),expand=expansion(mult=c(.02,.04)))+
  scale_y_continuous(limits=c(0,157),breaks=c(0,25,50,75,100,125,157))+
  labs(x='Mean land cover (%)',
       y='Number of supported associations\n(positive + negative)',colour='Land cover',alpha='Atlas period',
       caption='Counts out of 157 modelled species. Cover is the arithmetic mean percentage across final-model cells.\nLines connect atlas periods within each land-cover type.')+
  guides(alpha=guide_legend(order=1,override.aes=list(colour='black')),
         colour=guide_legend(order=2,nrow=1,override.aes=list(alpha=1)))+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position='bottom',
    legend.box='vertical',plot.caption=element_text(hjust=0,size=9))
# Heterogeneity uses the identical final-model cell lists, on its native index scale.
e <- new.env()
load('data-new/derived/hmsc-inputs/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData',e)
cells <- read.csv(file.path(out,'fig1b-habitat-area-model-cells.csv'))
folders <- list.dirs('intermediate-outputs/hmsc-output',recursive=FALSE,full.names=TRUE)
folders <- folders[grepl('2026-07-27.*Atlas[123]_MinOccs5_CoverageGoodAverage$',folders)]
hh <- bind_rows(lapply(folders,function(folder){
  a <- as.integer(sub('.*Atlas([123])_.*','\\1',folder))
  period <- c('1970s','1990s','2010s')[a]
  surveys <- cells$survey[cells$period==period]
  stopifnot(all(surveys %in% rownames(e$X)))
  values <- e$X[surveys,'hh']
  probs <- read_excel(file.path(folder,'Results',paste0(basename(folder),'parameter_estimates_Beta_.xlsx')),sheet='Pr(x>0)')
  stopifnot(nrow(probs)==157,!anyNA(probs$hh),!anyNA(values))
  tibble(period=period,variable='Land-use heterogeneity',mean_index=mean(values),
    n_associations=sum(probs$hh>=.95 | probs$hh<=.05))
})) |> mutate(period=factor(period,levels=c('1970s','1990s','2010s'))) |> arrange(period)
q <- ggplot(hh,aes(mean_index,n_associations,group=1))+
  geom_path(colour='orchid3',linewidth=.65,alpha=.6)+
  geom_point(aes(alpha=period),colour='orchid3',size=4)+
  scale_alpha_manual(values=c('1970s'=.4,'1990s'=.7,'2010s'=1))+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  scale_y_continuous(limits=c(0,157),breaks=c(0,25,50,75,100,125,157))+
  labs(x='Mean heterogeneity index',y=NULL,title='Land-use heterogeneity')+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position='none',
    plot.title=element_text(size=12,face='bold'))
legend <- get_legend(p+labs(caption=NULL))
body <- plot_grid(p+labs(caption=NULL,title='Land cover')+theme(legend.position='none',plot.title=element_text(size=12,face='bold')),
  q,nrow=1,rel_widths=c(2,1),align='hv',axis='tblr')
footer <- ggdraw()+draw_label('Counts out of 157 modelled species; means calculated across the exact final-model cells.\nLines connect atlas periods. Heterogeneity is shown on its native 0–1 scale.',size=10)
figure <- plot_grid(body,legend,footer,ncol=1,rel_heights=c(1,.23,.09))
ggsave(file.path(out,'fig1b-total-associations-vs-mean-cover.png'),figure,width=13,height=7,dpi=300,bg='white')
write.csv(select(totals,-n_cover_values),file.path(out,'fig1b-total-associations-vs-mean-cover.csv'),row.names=FALSE)
write.csv(hh,file.path(out,'fig1b-total-associations-vs-heterogeneity.csv'),row.names=FALSE)
print(hh)
