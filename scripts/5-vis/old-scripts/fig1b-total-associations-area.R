# Simplified view of the verified final-model-cell area diagnostic.
suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
out <- 'manuscript-outputs/figures/exploratory'
dat <- read.csv(file.path(out,'fig1b-habitat-area-diagnostic.csv'))
totals <- dat |> group_by(period,variable) |> summarise(
  n_associations=sum(supported_species),
  estimated_habitat_km2=first(estimated_habitat_km2),
  n_area_values=n_distinct(estimated_habitat_km2),.groups='drop') |>
  mutate(period=factor(period,levels=c('1970s','1990s','2010s'))) |>
  arrange(variable,period)
stopifnot(all(totals$n_area_values==1),all(totals$n_associations<=157),nrow(totals)==15)
cols <- c(Urban='grey55',Cropland='goldenrod',Pasture='darkorange',Forest='#008b45','Grass/shrubland'='#23c776')
p <- ggplot(totals,aes(estimated_habitat_km2,n_associations,colour=variable,group=variable))+
  geom_path(linewidth=.65,alpha=.6)+geom_point(aes(alpha=period),size=4)+
  scale_colour_manual(values=cols)+scale_alpha_manual(values=c('1970s'=.4,'1990s'=.7,'2010s'=1))+
  scale_x_continuous(labels=scales::label_comma(),limits=c(0,NA),expand=expansion(mult=c(.02,.04)))+
  scale_y_continuous(limits=c(0,157),breaks=c(0,25,50,75,100,125,157))+
  labs(x='Estimated habitat area in modelled cells (km²)',
       y='Number of supported associations\n(positive + negative)',colour='Land cover',alpha='Atlas period',
       caption='Counts out of 157 modelled species. Area = sum of cover fraction × cell area, using only final-model cells.\nLines connect atlas periods within each land-cover type.')+
  guides(alpha=guide_legend(order=1,override.aes=list(colour='black')),
         colour=guide_legend(order=2,nrow=1,override.aes=list(alpha=1)))+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position='bottom',
    legend.box='vertical',plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(out,'fig1b-total-associations-vs-habitat-area.png'),p,width=10,height=7,dpi=300,bg='white')
write.csv(select(totals,-n_area_values),file.path(out,'fig1b-total-associations-vs-habitat-area.csv'),row.names=FALSE)
