# Direct implementation of the final JSDM Methods paragraph in Antao et al. (2022).
# Run from project root with Rscript scripts/4-data-visualisation/atlas-vp-paper-lm.R
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(cowplot)})
vars <- c(tmean_breeding='Temperature',prec_breeding='Precipitation',hh='Land-use heterogeneity',perc_urban='Urban cover',perc_cropland='Cropland cover',perc_pasture='Pasture cover',perc_forest='Forest cover',perc_grass_shrub='Grass / shrub cover')
out <- 'intermediate-outputs/vp-paper-lm'
fig <- 'manuscript-outputs/figures/exploratory/vp-paper-lm'
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(fig,recursive=TRUE,showWarnings=FALSE)
folders<-sort(list.dirs('intermediate-outputs/hmsc-output',recursive=FALSE))
folders<-folders[grepl('^2026-07-27.*_Atlas[123]_MinOccs5_CoverageGoodAverage$',basename(folders))]
stopifnot(length(folders)==3)
data<-bind_rows(lapply(folders,function(folder){
  atlas<-as.integer(sub('.*Atlas([123])_.*','\\1',basename(folder)))
  a<-as.matrix(read.csv(file.path(folder,'Results',paste0(basename(folder),'parameter_estimates_VP_.csv')),row.names=1,check.names=FALSE))
  vp<-readRDS(file.path(folder,paste0(basename(folder),'VP.rds')))$vals
  mf<-readRDS(file.path(folder,paste0(basename(folder),'MF.rds')))
  stopifnot(identical(colnames(a),colnames(vp)),max(abs(a[rownames(vp),]-vp))<1e-10,max(abs(a['TjurR2',]-mf$TjurR2))<1e-10)
  bind_rows(lapply(names(vars),function(v)data.frame(species=colnames(a),variable=v,atlas=atlas,
    decade_since_1970s=c(0,2,4)[atlas],vp=as.numeric(a[v,]),tjur_r2=as.numeric(a['TjurR2',]),weighted_vp=as.numeric(a[v,]*a['TjurR2',]))))
}))
stopifnot(!anyNA(data),!anyDuplicated(data[c('species','variable','atlas')]))
# One pooled linear regression per covariate. Species-atlas estimates are the observations.
# This follows the paper's lm specification: no species term, paired test or p adjustment.
models<-setNames(lapply(names(vars),function(v)lm(weighted_vp ~ decade_since_1970s,data=filter(data,variable==v))),names(vars))
results<-bind_rows(lapply(names(models),function(v){
  m<-models[[v]]; s<-summary(m); ci<-confint(m)
  data.frame(variable=v,label=vars[[v]],n=nobs(m),n_species=n_distinct(filter(data,variable==v)$species),
    intercept=unname(coef(m)[1]),slope_per_decade=unname(coef(m)[2]),se=s$coefficients[2,2],
    ci_low=ci[2,1],ci_high=ci[2,2],t=s$coefficients[2,3],df=df.residual(m),p=s$coefficients[2,4],r_squared=s$r.squared)
})) |> mutate(stars=case_when(p<.001~'***',p<.01~'**',p<.05~'*',TRUE~''),direction=if_else(slope_per_decade>0,'Increase','Decrease'))
write.csv(data,file.path(out,'species-atlas-weighted-vp.csv'),row.names=FALSE)
write.csv(results,file.path(out,'linear-regression-results.csv'),row.names=FALSE)
saveRDS(models,file.path(out,'linear-regression-models.rds'))
cols<-c(Increase='#B2182B',Decrease='#2166AC')
panels<-lapply(names(vars),function(v){
  d<-filter(data,variable==v);r<-filter(results,variable==v);m<-models[[v]]
  nd<-data.frame(decade_since_1970s=seq(0,4,length.out=100))
  pred<-cbind(nd,as.data.frame(predict(m,nd,interval='confidence')))
  ggplot(d,aes(decade_since_1970s,weighted_vp))+
    geom_violin(aes(group=atlas),fill='grey92',colour='grey75',width=.9,linewidth=.3)+
    geom_point(position=position_jitter(width=.16,height=0,seed=19),colour='grey35',alpha=.4,size=.7)+
    geom_ribbon(data=pred,aes(x=decade_since_1970s,ymin=lwr,ymax=upr),inherit.aes=FALSE,fill=cols[[r$direction]],alpha=.18)+
    geom_line(data=pred,aes(decade_since_1970s,fit),inherit.aes=FALSE,colour=cols[[r$direction]],linewidth=.9)+
    annotate('text',x=0,y=Inf,hjust=0,vjust=1.2,size=3.1,label=sprintf('Slope = %+.4f / decade\np = %s %s',r$slope_per_decade,format.pval(r$p,digits=3,eps=.001),r$stars))+
    scale_x_continuous(breaks=c(0,2,4),labels=c('1970s','1990s','2010s'),limits=c(-.55,4.55))+
    scale_y_continuous(expand=expansion(mult=c(.02,.2)))+
    labs(title=vars[[v]],x=NULL,y='Tjur-weighted VP')+
    theme_minimal(base_size=10)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),plot.title=element_text(face='bold',hjust=.5),plot.margin=margin(8,8,8,8))
})
header<-ggdraw()+draw_label('Temporal changes in environmental variance components',fontface='bold',size=18)
subtitle<-ggdraw()+draw_label('Paper-method regressions: species-level VP × Tjur’s R² ~ decade',size=12)
foot<-ggdraw()+draw_label('157 species × 3 atlases per regression. Lines and ribbons: lm fit and 95% confidence interval.\nUnadjusted OLS p-values: * p < 0.05; ** p < 0.01; *** p < 0.001. Panel y-scales differ.',size=10)
plot<-plot_grid(header,subtitle,plot_grid(plotlist=panels,ncol=4),foot,ncol=1,rel_heights=c(.065,.04,1,.075))
ggsave(file.path(fig,'vp-paper-lm-trajectories.png'),plot,width=15,height=9,dpi=220,bg='white')
heat<-results |> mutate(variable=factor(variable,levels=names(vars)))
lim<-max(abs(heat$slope_per_decade))
h<-ggplot(heat,aes(variable,'All modelled birds',fill=slope_per_decade))+
  geom_tile(colour='white',linewidth=1)+geom_text(aes(label=sprintf('%+.4f\n%s',slope_per_decade,stars)),size=4)+
  scale_x_discrete(labels=c('Temperature','Precipitation','Land-use\nheterogeneity','Urban','Cropland','Pasture','Forest','Grass /\nshrub'))+
  scale_fill_gradient2(low='#4393C3',mid='white',high='#D6604D',limits=c(-lim,lim),breaks=c(-.002,0,.002),name='Slope / decade',guide=guide_colourbar(barwidth=grid::unit(6,'cm')))+
  labs(title='Temporal slopes of Tjur-weighted VP',subtitle='Linear regressions following the final JSDM Methods paragraph of Antão et al. (2022)',x=NULL,y=NULL,caption='* p < 0.05; ** p < 0.01; *** p < 0.001. Unadjusted ordinary least-squares tests; 471 species-atlas observations per variable.')+
  theme_minimal(base_size=12)+theme(panel.grid=element_blank(),legend.position='bottom',plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(fig,'vp-paper-lm-slopes.png'),h,width=13,height=3.6,dpi=220,bg='white')
writeLines(c('# VP trends: direct paper-method linear regressions','',
'Source: Antão et al. (2022), Climate change reshuffles northern species within their niches; final paragraph of Joint species distribution modelling, Methods, and Figure 4. https://doi.org/10.1038/s41558-022-01381-x','',
'For each of the eight environmental variables, fit lm(weighted_vp ~ decade_since_1970s), where weighted_vp is the species variance-partition fraction multiplied by its own atlas-specific full-model Tjur R². This is multiplication of the response before regression, not the weights argument to lm. Each fit pools 157 species across three atlases (471 observations, residual df = 469). Model outputs are the 2026-07-27 MinOccs5 CoverageGoodAverage fits; exported VP and Tjur values were checked against VP.rds and MF.rds.','',
'Time is coded 0, 2, 4 for the 1970s, 1990s and 2010s. Slopes are therefore changes in weighted VP per calendar decade based on these representative decades. Coding atlas as 1, 2, 3 instead doubles the slope and standard error but leaves the t statistic, p-value and R² unchanged.','',
'Outputs are ordinary lm slopes, standard errors, 95% confidence intervals, t tests, unadjusted p-values and regression R². No species term, bootstrap, paired test or multiple-testing adjustment is added. This project has one modelled taxonomic assemblage and region, so there is one trend per environmental variable; the paper’s comparisons among independently fitted taxonomic groups and zones are not applicable to these model outputs.','',
'This reproduces the requested point-estimate regression. Its ordinary errors treat species-atlas rows as independent; HMSC posterior uncertainty and repeated-species dependence are not incorporated.','',
'Run Rscript scripts/4-data-visualisation/atlas-vp-paper-lm.R from the project root. CSVs and fitted lm objects are stored here; PNG figures are in manuscript-outputs/figures/exploratory/vp-paper-lm/.','',
'| Variable | Slope / decade | 95% CI | p |','|---|---:|---:|---:|',
 vapply(seq_len(nrow(results)),function(i){r<-results[i,];sprintf('| %s | %+.6f | [%+.6f, %+.6f] | %.5g |',r$label,r$slope_per_decade,r$ci_low,r$ci_high,r$p)},character(1))),file.path(out,'README.md'))
print(results[,c('label','slope_per_decade','ci_low','ci_high','p','stars')],row.names=FALSE)
