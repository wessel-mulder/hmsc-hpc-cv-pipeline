# Temporal VP analysis adapted from Antao et al. (2022), Fig. 4.
# Run from the project root: Rscript scripts/4-data-visualisation/atlas-vp-temporal-change.R
# Uses existing point estimates; no HMSC refitting. Plots are PNG only.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
model_date <- '2026-07-27'
base <- 'intermediate-outputs/hmsc-output'
out <- 'intermediate-outputs/vp-temporal-change'
fig <- 'manuscript-outputs/figures/exploratory/vp-temporal-change'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
variables <- c(tmean_breeding='Temperature',prec_breeding='Precipitation',hh='Land-use heterogeneity',perc_urban='Urban cover',perc_cropland='Cropland cover',perc_pasture='Pasture cover',perc_forest='Forest cover',perc_grass_shrub='Grass / shrub cover')
components <- c(names(variables),'Random: site')
periods <- c('1970s','1990s','2010s')
folders <- sort(list.dirs(base,recursive=FALSE))
folders <- folders[grepl(paste0('^',model_date,'.*_Atlas[123]_MinOccs5_CoverageGoodAverage$'),basename(folders))]
stopifnot(length(folders)==3)
provenance <- list(); matrices <- list()
for (folder in folders) {
  atlas <- as.integer(sub('.*Atlas([123])_.*','\\1',basename(folder)))
  path <- file.path(folder,'Results',paste0(basename(folder),'parameter_estimates_VP_.csv'))
  a <- as.matrix(read.csv(path,row.names=1,check.names=FALSE))
  vp_path <- file.path(folder,paste0(basename(folder),'VP.rds'))
  vp <- readRDS(vp_path)$vals
  mf_path <- file.path(folder,paste0(basename(folder),'MF.rds'))
  mf <- readRDS(mf_path)
  stopifnot(setequal(rownames(a),c(components,'TjurR2')), !anyDuplicated(colnames(a)), all(is.finite(a)),
            max(abs(colSums(a[components,,drop=FALSE])-1))<1e-10,
            identical(colnames(a),colnames(vp)), max(abs(a[rownames(vp),]-vp))<1e-10,
            max(abs(a['TjurR2',]-mf$TjurR2))<1e-10, all(a>=0),all(a<=1))
  matrices[[atlas]] <- a
  provenance[[atlas]] <- data.frame(atlas=atlas,period=periods[atlas],source=c(path,vp_path,mf_path),md5=unname(tools::md5sum(c(path,vp_path,mf_path))))
}
stopifnot(all(vapply(matrices,function(a)setequal(colnames(a),colnames(matrices[[1]])),logical(1))))
species <- colnames(matrices[[1]])
matrices <- lapply(matrices,function(a)a[,species,drop=FALSE])
n <- length(species)
# Sample whole species trajectories together, preserving pairing across all periods and variables.
set.seed(20260909)
B <- 10000L
boot_index <- replicate(B,sample.int(n,n,replace=TRUE))
boot_mean <- function(y) colMeans(matrix(y[boot_index],nrow=n))
ci <- function(y) unname(quantile(boot_mean(y),c(.025,.975)))
long <- bind_rows(lapply(1:3,function(i){
  a<-matrices[[i]]
  bind_rows(lapply(components,function(v)data.frame(species=species,variable=v,atlas=i,period=periods[i],atlas_step=i-1,relative_vp=as.numeric(a[v,]),tjur_r2=as.numeric(a['TjurR2',]),weighted_vp=as.numeric(a[v,]*a['TjurR2',]))))
}))
write.csv(long,file.path(out,'species-atlas-vp-values.csv'),row.names=FALSE)
write.csv(bind_rows(provenance),file.path(out,'input-provenance.csv'),row.names=FALSE)
metrics <- c('weighted_vp','relative_vp')
summary_rows <- list(); atlas_rows <- list(); species_rows <- list(); contrast_rows <- list()
for(metric in metrics) for(v in components) {
  z <- sapply(1:3,function(i) {
    a<-matrices[[i]]; if(metric=='weighted_vp')as.numeric(a[v,]*a['TjurR2',]) else as.numeric(a[v,])
  })
  # For x = 0,1,2, the OLS slope is exactly (y3-y1)/2 for each species.
  slopes <- (z[,3]-z[,1])/2
  s_ci <- ci(slopes); tst <- t.test(slopes)
  pooled <- lm(as.vector(z) ~ rep(0:2,each=n))
  stopifnot(abs(coef(pooled)[2]-mean(slopes))<1e-12)
  summary_rows[[length(summary_rows)+1]] <- data.frame(metric=metric,variable=v,n_species=n,
    mean_slope=mean(slopes),bootstrap_low=s_ci[1],bootstrap_high=s_ci[2],paired_t_p=tst$p.value,
    pooled_ols_slope=unname(coef(pooled)[2]),pooled_ols_p=coef(summary(pooled))[2,4],
    mean_atlas1=mean(z[,1]),mean_atlas2=mean(z[,2]),mean_atlas3=mean(z[,3]),
    species_fraction_increasing=mean(slopes>0))
  species_rows[[length(species_rows)+1]] <- data.frame(species=species,variable=v,metric=metric,slope=slopes,
    delta_1_to_2=z[,2]-z[,1],delta_2_to_3=z[,3]-z[,2],delta_1_to_3=z[,3]-z[,1],
    monotonic_increase=z[,2]>z[,1]&z[,3]>z[,2],monotonic_decrease=z[,2]<z[,1]&z[,3]<z[,2])
  for(i in 1:3){a_ci<-ci(z[,i]);atlas_rows[[length(atlas_rows)+1]]<-data.frame(metric=metric,variable=v,atlas=i,mean=mean(z[,i]),median=median(z[,i]),bootstrap_low=a_ci[1],bootstrap_high=a_ci[2])}
  for(pair in list(c(1,2),c(2,3),c(1,3))){delta<-z[,pair[2]]-z[,pair[1]];d_ci<-ci(delta)
    contrast_rows[[length(contrast_rows)+1]]<-data.frame(metric=metric,variable=v,contrast=paste(pair,collapse='_to_'),mean_change=mean(delta),bootstrap_low=d_ci[1],bootstrap_high=d_ci[2],paired_t_p=t.test(delta)$p.value)}
}
slopes <- bind_rows(summary_rows) |> group_by(metric) |> mutate(paired_t_q=replace(rep(NA_real_,n()),variable %in% names(variables),p.adjust(paired_t_p[variable %in% names(variables)],'BH'))) |> ungroup()
contrasts <- bind_rows(contrast_rows) |> group_by(metric) |> mutate(paired_t_q=replace(rep(NA_real_,n()),variable %in% names(variables),p.adjust(paired_t_p[variable %in% names(variables)],'BH'))) |> ungroup()
atlas_summary <- bind_rows(atlas_rows)
write.csv(slopes,file.path(out,'overall-trends.csv'),row.names=FALSE)
write.csv(contrasts,file.path(out,'paired-period-changes.csv'),row.names=FALSE)
write.csv(atlas_summary,file.path(out,'atlas-summary.csv'),row.names=FALSE)
write.csv(bind_rows(species_rows),file.path(out,'species-trends.csv'),row.names=FALSE)
# Exact symmetric product decomposition of the first-to-last weighted change.
# d(VP * R2) = mean(R2 endpoints)*dVP + mean(VP endpoints)*dR2.
decomposition <- bind_rows(lapply(components,function(v){
  a<-matrices[[1]]; b<-matrices[[3]]
  vp_part<-(a['TjurR2',]+b['TjurR2',])/2*(b[v,]-a[v,])
  fit_part<-(a[v,]+b[v,])/2*(b['TjurR2',]-a['TjurR2',])
  total<-b[v,]*b['TjurR2',]-a[v,]*a['TjurR2',]
  stopifnot(max(abs(vp_part+fit_part-total))<1e-12)
  data.frame(variable=v,weighted_change=mean(total),vp_fraction_contribution=mean(vp_part),tjur_change_contribution=mean(fit_part))
}))
write.csv(decomposition,file.path(out,'weighted-change-decomposition.csv'),row.names=FALSE)
# Structural checks: changes sum and the total weighted component change equals R2 change.
stopifnot(max(abs((contrasts |> filter(contrast=='1_to_3') |> pull(mean_change))-2*slopes$mean_slope))<1e-10)
stopifnot(abs(sum(decomposition$weighted_change)-mean(matrices[[3]]['TjurR2',]-matrices[[1]]['TjurR2',]))<1e-10)
label_var <- function(v)factor(unname(variables[v]),levels=unname(variables))
colors <- c(Increase='#B2182B',Decrease='#2166AC')
main <- slopes |> filter(metric=='weighted_vp',variable %in% names(variables)) |> mutate(label=label_var(variable),direction=if_else(mean_slope>0,'Increase','Decrease'),annotation=sprintf('%+.2f pp / atlas step\nBH q = %s',100*mean_slope,format.pval(paired_t_q,digits=2,eps=.001)))
points <- long |> filter(variable %in% names(variables)) |> mutate(label=label_var(variable))
means <- atlas_summary |> filter(metric=='weighted_vp',variable %in% names(variables)) |> mutate(label=label_var(variable))
fit_lines <- bind_rows(lapply(seq_len(nrow(main)),function(i){d<-main[i,];data.frame(label=d$label,atlas=c(1,3),value=mean(c(d$mean_atlas1,d$mean_atlas2,d$mean_atlas3))+c(-1,1)*d$mean_slope,direction=d$direction)}))
p <- ggplot(points,aes(atlas,100*weighted_vp)) +
  geom_violin(aes(group=atlas),fill='#E7E9EC',colour='#BBC0C6',width=.65,linewidth=.3,trim=TRUE)+
  geom_point(position=position_jitter(width=.12,height=0,seed=17),size=.55,alpha=.35,colour='#40464D')+
  geom_line(data=fit_lines,aes(atlas,100*value,colour=direction),inherit.aes=FALSE,linewidth=1)+
  geom_errorbar(data=means,aes(x=atlas,ymin=100*bootstrap_low,ymax=100*bootstrap_high),inherit.aes=FALSE,width=.08,linewidth=.65)+
  geom_point(data=means,aes(atlas,100*mean),inherit.aes=FALSE,size=2.1,shape=21,fill='white',stroke=.8)+
  geom_text(data=main,aes(x=1,y=Inf,label=annotation),inherit.aes=FALSE,hjust=0,vjust=1.3,size=3.2)+
  facet_wrap(~label,ncol=4,scales='free_y')+scale_x_continuous(breaks=1:3,labels=periods,limits=c(.65,3.35))+
  scale_y_continuous(expand=expansion(mult=c(.02,.19)))+scale_colour_manual(values=colors,guide='none')+
  labs(title='Environmental VP through the three atlases',subtitle=paste0(n,' matched species | VP fraction × species-level Tjur’s R² | trend adapted from Antão et al. (2022)'),x=NULL,y='Tjur-weighted VP (× 100)',caption='Dots: species; white circles: species means with 95% paired-species bootstrap intervals. Red / blue: increasing / decreasing fitted trend.\nPanel y-scales differ. q: paired-species slope test, BH-adjusted across 8 predictors; does not include HMSC posterior uncertainty.')+
  theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(face='bold',size=19),plot.subtitle=element_text(size=11),plot.caption=element_text(hjust=0,size=9),panel.spacing=grid::unit(1.3,'lines'))
ggsave(file.path(fig,'atlas-vp-temporal-trajectories.png'),p,width=15,height=9,dpi=220,bg='white')
forest <- slopes |> filter(variable %in% names(variables)) |> mutate(label=factor(unname(variables[variable]),levels=rev(unname(variables))),metric_label=if_else(metric=='weighted_vp','Tjur-weighted VP','Relative VP (unweighted)'),direction=if_else(mean_slope>0,'Increase','Decrease'))
p2 <- ggplot(forest,aes(100*mean_slope,label,colour=direction)) + geom_vline(xintercept=0,linetype='dashed',colour='grey55')+
  geom_segment(aes(x=100*bootstrap_low,xend=100*bootstrap_high,yend=label),linewidth=.85)+geom_point(size=2.8)+
  facet_wrap(~metric_label,scales='free_x')+scale_colour_manual(values=colors,guide='none')+
  labs(title='Which environmental contributions changed?',subtitle='Mean within-species slope and 95% paired-species bootstrap interval',x='Change per atlas step (percentage points)',y=NULL,caption='One step = 1970s to 1990s or 1990s to 2010s; not one calendar decade. The two panels have different x-scales.\nIntervals describe variation among 157 matched species, conditional on fitted VP and Tjur estimates; 10,000 bootstrap resamples.')+
  theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(face='bold',size=18),plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(fig,'atlas-vp-temporal-slopes.png'),p2,width=12,height=6.5,dpi=220,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'session-info.txt'))
print(slopes |> filter(metric=='weighted_vp',variable %in% names(variables)) |> select(variable,mean_atlas1,mean_atlas2,mean_atlas3,mean_slope,bootstrap_low,bootstrap_high,paired_t_q))
message('Wrote analysis to ',out,' and PNGs to ',fig)
# Pairwise changes reveal reversals that a three-period linear slope cannot describe.
heat <- contrasts |> filter(metric=='weighted_vp',variable %in% names(variables)) |>
  mutate(variable_label=label_var(variable),comparison=factor(contrast,levels=c('1_to_3','2_to_3','1_to_2'),labels=c('1970s to 2010s','1990s to 2010s','1970s to 1990s')),
         star=case_when(paired_t_q<.001~'***',paired_t_q<.01~'**',paired_t_q<.05~'*',TRUE~''),
         annotation=sprintf('%+.2f%s',100*mean_change,star))
lim <- max(abs(100*heat$mean_change))
p3 <- ggplot(heat,aes(variable_label,comparison,fill=100*mean_change))+geom_tile(colour='white',linewidth=1)+
  geom_text(aes(label=annotation),size=4)+
  scale_fill_gradient2(low='#4393C3',mid='white',high='#D6604D',midpoint=0,limits=c(-lim,lim),name='Change (pp)')+
  scale_x_discrete(labels=c('Temperature','Precipitation','Land-use\nheterogeneity','Urban','Cropland','Pasture','Forest','Grass /\nshrub'))+
  labs(title='Changes between individual atlas periods',subtitle='Mean change in Tjur-weighted VP (percentage points); 157 paired species',x=NULL,y=NULL,caption='Stars: paired-species tests, BH-adjusted across all 24 predictor-by-contrast tests (* q < 0.05; ** q < 0.01; *** q < 0.001).\nThese exploratory tests condition on fitted estimates and do not propagate HMSC posterior uncertainty.')+
  theme_minimal(base_size=12)+theme(panel.grid=element_blank(),axis.text.x=element_text(size=11),plot.title=element_text(face='bold',size=18),plot.caption=element_text(hjust=0,size=9),legend.position='bottom')
ggsave(file.path(fig,'atlas-vp-period-changes.png'),p3,width=13,height=4.7,dpi=220,bg='white')
rows <- main |> mutate(result=if_else(paired_t_q<.05,'Increase supported by paired test','No clear overall trend'))
report_table <- c('| Predictor | 1970s | 1990s | 2010s | Slope per atlas step [95% bootstrap interval] | BH q |',
                  '|---|---:|---:|---:|---:|---:|',
                  vapply(seq_len(nrow(rows)),function(i){r<-rows[i,];sprintf('| %s | %.3f | %.3f | %.3f | %+.3f [%+.3f, %+.3f] | %.4g |',as.character(r$label),100*r$mean_atlas1,100*r$mean_atlas2,100*r$mean_atlas3,100*r$mean_slope,100*r$bootstrap_low,100*r$bootstrap_high,r$paired_t_q)},character(1)))
writeLines(c(
'# Environmental variance partitioning through the three atlases',
'',
'## Reference and adaptation',
'',
'Antão et al. (2022), *Climate change reshuffles northern species within their niches*, Nature Climate Change 12:587-592, Figure 4 and Methods (PDF pages 5 and 8). DOI: https://doi.org/10.1038/s41558-022-01381-x.',
'The paper regresses species-level variance components, weighted by Tjur’s R², against decade within taxon and zone. Its public extraction script explicitly multiplies each species VP by that species Tjur’s R² (and rounds the result); this implementation retains full precision. Source: https://github.com/benweigel/RECChange/blob/main/03%20Model%20parameter%20extracts.R. That extraction script does not contain the Figure 4 trend fits; the regression specification is drawn from the paper.',
'',
'## Project data and response',
'',
'Uses the three 2026-07-27 Atlas1/2/3 MinOccs5 CoverageGoodAverage model outputs under `intermediate-outputs/hmsc-output/`. All 157 species are present in all three fitted species lists, with identical names and order after explicit matching. The eight environmental predictors are analysed individually; the site random effect is retained in diagnostic tables but excluded from environmental multiple-testing families.',
'For species s, predictor v and atlas t, W[s,v,t] = VP[s,v,t] * TjurR2[s,t]. The Tjur estimates are the stored full-model fitted discrimination values, not cross-validation scores. Each CSV VP value was verified against its saved VP.rds, and each Tjur value against MF.rds; component sums equal one. No models are refitted, and no predictors are rescaled.',
'W is a Tjur-weighted VP index, not a uniquely identified causal contribution or an exact decomposition of observed binary-response variance. Unweighted VP is a relative allocation among the model’s environmental and random components.',
'',
'## Temporal model and uncertainty',
'',
'Atlas order is coded 0, 1, 2 for the 1970s, 1990s and 2010s. Slopes are per atlas step (roughly twenty years), not per calendar decade. A slope is fitted within each species and averaged across the matched species. With these equally spaced three periods the slope is exactly (W2010s - W1970s)/2: the middle period does not determine the slope. The trajectories and all three paired contrasts are therefore also reported to expose reversals.',
'For this balanced panel, the average within-species slope equals the coefficient of pooled lm(W ~ atlas_step). The latter slope and its ordinary, unpaired p-value are exported for comparison with the paper-style regression; the ordinary p-value is not used for the plotted inference.',
'Primary intervals use 10,000 bootstrap resamples of whole species trajectories, with the same resampling indices across variables and periods (seed 20260909). Plotted q-values come from two-sided one-sample t-tests of species slopes against zero, Benjamini-Hochberg adjusted across eight environmental predictors. Relative-VP tests form a separate sensitivity family. Pairwise contrast tests are adjusted across 24 comparisons (eight predictors times three contrasts). The site random-effect tests are diagnostic and have no adjusted q-value.',
'Pairing handles repeated measurements of each species. The bootstrap and t-tests still treat species as independent/exchangeable replicates; they do not model phylogenetic dependence or dependence induced by a jointly fitted community model. They describe cross-species variation conditional on stored estimates and do not propagate HMSC posterior uncertainty, model-fit uncertainty or grid sampling uncertainty. Interpret significance as exploratory, not as a posterior probability of temporal change.',
'',
'## Results',
'',
'Atlas means below are W multiplied by 100; slopes and intervals are percentage points per atlas step.',
'',report_table,'',
'Urban, cropland, pasture, forest and grass/shrub cover show positive weighted slopes with q < 0.05. Forest has the largest mean increase. Temperature, precipitation and land-use heterogeneity lack a clear weighted overall slope at this threshold. The same five cover variables increase in unweighted VP; heterogeneity also increases in unweighted VP. Thus the main cover-variable result is not generated solely by changing Tjur weights.',
'',
'## Interpretation and checks',
'',
'Increasing VP indicates greater allocation of fitted variation to that predictor under the model and environmental distribution of that atlas. It does not identify an increasing causal effect of land-use change. VP components share a denominator, environmental correlations can change, and the study-grid coverage differs across atlases. All these can affect comparisons. The spatial random component is included in the exports to inspect allocation changes.',
'The endpoint weighting decomposition is exact: change(VP * R2) = mean(endpoint R2) * change(VP) + mean(endpoint VP) * change(R2), calculated for each species and then averaged. This is an algebraic diagnostic, not a causal decomposition. The script checks that contributions sum to each weighted change, and that changes summed over all environmental and spatial components equal the change in mean Tjur R².',
'',
'## Files and reproduction',
'',
'Run `Rscript scripts/4-data-visualisation/atlas-vp-temporal-change.R` from the repository root.',
'Figures (PNG only) are in `manuscript-outputs/figures/exploratory/vp-temporal-change/`: trajectories, weighted/unweighted slope comparison, and paired period-change heatmap.',
'Analysis tables in this folder: `species-atlas-vp-values.csv`, `species-trends.csv`, `overall-trends.csv`, `atlas-summary.csv`, `paired-period-changes.csv`, and `weighted-change-decomposition.csv`. VP values are stored as fractions; changes, slopes and their intervals use the corresponding unmultiplied scale. Multiply these effect estimates and intervals by 100 to recover plotted units. Counts, indices, p-values and q-values retain their usual units. `input-provenance.csv` records source paths and MD5 checksums; `session-info.txt` records the R environment.'
),file.path(out,'README.md'))
