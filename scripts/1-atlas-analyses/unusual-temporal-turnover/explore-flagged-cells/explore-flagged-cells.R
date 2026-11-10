# Run from repository root. Read the completed null run; never rerun/mutate it.
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
base <- 'scripts/1-atlas-analyses/unusual-temporal-turnover'
out <- file.path(base,'explore-flagged-cells')
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'data'),showWarnings=FALSE)
# Pin the audited run so the snapshot-specific narrative stays reproducible.
run <- file.path(base,'simulation-results/independent-change/20261104T125100_nsim-999999_seed-20261104')
stopifnot(readLines(file.path(run,'run-status.txt'))[1]=='complete')
e <- new.env(); load(file.path(run,'simulation-results.RData'),e)
stopifnot(all(as.numeric(file.info(e$metadata$input_paths)$mtime)==e$metadata$input_mtime))
env <- read.csv(e$metadata$input_paths[2],check.names=FALSE)
stopifnot(!anyDuplicated(env$survey),length(e$sites)==1698L,
          all(vapply(e$matrices,ncol,integer(1))==201L))
env <- env[match(unlist(lapply(1:3,function(a)paste0(e$sites,'_',a))),env$survey),]
stopifnot(nrow(env)==3*length(e$sites),!anyNA(env$survey))
# All original classes retained in exported data; aggregate forest is additional.
env$forest <- rowSums(env[paste0('pct_class_',40:45)])
env$land <- 100-env$pct_class_0-env$pct_class_77
env$land_diversity <- 1-rowSums((as.matrix(env[paste0('pct_class_',c(11,22:24,33,40:45,55,66))])/env$land)^2)
cover <- c(forest='Total forest',pct_class_41='Needleleaf evergreen forest',
  pct_class_44='Broadleaf deciduous forest',pct_class_45='Mixed forest',
  pct_class_22='Annual crops',pct_class_33='Pasture/rangeland',
  pct_class_55='Unmanaged grass/shrub',pct_class_11='Urban',pct_class_77='Freshwater')
for(v in setdiff(names(cover),'pct_class_77'))env[[paste0(v,'_terrestrial')]]<-100*env[[v]]/env$land
metrics <- c(cover,land='Land area',hh='All-class diversity',land_diversity='Terrestrial diversity')
terr <- paste0(setdiff(names(cover),'pct_class_77'),'_terrestrial')
names_terr <- setNames(paste0(unname(cover[setdiff(names(cover),'pct_class_77')]),' / land'),terr)
metrics <- c(metrics,names_terr)
all_env <- c(paste0('pct_class_',c(0,11,22:24,33,40:45,55,66,77)),
             'forest','land','hh','land_diversity',terr)
tab <- e$site_summary
tab$flag <- tab$q_holm<0.05 & tab$excess>0
tab$group <- ifelse(tab$flag,'Supported excess','Other cells')
tab$label <- factor(tab$pair,levels=c('1-2','2-3','1-3'),
                    labels=c('Atlas 1 → 2','Atlas 2 → 3','Atlas 1 → 3'))
for(a in 1:3)tab[[paste0('richness_atlas',a)]]<-rowSums(e$matrices[[a]])[match(tab$site,e$sites)]
for(n in c('start_richness','end_richness','gains','losses','persistent','union_species'))tab[[n]]<-NA_real_
for(v in all_env)for(s in c('start','end','change'))tab[[paste(s,v,sep='_')]]<-NA_real_
species <- list()
for(pair in names(e$pairs)) {
  ab<-e$pairs[[pair]]; ii<-which(tab$pair==pair); ids<-match(tab$site[ii],e$sites)
  x<-e$matrices[[ab[1]]][ids,,drop=FALSE]; y<-e$matrices[[ab[2]]][ids,,drop=FALSE]
  tab$start_richness[ii]<-rowSums(x); tab$end_richness[ii]<-rowSums(y)
  tab$gains[ii]<-rowSums(x==0 & y==1); tab$losses[ii]<-rowSums(x==1 & y==0)
  tab$persistent[ii]<-rowSums(x==1 & y==1); tab$union_species[ii]<-rowSums(x+y>0)
  for(v in all_env) {
    a<-env[[v]][match(paste0(tab$site[ii],'_',ab[1]),env$survey)]
    b<-env[[v]][match(paste0(tab$site[ii],'_',ab[2]),env$survey)]
    tab[[paste0('start_',v)]][ii]<-a; tab[[paste0('end_',v)]][ii]<-b
    tab[[paste0('change_',v)]][ii]<-b-a
  }
  for(j in which(tab$flag[ii])) {
    changed<-which(x[j,]!=y[j,]); species[[length(species)+1L]]<-data.frame(
      pair=pair,site=tab$site[ii[j]],species=colnames(x)[changed],
      transition=ifelse(y[j,changed]==1,'Gain','Loss'))
  }
}
tab$richness_change<-tab$end_richness-tab$start_richness
stopifnot(max(abs(tab$observed_tbi-(tab$gains+tab$losses)/pmax(1,tab$union_species)))<1e-12,
          all(tab$richness_change==tab$gains-tab$losses),!anyNA(tab))
# Coverage is categorical; Atlas 3 effort is unavailable, not measured as good.
for(a in 1:2) {
  cov<-readxl::read_excel(e$metadata$input_paths[4],sheet=paste('kvadrater atlas',a,sep=''))
  field<-if(a==1)'UNDS_KODE_' else 'UNDS_KOD_1'
  tab[[paste0('coverage_atlas',a)]]<-c('Good','Average')[match(cov[[field]][match(tab$site,cov$KVADRATKOD)],c(1,2))]
}
write.csv(tab,file.path(out,'data/all-cell-pair-profiles.csv'),row.names=FALSE)
write.csv(do.call(rbind,species),file.path(out,'data/flagged-cell-species-changes.csv'),row.names=FALSE)
focal <- c(start_richness='Starting richness',end_richness='Ending richness',
           richness_change='Richness change',gains='Species gains',losses='Species losses',
           persistent='Persistent species')
vars <- c(focal,setNames(paste0('Starting ',metrics),paste0('start_',names(metrics))),
                    setNames(paste0('Change in ',metrics),paste0('change_',names(metrics))))
# Midrank reference percentile: tied values receive the midpoint of their mass.
percentile <- function(x,ref) {
  x<-round(x,8); ref<-round(ref,8) # numerical extraction noise is not a rank difference
  100*(sum(ref<x)+0.5*sum(ref==x))/length(ref)
}
distance_km <- function(lon,lat,xlon,xlat) {
  r<-pi/180; a<-sin((xlat-lat)*r/2)^2+cos(lat*r)*cos(xlat*r)*sin((xlon-lon)*r/2)^2
  6371*2*asin(pmin(1,sqrt(a)))
}
profiles<-list(); controls<-list()
for(pair in names(e$pairs)) {
  z<-tab[tab$pair==pair,]; cases<-which(z$flag)
  for(i in cases) {
    d<-distance_km(z$lon[i],z$lat[i],z$lon,z$lat)
    candidates<-which(!z$flag & d<=75)
    # Guarantee a reasonable candidate pool on islands/coasts; report expansion.
    expanded<-length(candidates)<40
    if(expanded)candidates<-head(which(!z$flag)[order(d[!z$flag])],100)
    for(v in names(vars)) {
      # Richness contrasts: geography/land only. Habitat baselines also match richness.
      # Habitat changes additionally match the same habitat's starting value.
      covars<-'start_land'
      if(!v %in% c('start_richness','end_richness','richness_change','gains','losses','persistent'))
        covars<-c(covars,'start_richness')
      if(grepl('^change_',v))covars<-unique(c(covars,sub('^change_','start_',v)))
      score<-(d[candidates]/75)^2
      for(cv in covars) {
        s<-sd(z[[cv]])
        if(is.finite(s)&&s>1e-10)score<-score+((z[[cv]][candidates]-z[[cv]][i])/s)^2
      }
      chosen<-candidates[order(score,z$site[candidates])][1:20]
      ref<-z[[v]][!z$flag]; local<-z[[v]][chosen]
      profiles[[length(profiles)+1L]]<-data.frame(pair=pair,site=z$site[i],variable=v,
        description=unname(vars[v]),value=z[[v]][i],global_percentile=percentile(z[[v]][i],ref),
        reference_median=median(ref),local_percentile=percentile(z[[v]][i],local),
        local_median=median(local),local_difference=z[[v]][i]-median(local),
        local_n=20,max_distance_km=max(d[chosen]),expanded_pool=expanded)
      controls[[length(controls)+1L]]<-data.frame(pair=pair,site=z$site[i],variable=v,
        control_site=z$site[chosen],distance_km=d[chosen],score=score[match(chosen,candidates)])
    }
  }
}
profiles<-do.call(rbind,profiles); controls<-do.call(rbind,controls)
write.csv(profiles,file.path(out,'data/flagged-cell-percentiles.csv'),row.names=FALSE)
write.csv(controls,file.path(out,'data/local-reference-cells.csv'),row.names=FALSE)
summary<-do.call(rbind,lapply(split(profiles,interaction(profiles$pair,profiles$variable,drop=TRUE)),function(z)
  data.frame(pair=z$pair[1],variable=z$variable[1],description=z$description[1],n=nrow(z),
    flagged_median=median(z$value),other_median=z$reference_median[1],
    median_percentile=median(z$global_percentile),low_10=sum(z$global_percentile<=10),
    high_90=sum(z$global_percentile>=90),median_local_percentile=median(z$local_percentile),
    median_local_difference=median(z$local_difference))))
write.csv(summary,file.path(out,'data/group-descriptive-summary.csv'),row.names=FALSE)
# Continuous associations: descriptive Spearman correlations, without p-values.
assoc<-do.call(rbind,lapply(names(e$pairs),function(pair)do.call(rbind,lapply(names(vars),function(v)
  do.call(rbind,lapply(c('All cells','Other cells'),function(group) {
    z<-tab[tab$pair==pair & (group=='All cells' | !tab$flag),]
    rho<-if(sd(round(z[[v]],8))>0)cor(z$excess,round(z[[v]],8),method='spearman')else NA_real_
    data.frame(pair=pair,variable=v,description=unname(vars[v]),group=group,n=nrow(z),rho=rho)
  }))))))
write.csv(assoc,file.path(out,'data/continuous-associations.csv'),row.names=FALSE)
sens<-do.call(rbind,lapply(names(e$pairs),function(pair) {
  z<-tab[tab$pair==pair,]
  do.call(rbind,lapply(c('Holm < 0.01','Holm < 0.05','Raw p < 0.05','Top 5% excess'),function(rule) {
    hit<-switch(rule,'Holm < 0.01'=z$q_holm<0.01 & z$excess>0,
      'Holm < 0.05'=z$flag,'Raw p < 0.05'=z$p_two_sided<0.05 & z$excess>0,
      'Top 5% excess'=z$excess>=quantile(z$excess,.95))
    data.frame(pair=pair,rule=rule,n=sum(hit),start_richness=median(z$start_richness[hit]),
      richness_change=median(z$richness_change[hit]),forest_change_pp=median(z$change_forest[hit]),
      crops_change_pp=median(z$change_pct_class_22[hit]))
  }))
}))
write.csv(sens,file.path(out,'data/selection-sensitivity.csv'),row.names=FALSE)
effort<-do.call(rbind,lapply(names(e$pairs),function(pair)do.call(rbind,lapply(c(FALSE,TRUE),function(hit) {
  z<-tab[tab$pair==pair & tab$flag==hit,]
  data.frame(pair=pair,flagged=hit,n=nrow(z),average_atlas1=sum(z$coverage_atlas1=='Average'),
    average_atlas2=sum(z$coverage_atlas2=='Average'),
    good_to_average=sum(z$coverage_atlas1=='Good' & z$coverage_atlas2=='Average'))
}))))
write.csv(effort,file.path(out,'data/survey-coverage-summary.csv'),row.names=FALSE)
recovery<-do.call(rbind,lapply(tab$site[tab$pair=='1-2' & tab$flag],function(site) {
  j<-match(site,e$sites); x<-e$matrices[[1]][j,]; y<-e$matrices[[2]][j,]; z<-e$matrices[[3]][j,]
  lost<-x==1 & y==0
  data.frame(site=site,lost_in_atlas2=sum(lost),recorded_again_atlas3=sum(lost & z==1),
    fraction_recorded_again=sum(lost & z==1)/sum(lost))
}))
write.csv(recovery,file.path(out,'data/atlas2-losses-recorded-again.csv'),row.names=FALSE)
theme_set(theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face='bold'),strip.text=element_text(face='bold'),
  plot.caption=element_text(hjust=0),legend.position='bottom'))
save_png<-function(name,p,w,h)ggsave(file.path(out,'figures',paste0(name,'.png')),p,
                                    width=w,height=h,dpi=180,bg='white',limitsize=FALSE)
colours<-c('Other cells'='#B9B9B9','Supported excess'='#D95F0E')
# All-cell richness distribution by atlas; individual flagged trajectories for each interval.
rich<-do.call(rbind,lapply(1:3,function(a)data.frame(site=e$sites,atlas=a,richness=rowSums(e$matrices[[a]]))))
bands<-do.call(rbind,lapply(split(rich,rich$atlas),function(z)data.frame(atlas=z$atlas[1],
   lo=quantile(z$richness,.1),hi=quantile(z$richness,.9),med=median(z$richness))))
traj<-merge(tab[tab$flag,c('pair','site','label')],rich,by='site')
trends<-do.call(rbind,lapply(levels(tab$label),function(l){z<-bands;z$label<-l;z}))
trends$label<-factor(trends$label,levels=levels(tab$label))
traj$label<-factor(traj$label,levels=levels(tab$label))
p<-ggplot()+geom_ribbon(data=trends,aes(atlas,ymin=lo,ymax=hi),fill='grey88')+
  geom_line(data=trends,aes(atlas,med),colour='grey35',linetype=2)+
  geom_line(data=traj,aes(atlas,richness,group=site),colour='#D95F0E',alpha=.55)+
  geom_point(data=traj,aes(atlas,richness),colour='#D95F0E',size=1.4)+
  facet_wrap(~label,nrow=1)+scale_x_continuous(breaks=1:3)+labs(
    title='Richness trajectories of cells with supported excess turnover',x='Atlas',y='Recorded species richness',
    subtitle='Each orange line is one flagged cell; panels identify the interval in which it was flagged',
    caption='Grey band: 10th–90th percentiles across all 1,698 cells; dashed line: median. All 201 species.')
save_png('01-richness-trajectories',p,13,5.5)
p<-ggplot(tab,aes(losses,gains))+geom_abline(slope=1,intercept=0,colour='grey55',linetype=2)+
  geom_point(data=tab[!tab$flag,],colour='grey70',alpha=.3,size=1)+
  geom_point(data=tab[tab$flag,],colour='#D95F0E',size=2)+
  facet_wrap(~label,nrow=1)+coord_equal()+labs(title='Do flagged cells lose species, gain species, or both?',
    subtitle='Orange: supported excess turnover; grey: other included cells',
    x='Species losses',y='Species gains',caption='Above the dashed line: richness increases; below: richness decreases. Counts describe recorded occurrence changes.')
save_png('02-gains-versus-losses',p,13,5.5)
heat<-function(variables,title,name) {
  z<-profiles[profiles$variable %in% variables,]
  z$pair<-factor(z$pair,levels=names(e$pairs))
  z$variable<-factor(z$variable,levels=rev(variables),labels=unname(vars[rev(variables)]))
  z$site<-factor(z$site,levels=unique(tab$site[tab$flag]))
  p<-ggplot(z,aes(site,variable,fill=global_percentile))+geom_tile(colour='white',linewidth=.2)+
    geom_point(data=z[z$global_percentile<=10 | z$global_percentile>=90,],shape=4,size=1.5)+
    facet_grid(~pair,scales='free_x',space='free_x')+
    scale_fill_gradient2(low='#2C7FB8',mid='#FFF7E6',high='#D95F0E',midpoint=50,limits=c(0,100),name='Percentile')+
    labs(title=title,x=NULL,y=NULL,caption='Reference: other cells in the same interval. × marks bottom/top 10%. Ties use midrank percentiles; changes are signed.')+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=8),panel.grid=element_blank())
  save_png(name,p,15,5+length(variables)*.13)
}
heat(names(focal),'Richness and turnover components: where do flagged cells sit?','03-richness-percentiles')
heat(paste0('start_',names(metrics)[1:12]),'Starting habitat: national reference percentiles','04-starting-habitat-percentiles')
heat(paste0('change_',names(metrics)[1:12]),'Habitat change: national reference percentiles','05-habitat-change-percentiles')
heat(c(paste0('start_',terr),paste0('change_',terr)),'Terrestrial-area cover: starting values and changes','06-terrestrial-cover-percentiles')
# Comparison of national/local percentiles: individual cells remain visible.
show<-c('start_richness','richness_change','start_forest','change_forest','start_pct_class_22',
         'change_pct_class_22','start_pct_class_11','change_pct_class_11')
z<-profiles[profiles$variable %in% show,]; z$description<-factor(z$description,levels=unname(vars[show]))
p<-ggplot(z,aes(global_percentile,local_percentile,colour=pair))+geom_abline(slope=1,intercept=0,colour='grey70')+
  geom_hline(yintercept=c(10,90),linetype=3,colour='grey70')+
  geom_vline(xintercept=c(10,90),linetype=3,colour='grey70')+
  geom_point(size=2,alpha=.7)+facet_wrap(~description,ncol=4)+
  scale_colour_manual(values=c('1-2'='#D95F0E','2-3'='#0072B2','1-3'='#009E73'))+
  scale_x_continuous(limits=c(0,100))+scale_y_continuous(limits=c(0,100))+
  labs(title='Are national outliers also unusual among similar nearby cells?',x='National reference percentile',
    y='Local reference percentile',colour='Atlas pair',
    caption='Each point: one flagged cell. Local references: 20 other cells matched on geography and land area;\nhabitat comparisons also use starting richness, and habitat changes additionally use the same starting cover.')
save_png('07-national-versus-local',p,14,8)
# Denmark map with the established Bornholm inset.
maps<-lapply(levels(tab$label),function(l) {
  z<-tab[tab$label==l,]
  panel<-function(xlim,ylim,border=FALSE)ggplot(z,aes(lon,lat))+
    geom_point(colour='grey80',size=1.1)+geom_point(data=z[z$flag,],
      fill='#D95F0E',shape=21,size=2.5,colour='white',stroke=.2)+
    coord_quickmap(xlim=xlim,ylim=ylim,expand=FALSE)+theme_void()+theme(legend.position='none',
      panel.border=if(border)element_rect(fill=NA,colour='grey35')else element_blank())
  main<-panel(c(7.95,13),c(54.4,57.85)); inset<-panel(c(14.65,15.2),c(54.95,55.32),TRUE)
  cowplot::ggdraw(main)+draw_plot(inset,x=.79,y=.68,width=.18,height=.22)+
    draw_label(paste0(l,' (n = ',sum(z$flag),')'),x=.45,y=.98,size=12,fontface='bold')
})
map_body<-plot_grid(plotlist=maps,nrow=1)
save_png('08-flagged-cell-locations',ggdraw(map_body)+
  draw_label('Orange: Holm-supported excess turnover; grey: other included cells',x=.5,y=.015,size=10),13,5.3)
show<-c('start_richness','change_forest','change_pct_class_22','change_pct_class_11')
z<-do.call(rbind,lapply(show,function(v)data.frame(pair=tab$pair,variable=unname(vars[v]),
  value=tab[[v]],excess=tab$excess,flag=tab$flag)))
p<-ggplot(z,aes(value,excess))+geom_hline(yintercept=0,colour='grey75')+
  geom_point(data=z[!z$flag,],size=.6,alpha=.18,colour='grey45')+
  geom_point(data=z[z$flag,],size=1.5,colour='#D95F0E')+
  facet_wrap(vars(pair,variable),ncol=4,scales='free_x')+
  labs(title='Continuous excess turnover across all included cells',x='Richness (species) or cover change (percentage points)',
    y='Observed TBI − null mean',caption='Orange: Holm-supported excess. Habitat changes use whole-cell area. These descriptive associations are not causal estimates.')
save_png('09-continuous-excess-turnover',p,14,10)
# Coverage contrast is essential to interpreting the Atlas 2 richness trough.
z<-tab[tab$pair=='1-2',]; z$coverage_transition<-paste(z$coverage_atlas1,'→',z$coverage_atlas2)
p<-ggplot(z,aes(coverage_transition,richness_change))+
  geom_hline(yintercept=0,colour='grey65')+
  geom_boxplot(outlier.shape=NA,fill='grey90',colour='grey50',width=.5)+
  geom_point(data=z[z$flag,],colour='#D95F0E',position=position_jitter(width=.08,height=0,seed=7),size=2)+
  labs(title='Atlas 1 → 2: richness change and recorded survey coverage',x='Coverage in Atlas 1 → Atlas 2',
    y='Recorded richness change (species)',caption='Boxes: all included cells within each coverage category; orange: supported excess turnover.\nCoverage is categorical; this comparison does not estimate or correct detection probability.')
save_png('10-richness-and-survey-coverage',p,10,6)
# Markdown report and a profile for every flagged cell, regenerated from the data.
fmt<-function(x,d=1)formatC(x,digits=d,format='f')
lines<-c('# Exploring cells with supported excess temporal turnover','',
  paste('Generated:',format(Sys.Date()),'| Source run:',basename(run)),
  '', '## Main findings','',
  'All 30 flagged cells lost recorded species richness during their flagged interval. The clearest shared feature is an unusually large recorded loss, rather than a common extreme habitat profile. Richness and transition counts contribute directly to TBI, so these are descriptions of the selected outcome, not independent explanatory tests.',
  '', '| Interval | Flagged cells | Median starting richness (others) | Median ending richness (others) | Median richness change (others) |',
  '|---|---:|---:|---:|---:|')
for(pair in names(e$pairs)) {
  z<-tab[tab$pair==pair,]; a<-z[z$flag,]; b<-z[!z$flag,]
  lines<-c(lines,sprintf('| %s | %d | %s (%s) | %s (%s) | %s (%s) |',pair,nrow(a),
    fmt(median(a$start_richness)),fmt(median(b$start_richness)),
    fmt(median(a$end_richness)),fmt(median(b$end_richness)),
    fmt(median(a$richness_change)),fmt(median(b$richness_change))))
}
lines<-c(lines,'',
  'All 22 Atlas 1→2 cells increase in richness again in Atlas 3. Nineteen finish Atlas 2 in the bottom tenth of the other cells’ richness distribution. Every flagged cell, in every interval, has losses in the top tenth and richness change in the bottom tenth of its interval’s reference distribution.',
  '',sprintf('Among Atlas 1→2 flagged cells, the median fraction of lost species recorded again in Atlas 3 is %s%%. Pooling species–cell losses gives %s%% (%d/%d); these are different denominators.',
    fmt(100*median(recovery$fraction_recorded_again)),fmt(100*sum(recovery$recorded_again_atlas3)/sum(recovery$lost_in_atlas2)),sum(recovery$recorded_again_atlas3),sum(recovery$lost_in_atlas2)),
  '', '### Survey coverage changes the interpretation','',
  'Nineteen of 22 Atlas 1→2 flagged cells (86%) have Average Atlas 2 coverage, compared with 484/1,676 other cells (29%). Seventeen of 22 (77%) shift from Good to Average, compared with 297/1,676 other cells (18%). The other three flagged cells retain Good coverage in both atlases. All six Atlas 2→3 flagged cells have Good Atlas 2 coverage; Atlas 3 effort is unavailable in this workbook.',
  '', 'The coverage association, richness rebound, and re-recording of lost species make uneven recording a concrete candidate explanation for the Atlas 2 trough. They do not establish that the changes are survey artefacts: real losses and subsequent recoveries could contribute. The current null does not model detection or effort.',
  '', '### Habitat patterns','',
  'The 22 Atlas 1→2 cells span contrasting landscapes. Median starting total forest cover is 5.9% (other cells: 4.8%); median forest change is zero in both groups. Median annual-crop change is −5.2 percentage points (others: −4.3). These medians do not reveal a shared extreme forest/cropland trajectory; individual exceptions are retained in the heatmaps and cell profiles.',
  '', 'The six Atlas 2→3 cells are relatively forested: median starting forest cover is 25.7%, versus 5.7% elsewhere. Their median starting forest percentile is 81 nationally and 78 among the local references. Only one exceeds the national 90th percentile for total forest. Median forest change is −1.1 percentage points versus zero elsewhere. This is a small-group lead, not evidence that forest change explains their turnover.',
  '', 'The two Atlas 1→3 cells, BC45 and FG85, have contrasting trajectories and should be read individually. BC45 has 33→19→10 species; FG85 has 54→71→41.',
  '', 'Across all cells, the largest absolute descriptive Spearman association between excess turnover and the habitat variables considered here is about 0.15. This is weak evidence of a shared monotonic relationship; it does not rule out nonlinear, regional, or habitat-specific associations.',
  '', '## Figures','',
  '1. [Richness trajectories](figures/01-richness-trajectories.png)',
  '2. [Gains versus losses](figures/02-gains-versus-losses.png)',
  '3. [Richness and transition percentiles](figures/03-richness-percentiles.png)',
  '4. [Starting habitat percentiles](figures/04-starting-habitat-percentiles.png)',
  '5. [Habitat change percentiles](figures/05-habitat-change-percentiles.png)',
  '6. [Terrestrial-area cover sensitivity](figures/06-terrestrial-cover-percentiles.png)',
  '7. [National versus local comparisons](figures/07-national-versus-local.png)',
  '8. [Flagged cell locations, with Bornholm inset](figures/08-flagged-cell-locations.png)',
  '9. [Continuous excess turnover](figures/09-continuous-excess-turnover.png)',
  '10. [Richness change and survey coverage](figures/10-richness-and-survey-coverage.png)',
  '', '## Definitions and reproducibility','',
  'Run from repository root: `Rscript scripts/1-atlas-analyses/unusual-temporal-turnover/explore-flagged-cells/explore-flagged-cells.R`.',
  '', 'The script loads the completed 999,999-draw simulation, checks all four current input modification times against its metadata, and uses its saved occurrence matrices and 1,698 common cells / 201 species. It checks reconstructed TBI and richness changes. It neither reruns nor modifies the simulation. Existing land/coverage inclusion rules are preserved: at least 25% terrestrial area in every atlas; Good/Average effort in Atlases 1 and 2; adequate Atlas 3 effort is an assumption.',
  '', 'Supported excess means positive observed TBI minus null mean and two-sided Holm-adjusted p < 0.05 over all 5,094 cell-by-interval tests. The independent-change null fixes each species’ earlier occupancy and gain/loss totals, but independently reallocates changes without geography, habitat, or effort constraints. No cells are flagged for excess in more than one interval.',
  '', 'TBI = (gains + losses) / species recorded in either atlas. Richness change = gains − losses. A cell’s baseline species composition influences its null expectation; excess is not raw turnover. “Other cells” includes all cells without supported excess, including cells with supported low turnover.',
  '', 'National reference percentiles compare each flagged cell with all other cells in its interval (1,676 / 1,692 / 1,696). Percentile = 100 × (number strictly below + half the number tied) / reference size. Bottom/top tenth use ≤10/≥90. Values are rounded to eight decimals for percentile and correlation ties to avoid numerical extraction noise. Zero-heavy cover variables can have broad tied ranks; an extreme percentile need not imply a large absolute change.',
  '', 'Whole-cell cover percentages include ocean and freshwater in the denominator. Total forest sums HILDA classes 40–45; annual crops means class 22 only. Terrestrial cover is 100 × period-mean cover / (100 − period-mean ocean − period-mean freshwater). Signed changes are later minus earlier, in percentage points. They measure net compositional change, not mapped gross habitat transitions. HILDA periods are 1971–1974, 1993–1996 and 2014–2017.',
  '', 'Existing hh is the mean annual Simpson diversity over all original classes including water. The additional terrestrial diversity is 1 − sum(p²) over terrestrial classes, calculated from period-mean cover renormalised to land. It is a descriptive sensitivity measure with a different temporal aggregation from hh. Neither is fragmentation.',
  '', 'Local references are exploratory comparisons, not a causal matching design or an exchangeability model. For each flagged cell and variable, consider other cells within 75 km. If fewer than 40 are available, use the nearest 100 other cells. Select 20 by smallest score: (distance/75 km)² plus squared differences in covariates divided by their whole-grid SDs. Richness/transition comparisons match on starting land area; habitat comparisons also match on starting richness; habitat changes additionally match on the same starting habitat variable. There are no strict covariate calipers, matching is variable-specific, and controls may be reused. Local land-area ranks are constrained by matching on land itself and should not be read as independent evidence. Exact controls, scores and distances are exported.',
  '', 'No new significance tests are performed. Spatial dependence, outcome-based selection, overlapping atlas intervals, imperfect matching and n=2/n=6 groups limit inference. Spearman correlations are descriptive; differences in selected cells do not independently validate the selection criterion.',
  '', 'Selection sensitivity compares Holm <0.05, raw p <0.05 and the top 5% of excess. Holm <0.01 selects zero because the smallest attainable initial Holm value with 999,999 draws and doubled inclusive Monte Carlo tails is 5,094 × 2/1,000,000 = 0.010188. This stricter cutoff is limited by Monte Carlo resolution, rather than demonstrating absence of stronger support.',
  '', '## Exported data','',
  '- `data/all-cell-pair-profiles.csv`: all 5,094 rows, raw bird counts and every original cover class at start/end/change, plus aggregate/sensitivity variables and coverage.',
  '- `data/flagged-cell-percentiles.csv`: flagged values, national/local percentiles and reference medians.',
  '- `data/local-reference-cells.csv`: exact controls and distance/score for each variable and case.',
  '- `data/group-descriptive-summary.csv`: medians and extreme-percentile counts by interval/variable.',
  '- `data/flagged-cell-species-changes.csv`: each changed species in each flagged cell.',
  '- `data/atlas2-losses-recorded-again.csv`: species lost in Atlas 2 and re-recorded in Atlas 3, for Atlas 1→2 flagged cells.',
  '- `data/survey-coverage-summary.csv`, `data/continuous-associations.csv`, `data/selection-sensitivity.csv`.',
  '', '## Individual cell profiles','',
  'Habitat entries below show the most extreme national reference percentiles among starting cover and cover changes (whole-cell metrics); percentile extremeness does not measure effect magnitude. Full profiles, including terrestrial denominators and local comparisons, are in the CSVs.')
sp<-do.call(rbind,species)
for(i in which(tab$flag)) {
  z<-tab[i,]; q<-profiles[profiles$pair==z$pair & profiles$site==z$site &
      profiles$variable %in% c(paste0('start_',names(cover)),paste0('change_',names(cover))),]
  q<-head(q[order(-abs(q$global_percentile-50)),],3)
  changes<-sp[sp$pair==z$pair & sp$site==z$site,]
  lines<-c(lines,'',paste0('### ',z$site,' | Atlas ',z$pair),'',
    sprintf('- Richness: %d → %d → %d. Flagged interval: %d gains, %d losses, %d persistent species.',
      z$richness_atlas1,z$richness_atlas2,z$richness_atlas3,z$gains,z$losses,z$persistent),
    sprintf('- TBI %.3f; null mean %.3f; excess +%.3f; Holm p %.6f.',z$observed_tbi,z$null_mean,z$excess,z$q_holm),
    paste0('- Coverage: ',z$coverage_atlas1,' → ',z$coverage_atlas2,' → unavailable.'),
    paste0('- ',q$description,': ',fmt(q$value,2),ifelse(grepl('^change_',q$variable),' percentage points','%'),'; reference percentile ',fmt(q$global_percentile),
      '; local percentile ',fmt(q$local_percentile),'.'),
    paste0('- Lost species: ',paste(gsub('_',' ',changes$species[changes$transition=='Loss']),collapse=', '),'.'),
    paste0('- Gained species: ',if(any(changes$transition=='Gain'))paste(gsub('_',' ',changes$species[changes$transition=='Gain']),collapse=', ')else 'none','.'))
}
writeLines(lines,file.path(out,'README.md'))
message('Finished: ',out,'; ',sum(tab$flag),' flagged cell-pairs; ',length(unique(tab$site[tab$flag])),' distinct cells.')
