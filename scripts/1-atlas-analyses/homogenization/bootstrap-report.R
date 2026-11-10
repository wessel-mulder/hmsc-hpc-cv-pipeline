# Read completed draws; this script never runs bootstrap simulations.
source(file.path(out,'bootstrap-functions.R'),local=TRUE)
pointer <- file.path(out,'bootstrap-results/latest-run.txt')
bootstrap_available <- file.exists(pointer)
if(bootstrap_available) {
  run_name<-readLines(pointer,warn=FALSE)
  stopifnot(length(run_name)==1,identical(basename(run_name),run_name))
  bootstrap_run_dir<-file.path(out,'bootstrap-results',run_name)
  if(!file.exists(file.path(bootstrap_run_dir,'run-status.txt')) ||
     !identical(readLines(file.path(bootstrap_run_dir,'run-status.txt'),warn=FALSE),'complete'))
    stop('Latest bootstrap run is incomplete; rerun the explicit bootstrap runner.')
  boot<-new.env(); load(file.path(bootstrap_run_dir,'bootstrap-results.RData'),envir=boot)
  if(!identical(boot$metadata$input_code_hashes,bootstrap_hashes(root,out)))
    stop('Saved bootstrap results are incompatible with current inputs or bootstrap code. Run run-spatial-bootstrap.R again.')
  stopifnot(identical(boot$sites,sites),max(abs(boot$beta-beta))<1e-12)
  bs <- boot$summary_table
  rownames(bs)<-NULL
  interval_data<-rbind(data.frame(bs,method='Percentile',low=bs$percentile_low,high=bs$percentile_high),
    data.frame(bs,method='Basic',low=bs$basic_low,high=bs$basic_high))
  interval_data$width_km<-factor(interval_data$width_km,levels=c(10,20,40,80))
  interval_data$origin<-factor(interval_data$origin,levels=c('unshifted','east','north','both'))
  for(batch_i in 1:2) {
    z<-interval_data[interval_data$batch==batch_i,]
    p<-ggplot(z,aes(origin,bootstrap_mean,colour=method))+
      geom_hline(yintercept=0,colour='grey60',linetype='dashed')+
      geom_hline(data=data.frame(pair=names(boot$observed),value=unname(boot$observed)),
        aes(yintercept=value),colour='#287D8E',linetype='dotted',inherit.aes=FALSE)+
      geom_linerange(aes(ymin=low,ymax=high),position=position_dodge(width=.4))+
      geom_point(position=position_dodge(width=.4),size=1.5)+
      facet_grid(pair~width_km,scales='free_y')+
      scale_colour_manual(values=c(Basic='#C86642',Percentile='#287D8E'))+
      labs(x='Grid origin (half-block shifts)',y='Change in mean spatial Jaccard',colour='95% interval',
        title=paste('Exploratory bootstrap intervals: seed batch',batch_i),
        subtitle='Columns: block width (km); rows: atlas comparison. Points: bootstrap means.',
        caption='Dotted teal: observed change; dashed grey: zero. Separate configurations are not pooled.\nIntervals have not been validated by coverage simulation.')+
      theme_bw(base_size=11)+theme(legend.position='bottom',axis.text.x=element_text(angle=35,hjust=1),plot.caption=element_text(hjust=0))
    ggsave(file.path(bootstrap_run_dir,paste0('intervals-batch-',batch_i,'.png')),p,width=13,height=10,dpi=200,bg='white')
  }
  geometry_summary<-do.call(rbind,lapply(split(boot$block_audit,boot$block_audit$config),function(z)
    data.frame(config=z$config[1],blocks=nrow(z),min_cells=min(z$cells),median_cells=median(z$cells),
      max_cells=max(z$cells),singletons=sum(z$cells==1))))
  rownames(geometry_summary)<-NULL
  replicate_summary<-do.call(rbind,lapply(split(boot$draws,interaction(boot$draws$config,boot$draws$batch)),function(z)
    data.frame(config=z$config[1],batch=z$batch[1],min_total=min(z$total_cells),max_total=max(z$total_cells),
      mean_unique_cells=mean(z$unique_cells),mean_unique_blocks=mean(z$unique_blocks),
      mean_duplicate_pair_fraction=mean(z$duplicate_pair_fraction))))
  stability<-merge(bs[bs$batch==1,],bs[bs$batch==2,],by=c('config','pair'),suffixes=c('_1','_2'))
  stability<-data.frame(config=stability$config,pair=stability$pair,
    se_difference=stability$se_2-stability$se_1,
    percentile_low_difference=stability$percentile_low_2-stability$percentile_low_1,
    percentile_high_difference=stability$percentile_high_2-stability$percentile_high_1)
  regional_summary<-do.call(rbind,lapply(split(boot$draws,interaction(boot$draws$config,boot$draws$batch)),function(z)
    do.call(rbind,lapply(boot$metadata$quadrant_levels,function(q) {
      x<-z[[paste0('share_',q)]]
      data.frame(config=z$config[1],batch=z$batch[1],quadrant=q,
        original=unname(boot$metadata$original_quadrant_share[q]),mean=mean(x),
        low=unname(quantile(x,.025)),high=unname(quantile(x,.975)),omitted_fraction=mean(x==0))
    }))))
  diagnostic_long<-do.call(rbind,lapply(c('total_cells','unique_cells','duplicate_pair_fraction'),function(v)
    data.frame(boot$draws[c('config','width_km','origin','batch')],diagnostic=v,value=boot$draws[[v]])))
  p<-ggplot(diagnostic_long,aes(factor(width_km),value,fill=factor(batch)))+
    geom_boxplot(outlier.shape=NA,position=position_dodge(.8))+
    facet_grid(diagnostic~origin,scales='free_y')+labs(x='Block width (km)',y=NULL,fill='Seed batch',
      title='Resampled cell counts and repeated-history comparisons',
      subtitle='Duplicate-pair fraction counts comparisons between copies of the same original cell')+
    theme_bw(base_size=11)+theme(legend.position='bottom')
  ggsave(file.path(bootstrap_run_dir,'resampling-diagnostics.png'),p,width=13,height=9,dpi=200,bg='white')
  # Explicit origins avoid encoding configuration labels as ecological regions.
  regional_summary$origin<-sub('^[0-9]+km-','',regional_summary$config)
  regional_summary$width_km<-as.numeric(sub('km-.*','',regional_summary$config))
  p<-ggplot(regional_summary,aes(factor(width_km),mean,colour=factor(batch)))+
    geom_hline(aes(yintercept=original),colour='grey70',linetype='dashed')+
    geom_linerange(aes(ymin=low,ymax=high),position=position_dodge(.4))+
    geom_point(position=position_dodge(.4))+facet_grid(quadrant~origin)+
    labs(x='Block width (km)',y='Fraction of resampled cell histories',colour='Seed batch',
      title='Spatial representation in bootstrap replicates',
      subtitle='Mean and central 95% of shares; dashed lines: original shares. Quadrants are fixed audit regions.')+
    theme_bw(base_size=11)+theme(legend.position='bottom')
  ggsave(file.path(bootstrap_run_dir,'spatial-representation.png'),p,width=13,height=10,dpi=200,bg='white')
  influence_long<-do.call(rbind,lapply(c('12','23','13'),function(v)
    data.frame(boot$influences[c('config','block','cells_removed')],pair=paste(strsplit(v,'')[[1]],collapse='-'),
      change=boot$influences[[paste0('delta_',v)]]-boot$observed[paste(strsplit(v,'')[[1]],collapse='-')])))
  p<-ggplot(influence_long,aes(cells_removed,change)) +geom_hline(yintercept=0,colour='grey70')+
    geom_point(alpha=.45,size=.8)+facet_grid(pair~config,scales='free_x')+
    labs(x='Cells in omitted block',y='Change from full-data estimate',title='Leave-one-block-out influence')+theme_bw(base_size=9)+
    theme(strip.text.x=element_text(angle=90))
  ggsave(file.path(bootstrap_run_dir,'block-influence.png'),p,width=17,height=8,dpi=200,bg='white')
  influence_top<-do.call(rbind,lapply(split(influence_long,influence_long$pair),function(z)
    head(z[order(abs(z$change),decreasing=TRUE),],5)))
  # Maps: same Denmark extent/inset/headers as temporal-turnover maps.
  for(w in c(10,20,40,80)) {
    panels<-lapply(c('unshifted','east','north','both'),function(origin_i) {
      z<-boot$layouts[[paste0(w,'km-',origin_i)]]
      off<-c(if(origin_i %in% c('east','both')) w/2 else 0,if(origin_i %in% c('north','both')) w/2 else 0)
      keys<-z[!duplicated(z$block),c('block','bx','by')]
      polys<-lapply(seq_len(nrow(keys)),function(i) {
        x<-keys$bx[i]*w+off[1]; y<-keys$by[i]*w+off[2]
        sf::st_polygon(list(1000*rbind(c(x,y),c(x+w,y),c(x+w,y+w),c(x,y+w),c(x,y))))
      })
      shapes<-sf::st_transform(sf::st_sf(block=factor(keys$block),geometry=sf::st_sfc(polys,crs=25832)),4326)
      z$block<-factor(z$block,levels=keys$block)
      # Deterministic palette shared by mainland and inset within each layout.
      palette<-setNames(grDevices::hcl.colors(nrow(keys),'Dynamic'),keys$block)
      draw<-function(bbox,border=FALSE) ggplot()+
        geom_sf(data=shapes,aes(fill=block),colour='grey65',linewidth=.12,alpha=.18)+
        geom_point(data=z,aes(lon,lat,fill=block),shape=21,colour='white',stroke=.12,size=1.4)+
        scale_fill_manual(values=palette)+coord_sf(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE,datum=NA)+
        theme_void()+theme(legend.position='none',plot.background=element_rect(fill='white',colour=NA),
          panel.border=if(border) element_rect(colour='grey35',fill=NA,linewidth=.55) else element_blank(),
          plot.margin=ggplot2::margin(0,0,0,0))
      composed<-cowplot::ggdraw(draw(mainland_bbox))+
        cowplot::draw_plot(draw(bornholm_bbox,TRUE),x=.79,y=.68,width=.18,height=.22)
      cowplot::plot_grid(map_strip(paste(origin_i,'|',nrow(keys),'blocks')),composed,ncol=1,rel_heights=c(.07,1))
    })
    title<-cowplot::ggdraw()+cowplot::draw_label(paste(w,'km blocks: four grid origins'),fontface='bold',size=15)
    p<-cowplot::plot_grid(title,cowplot::plot_grid(plotlist=panels,ncol=2),ncol=1,rel_heights=c(.06,1))
    ggsave(file.path(bootstrap_run_dir,paste0('block-layouts-',w,'km.png')),p,width=12,height=12,dpi=200,bg='white')
  }
}
