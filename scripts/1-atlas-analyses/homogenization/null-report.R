# Read saved null results only; no compiled sampler or simulations are launched.
source(file.path(out,'null-functions.R'),local=TRUE)
null_root<-file.path(out,'null-results/independent-change')
null_pointer<-file.path(null_root,'latest-run.txt');null_available<-file.exists(null_pointer)
if(null_available) {
  run_name<-readLines(null_pointer,warn=FALSE)
  if(length(run_name)!=1||!identical(basename(run_name),run_name)) stop('Invalid null latest-run pointer')
  null_run_dir<-file.path(null_root,run_name)
  if(!file.exists(file.path(null_run_dir,'run-status.txt'))||
    !identical(readLines(file.path(null_run_dir,'run-status.txt'),warn=FALSE),'complete')) stop('Latest null run is incomplete')
  null<-new.env();load(file.path(null_run_dir,'null-results.RData'),null)
  if(!identical(null$metadata$input_code_hashes,null_hashes(root,out)))
    stop('Saved null results are incompatible with current inputs or calculation code; rerun the explicit null runner.')
  stopifnot(identical(null$sites,sites),max(abs(null$beta-beta))<1e-11)
  null_table<-null$inference;rownames(null_table)<-NULL
  null_table$observed_direction<-ifelse(null_table$observed_delta<0,'Homogenisation',
    ifelse(null_table$observed_delta>0,'Heterogenisation','No change'))
  null_table$null_departure<-ifelse(null_table$observed_minus_null<0,
    'Less spatial differentiation than null mean','More spatial differentiation than null mean')
  # Bin each null distribution independently of the observed reference line,
  # avoiding loss of histogram resolution when the observation lies far outside it.
  hist_data<-do.call(rbind,lapply(names(null$observed),function(pair) {
    h<-hist(null$draws$delta[null$draws$pair==pair],breaks=40,plot=FALSE)
    data.frame(pair=pair,left=head(h$breaks,-1),right=tail(h$breaks,-1),density=h$density)
  }))
  p<-ggplot(hist_data)+geom_rect(aes(xmin=left,xmax=right,ymin=0,ymax=density),
    fill='#B7D7E2',colour='white',linewidth=.15)+
    geom_vline(data=null_table,aes(xintercept=null_mean,colour='Null mean'),linetype='dashed',linewidth=.7)+
    geom_vline(data=null_table,aes(xintercept=observed_delta,colour='Observed'),linewidth=.8)+
    facet_wrap(~pair,ncol=1,scales='free')+
    scale_colour_manual(values=c('Null mean'='#287D8E',Observed='#C86642'),name=NULL)+
    labs(x='Change in mean spatial Jaccard (later minus fixed earlier)',y='Null density',
      title='Spatial homogenisation under independent allocation of changes',
      subtitle=paste(format(null$metadata$settings$simulations,big.mark=','),'draws per comparison; null shape and observed change'),
      caption='Negative changes are homogenisation. Position relative to the null is a separate comparison.\nAxes vary by comparison. The null preserves species gains, losses and endpoint occupancies.')+
    theme_minimal(base_size=11)+theme(legend.position='top',plot.caption=element_text(hjust=0),panel.grid.minor=element_blank())
  ggsave(file.path(null_run_dir,'null-distributions.png'),p,width=10,height=10,dpi=200,bg='white')
  # A second view resolves the shape without stretching the axis to the observation.
  pzoom<-ggplot(hist_data)+geom_rect(aes(xmin=left,xmax=right,ymin=0,ymax=density),
    fill='#B7D7E2',colour='white',linewidth=.15)+
    geom_vline(data=null_table,aes(xintercept=null_mean),colour='#287D8E',linetype='dashed')+
    facet_wrap(~pair,ncol=1,scales='free')+
    labs(x='Simulated change in mean spatial Jaccard',y='Null density',
      title='Null distributions at their own scale',
      subtitle='Observed changes are omitted from this zoom; see the preceding figure and table',
      caption='Dashed lines: null means. Axes vary by comparison; these are null ranges, not bootstrap intervals.')+
    theme_minimal(base_size=11)+theme(plot.caption=element_text(hjust=0),panel.grid.minor=element_blank())
  ggsave(file.path(null_run_dir,'null-distributions-zoom.png'),pzoom,width=10,height=9,dpi=200,bg='white')
}
