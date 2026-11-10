# Add two-sided site tests and run-specific PNG maps to a completed simulation.
analyse_site_turnover_run <- function(run_dir, alpha=0.05) {
  result_path <- file.path(run_dir,'simulation-results.RData')
  stopifnot(file.exists(result_path),is.finite(alpha),alpha>0,alpha<1)
  e <- new.env(parent=emptyenv()); load(result_path,envir=e)
  tab <- e$site_summary
  coord_ix <- match(tab$site,e$coords$site)
  stopifnot(!anyNA(coord_ix))
  tab$lon <- e$coords$lon[coord_ix]
  tab$lat <- e$coords$lat[coord_ix]
  B <- e$n_sim
  tab$excess <- tab$observed_tbi-tab$null_mean
  tab$p_lower <- (tab$below+tab$equal+1)/(B+1)
  tab$p_upper <- (tab$above+tab$equal+1)/(B+1)
  tab$p_two_sided <- pmin(1,2*pmin(tab$p_lower,tab$p_upper))
  tab$direction <- ifelse(tab$excess>0,'More than expected',
    ifelse(tab$excess<0,'Less than expected','At null mean'))
  # The family is all site-by-atlas-pair hypotheses (1698 x 3 here).
  tab$q_holm <- p.adjust(tab$p_two_sided,method='holm')
  tab$q_bh <- NULL
  tab$q_by <- NULL
  e$site_summary <- tab
  for (pair in names(e$result)) {
    ix <- match(paste(pair,e$sites),paste(tab$pair,tab$site))
    stopifnot(!anyNA(ix))
    e$result[[pair]]$site_p_lower <- tab$p_lower[ix]
    e$result[[pair]]$site_p_upper <- tab$p_upper[ix]
    e$result[[pair]]$site_p_two_sided <- tab$p_two_sided[ix]
    e$result[[pair]]$site_q_holm <- tab$q_holm[ix]
    e$result[[pair]]$site_excess <- tab$excess[ix]
    e$result[[pair]]$site_q_bh <- NULL
    e$result[[pair]]$site_q_by <- NULL
  }
  adj_methods <- c('Unadjusted'='p_two_sided','Holm–Bonferroni'='q_holm')
  e$site_inference_summary <- do.call(rbind,lapply(names(adj_methods),function(method) {
    qcol <- unname(adj_methods[method])
    do.call(rbind,lapply(names(e$result),function(pair) {
      z <- tab[tab$pair==pair,,drop=FALSE]
      hit <- z[[qcol]]<alpha
      data.frame(adjustment=method,pair=pair,alpha=alpha,
        significant=sum(hit),more=sum(hit & z$direction=='More than expected'),
        less=sum(hit & z$direction=='Less than expected'))
    }))
  }))
  e$metadata$site_test <- 'two-sided equal-tail Monte Carlo p; double smaller inclusive tail'
  e$metadata$site_p_adjustment_family <- 'all site-by-pair tests'
  e$metadata$site_alpha <- alpha
  save(list=ls(envir=e,all.names=TRUE),envir=e,file=result_path,compress='xz')

  suppressPackageStartupMessages(library(ggplot2))
  pair_labels <- c('1-2'='Atlas 1 → 2','2-3'='Atlas 2 → 3','1-3'='Atlas 1 → 3')
  effect_limit <- max(abs(tab$excess),na.rm=TRUE)
  mainland_bbox <- c(xmin=7.95,xmax=13.0,ymin=54.4,ymax=57.85)
  bornholm_bbox <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
  facet_strip <- function(label) cowplot::ggdraw() + cowplot::draw_label(
    label,x=0.5,y=0.52,hjust=0.5,vjust=0.5,fontface='bold',size=11,
    colour='black') + theme(plot.background=element_rect(fill='grey92',colour='white'),
    plot.margin=ggplot2::margin(0,0,1,0))
  make_map <- function(qcol,adjustment,file_name) {
    hits <- tab[tab[[qcol]]<alpha,,drop=FALSE]
    point_map <- function(pair, bbox, border=FALSE, point_size=1.55) {
      all_pair <- tab[tab$pair==pair,,drop=FALSE]
      hit_pair <- hits[hits$pair==pair,,drop=FALSE]
      ggplot() +
        geom_point(data=all_pair,aes(lon,lat),inherit.aes=FALSE,
          colour='#D4D4D4',size=point_size) +
        geom_point(data=hit_pair,aes(lon,lat,fill=excess),inherit.aes=FALSE,
          shape=21,colour='white',stroke=0.18,size=point_size+0.3) +
        scale_fill_gradient2(low='#2C7FB8',mid='#FFF7E6',high='#D95F0E',
          midpoint=0,limits=c(-effect_limit,effect_limit),
          name='Observed − null\nmean TBI') +
        coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],
          expand=FALSE) +
        theme_minimal(base_size=10) +
        theme(panel.grid=element_blank(),panel.background=element_blank(),
          axis.text=element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),
          legend.position='none',plot.background=element_rect(fill='white',colour=NA),
          panel.border=if(border) element_rect(colour='grey35',fill=NA,linewidth=0.55)
            else element_blank(),plot.margin=ggplot2::margin(0,0,0,0))
    }
    pair_maps <- lapply(names(pair_labels),function(pair) {
      main <- point_map(pair,mainland_bbox)
      inset <- point_map(pair,bornholm_bbox,border=TRUE,point_size=1.7) +
        theme_void()+theme(legend.position='none',
          plot.background=element_rect(fill='white',colour=NA),
          panel.border=element_rect(colour='grey35',fill=NA,linewidth=0.55),
          plot.margin=ggplot2::margin(0,0,0,0))
      map <- cowplot::ggdraw(main)+cowplot::draw_plot(inset,x=0.79,y=0.68,
        width=0.18,height=0.22)+theme(plot.margin=ggplot2::margin(0,0,0,0))
      cowplot::plot_grid(facet_strip(unname(pair_labels[pair])),map,ncol=1,
        rel_heights=c(0.07,1))
    })
    legend_plot <- ggplot(hits[hits$pair=='1-2',],aes(lon,lat,fill=excess))+
      geom_point(shape=21,size=2)+
      scale_fill_gradient2(low='#2C7FB8',mid='#FFF7E6',high='#D95F0E',
        midpoint=0,limits=c(-effect_limit,effect_limit),
        name='Observed − null\nmean TBI')+
      theme_void()+theme(legend.position='right',legend.title=element_text(face='bold'),
        legend.background=element_rect(fill='white',colour=NA))
    legend <- cowplot::get_legend(legend_plot)
    panels <- cowplot::plot_grid(plotlist=pair_maps,nrow=1,align='hv')
    body <- cowplot::plot_grid(panels,legend,nrow=1,rel_widths=c(1,0.13))
    title <- cowplot::ggdraw()+
      cowplot::draw_label('Sites with turnover different from the independent-change null',
        x=0,y=0.78,hjust=0,fontface='bold',size=14)+
      cowplot::draw_label(paste0(adjustment,' two-sided p < ',alpha,
        '; blue = less than expected, orange = more; grey = not flagged'),
        x=0,y=0.23,hjust=0,size=9.5)
    p <- cowplot::plot_grid(title,body,ncol=1,rel_heights=c(0.12,1))
    ggsave(file.path(run_dir,file_name),p,width=12.5,height=5.7,dpi=300,bg='white')
  }
  make_map('p_two_sided','Unadjusted','sitewise-turnover-unadjusted.png')
  make_map('q_holm','Holm family-wise','sitewise-turnover-holm.png')
  result <- e$result
  pair_order <- names(pair_labels)
  hist_data <- do.call(rbind,lapply(pair_order,function(pair) data.frame(
    pair=pair,pair_label=unname(pair_labels[pair]),
    mean_cell_TBI=result[[pair]]$global_mean_draws)))
  observed_global <- data.frame(pair=pair_order,pair_label=unname(pair_labels[pair_order]),
    observed_mean=vapply(pair_order,function(pair) result[[pair]]$observed_global_mean,
      numeric(1)),upper_p=vapply(pair_order,function(pair)
      result[[pair]]$global_upper_p,numeric(1)))
  observed_global$p_label <- ifelse(observed_global$upper_p<1e-6,
    'upper-tail p < 1e-6',paste0('upper-tail p = ',
      format.pval(observed_global$upper_p,digits=2,eps=1e-6)))
  observed_global$facet_label <- paste0(observed_global$pair_label,
    '  |  observed mean = ',format(round(observed_global$observed_mean,4),nsmall=4),
    '  |  ',observed_global$p_label)
  observed_global$facet_order <- factor(seq_along(pair_order),levels=seq_along(pair_order))
  hist_data$facet_order <- factor(match(hist_data$pair,pair_order),
    levels=seq_along(pair_order))
  facet_labs <- stats::setNames(observed_global$facet_label,
    as.character(seq_along(pair_order)))
  p_null <- ggplot(hist_data,aes(mean_cell_TBI,after_stat(density)))+
    geom_histogram(binwidth=0.0001,boundary=0,fill='#78A9A1',
      colour='white',linewidth=0.12)+
    geom_vline(data=observed_global,aes(xintercept=observed_mean),inherit.aes=FALSE,
      colour='#B33E32',linewidth=0.85)+
    facet_wrap(~facet_order,ncol=1,scales='free_x',
      labeller=ggplot2::as_labeller(facet_labs))+
    scale_x_continuous(expand=expansion(mult=c(0.04,0.04)))+
    labs(x='Grid-wide mean cell TBI',y='Null density',
      title='Independent-change null distributions across atlas pairs',
      subtitle='Each row shows simulated grid-wide mean TBI; red line marks the observed mean')+
    theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill='grey92',colour='white'),
      strip.text.y=element_text(face='bold',colour='black'),
      plot.title=element_text(face='bold',size=14),plot.subtitle=element_text(size=10),
      panel.spacing=grid::unit(0.5,'lines'),
      plot.background=element_rect(fill='white',colour=NA))
  ggsave(file.path(run_dir,'null-distributions-all-pairs.png'),p_null,
    width=9,height=8,dpi=300,bg='white')
  observed_tbi <- do.call(rbind,lapply(pair_order,function(pair) data.frame(
    pair=pair,observed_tbi=result[[pair]]$observed)))
  observed_tbi$facet_order <- factor(match(observed_tbi$pair,pair_order),
    levels=seq_along(pair_order))
  observed_means <- data.frame(facet_order=factor(seq_along(pair_order),
      levels=seq_along(pair_order)),
    mean_tbi=vapply(pair_order,function(pair)mean(result[[pair]]$observed),numeric(1)))
  observed_facet_labs <- stats::setNames(paste0(
    unname(pair_labels[pair_order]),'  |  observed mean TBI = ',
    format(round(observed_means$mean_tbi,4),nsmall=4)),
    as.character(seq_along(pair_order)))
  p_observed <- ggplot(observed_tbi,aes(observed_tbi,after_stat(density)))+
    geom_histogram(binwidth=0.025,boundary=0,fill='#78A9A1',
      colour='white',linewidth=0.15)+
    geom_vline(data=observed_means,aes(xintercept=mean_tbi),inherit.aes=FALSE,
      colour='#B33E32',linewidth=0.85)+
    facet_wrap(~facet_order,ncol=1,labeller=ggplot2::as_labeller(observed_facet_labs))+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),
      expand=expansion(mult=c(0.01,0.01)))+
    labs(x='Site-level TBI',y='Site density',
      title='Observed site-level TBI distributions across atlas pairs',
      subtitle='Each site is weighted equally; red line marks the atlas-wide mean')+
    theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill='grey92',colour='white'),
      strip.text=element_text(face='bold',colour='black'),
      plot.title=element_text(face='bold',size=14),plot.subtitle=element_text(size=10),
      panel.spacing=grid::unit(0.5,'lines'),
      plot.background=element_rect(fill='white',colour=NA))
  ggsave(file.path(run_dir,'observed-site-tbi-distributions.png'),p_observed,
    width=9,height=8,dpi=300,bg='white')
  message('Site tests saved for ',nrow(tab),' site-pair combinations; unadjusted flagged ',
    sum(tab$p_two_sided<alpha),', Holm–Bonferroni flagged ',sum(tab$q_holm<alpha),'.')
  invisible(list(site_summary=tab,summary=e$site_inference_summary))
}

analyse_latest_site_turnover <- function(results_root) {
  pointer <- file.path(results_root,'latest-run.txt')
  stopifnot(file.exists(pointer))
  run_dir <- file.path(results_root,readLines(pointer,warn=FALSE)[1])
  stopifnot(identical(readLines(file.path(run_dir,'run-status.txt'),warn=FALSE)[1],'complete'))
  analyse_site_turnover_run(run_dir)
}
