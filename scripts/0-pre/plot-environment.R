# Current reduced environmental matrix: PNG maps and breeding climate point cloud.
plot_environment <- function(input_file='data-new/0-pre/derived/environmental/X_environment.csv',
                             output_dir='data-new/0-pre/derived/environmental/figures/current') {
  suppressPackageStartupMessages({library(ggplot2);library(sf);library(cowplot)})
  x<-read.csv(input_file)
  g<-st_read('data-new/0-pre/raw/DOF_Shapefiles_/DK5km_ED50grid_approx_kvadrkod_DOF.shp',quiet=TRUE)
  stopifnot(!anyDuplicated(x$survey),nrow(x)==3*nrow(g))
  for(a in 1:3)stopifnot(setequal(x$kvadratkod[x$atlas==a],g$kvadratkod))
  codes<-c(0,11,22,23,24,33,40:45,55,66,77)
  names_class<-c('Ocean','Urban','Annual crops','Tree crops','Agroforestry','Pasture/rangeland',
                'Forest unknown/other','Evergreen needleleaf forest','Evergreen broadleaf forest',
                'Deciduous needleleaf forest','Deciduous broadleaf forest','Mixed forest',
                'Unmanaged grass/shrub','Sparse/no vegetation','Freshwater')
  labels<-c(tmean_breeding='Breeding temperature (°C)',prec_breeding='Breeding precipitation (mm/season)',
            hh='Habitat heterogeneity (Simpson diversity)',setNames(paste0(names_class,' (%)'),paste0('pct_class_',codes)))
  dir.create(file.path(output_dir,'maps'),recursive=TRUE,showWarnings=FALSE)
  geo<-g[match(x$kvadratkod,g$kvadratkod),];scaling<-list()
  for(v in names(labels)){
    years<-if(grepl('^(tmean|prec)_',v))c('1971–1974','1993–1996','2014–2016')else c('1971–1974','1993–1996','2014–2017')
    geo$period<-factor(paste0('Atlas ',x$atlas,' · ',years[x$atlas]),levels=paste0('Atlas ',1:3,' · ',years))
    available<-Reduce(intersect,lapply(1:3,function(a)x$kvadratkod[x$atlas==a&is.finite(x[[v]])]))
    relative<-rep(NA_real_,nrow(x))
    for(a in 1:3){
      idx<-which(x$atlas==a&x$kvadratkod %in% available);mu<-mean(x[[v]][idx]);s<-sd(x[[v]][idx])
      if(is.finite(s)&&s>1e-12)relative[idx]<-(x[[v]][idx]-mu)/s
      scaling[[length(scaling)+1L]]<-data.frame(variable=v,atlas=a,n=length(idx),mean=mu,sd=s)
    }
    theme_map<-theme_void(base_size=11)+theme(legend.position='bottom',plot.title=element_text(face='bold'),
                                             strip.text=element_text(face='bold'),plot.caption=element_text(hjust=0,size=8))
    geo$value<-x[[v]]
    lim<-if(grepl('^pct_',v))c(0,100)else range(geo$value,na.rm=TRUE)
    if(diff(lim)==0)lim<-lim+c(-.5,.5)
    p<-ggplot(geo)+geom_sf(aes(fill=value),colour=NA)+facet_wrap(~period,nrow=1)+coord_sf(datum=NA)+
      scale_fill_viridis_c(limits=lim,option=if(grepl('tmean',v))'C'else'viridis',name=labels[[v]],na.value='#c8c8c8')+
      labs(title=paste(labels[[v]],'— absolute values'),subtitle='One common scale across atlases')+theme_map+
      guides(fill=guide_colourbar(barwidth=grid::unit(7,'cm'),title.position='top'))
    geo$value<-relative;cap<-if(any(is.finite(relative)))max(abs(relative),na.rm=TRUE)else 1
    if(cap==0)cap<-1
    q<-ggplot(geo)+geom_sf(aes(fill=value),colour=NA)+facet_wrap(~period,nrow=1)+coord_sf(datum=NA)+
      scale_fill_gradient2(low='#2166ac',mid='#f7f7f7',high='#b2182b',limits=c(-cap,cap),name='SD from atlas mean',na.value='#c8c8c8')+
      labs(title='Relative values — standardised within each atlas',subtitle=paste(length(available),'shared non-missing cells'),
           caption='All input cells shown. Grey: missing or zero variance. Climate: March–June; Atlas 3 climate ends in 2016.\nPercentages and heterogeneity include ocean and freshwater; original land-cover classes remain separate.')+
      theme_map+guides(fill=guide_colourbar(barwidth=grid::unit(7,'cm'),title.position='top'))
    ggsave(file.path(output_dir,'maps',paste0(v,'-absolute-relative.png')),plot_grid(p,q,ncol=1),width=14,height=11,dpi=180,bg='white')
  }
  ok<-is.finite(x$tmean_breeding)&is.finite(x$prec_breeding)
  shared<-Reduce(intersect,lapply(1:3,function(a)x$kvadratkod[x$atlas==a&ok]))
  z<-x[x$kvadratkod %in% shared,];z$atlas<-factor(z$atlas,levels=1:3,
    labels=c('Atlas 1 · 1971–1974','Atlas 2 · 1993–1996','Atlas 3 · 2014–2016'))
  set.seed(417);z<-z[sample(nrow(z)),]
  p<-ggplot(z,aes(tmean_breeding,prec_breeding,colour=atlas))+geom_point(size=1.2,alpha=.3)+
    scale_colour_manual(values=c('#0072b2','#d55e00','#009e73'))+
    guides(colour=guide_legend(override.aes=list(alpha=1,size=3)))+
    labs(title='Breeding climate across the three atlases',subtitle=paste(length(shared),'shared cells per atlas'),
         x='March–June mean temperature (°C)',y='March–June precipitation (mm/season)',colour=NULL,
         caption='Each point is one cell in one atlas. Shared complete cells; no jitter. Atlas 3 climate: 2014–2016.')+
    theme_minimal(base_size=13)+theme(legend.position='bottom',panel.grid.minor=element_blank(),plot.title=element_text(face='bold'))
  ggsave(file.path(output_dir,'temperature-precipitation-breeding.png'),p,width=11,height=7,dpi=200,bg='white')
  write.csv(z[c('survey','kvadratkod','atlas','tmean_breeding','prec_breeding')],file.path(output_dir,'point-cloud-data.csv'),row.names=FALSE)
  write.csv(do.call(rbind,scaling),file.path(output_dir,'relative-map-scaling.csv'),row.names=FALSE)
  writeLines(c('Current final environmental schema: 18 predictor maps and one breeding climate point cloud. PNG only.',
               'Absolute maps use a common scale across atlases; relative maps use within-atlas z scores on shared complete sites.',
               'Relative maps remove atlas-wide mean/SD differences; use absolute maps to assess temporal shifts.',
               'Original 15 HILDA classes are separate; hh includes all classes. All 2165 cells mapped.',
               paste('Point cloud:',length(shared),'shared complete cells per atlas.'),
               'Earlier expanded-variable figures elsewhere in figures/ are legacy explorations and superseded.'),file.path(output_dir,'README.txt'))
  invisible(output_dir)
}
if(sys.nframe()==0L)plot_environment()
