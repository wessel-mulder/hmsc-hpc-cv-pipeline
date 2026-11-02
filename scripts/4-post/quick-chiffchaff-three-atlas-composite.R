# Recorded chiffchaff distributions; shared sites after coverage filtering; PNG only.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
d <- read.csv('data-new/derived/design/studyDesign.csv')
y <- read.csv('data-new/derived/occurrences/Y_occurrences.csv')
g <- st_transform(st_read('data-new/derived/atlas-grid/dof-5km-grid/DK5km_ED50grid_approx_kvadrkod_DOF.shp',quiet=TRUE),25832)
stopifnot(!anyDuplicated(d$survey), !anyDuplicated(y$survey), !anyDuplicated(g$kvadratkod), all(d$survey %in% y$survey), all(d$site %in% g$kvadratkod))
d$pa <- y$Phylloscopus_collybita[match(d$survey,y$survey)]
stopifnot(!anyNA(d$pa),all(d$pa %in% c(0,1)))

# Match the pipeline coverage rule: retain Good (1) and Average (2).
# Atlas 2 uses the corrected UNDS_KOD_1 field, as in S0_Data_Definitions.R.
coverage_file <- 'data-new/raw/effort/dækning atlas1+2.xlsx'
c1 <- readxl::read_excel(coverage_file,sheet='kvadrater atlas1')
c2 <- readxl::read_excel(coverage_file,sheet='kvadrater atlas2')
stopifnot(!anyDuplicated(c1$KVADRATKOD),!anyDuplicated(c2$KVADRATKOD))
good1 <- as.character(c1$KVADRATKOD[c1$UNDS_KODE_ %in% c(1,2)])
good2 <- as.character(c2$KVADRATKOD[c2$UNDS_KOD_1 %in% c(1,2)])
shared <- Reduce(intersect,list(d$site[d$atlas==1],d$site[d$atlas==2],d$site[d$atlas==3],good1,good2))
stopifnot(length(shared)>0)
d <- d[d$site %in% shared,]
stopifnot(all(table(d$site)==3),all(table(d$atlas)==length(shared)),all(d$site %in% good1),all(d$site %in% good2))

map <- function(a,inset=FALSE) {
 z <- d[d$atlas==a,]; s <- g[match(z$site,g$kvadratkod),]
 s$pa <- factor(z$pa,levels=c(0,1),labels=c('not-recorded','recorded'))
 ggplot(s)+geom_sf(aes(fill=pa),colour=NA)+
 scale_fill_manual(values=c('not-recorded'='#E6E6E6','recorded'='#285B47'),drop=FALSE,name=NULL)+
 coord_sf(xlim=if(inset)c(855000,905000) else c(400000,750000),ylim=if(inset)c(6100000,6160000) else c(6000000,6450000),expand=FALSE)+
 theme_void(base_size=14)+theme(legend.position='none',plot.background=element_rect(fill='white',colour=NA))
}
years <- c('1970s','1990s','2010s')
panels <- lapply(1:3,function(a) {
 main <- map(a)+labs(title=years[a])+theme(plot.title=element_text(size=20,hjust=.5),plot.margin=margin(12,10,0,10))
 inset <- map(a,TRUE)+theme(panel.border=element_rect(fill=NA,colour='#AAAAAA',linewidth=.35))
 ggdraw(main)+draw_plot(inset,x=.73,y=.69,width=.17,height=.19)
})
legend <- get_legend(map(1)+theme(legend.position='bottom',legend.text=element_text(size=14),legend.key.size=unit(.55,'cm')))
p <- plot_grid(plot_grid(plotlist=panels,nrow=1),legend,ncol=1,rel_heights=c(1,.10))
out <- 'manuscript-outputs/figures/exploratory/2026-10-13-chiffchaff-atlas1-2-3.png'
ggsave(out,p,width=15,height=6.5,dpi=240,bg='white')
print(aggregate(pa~atlas,d,function(x)c(cells=length(x),recorded=sum(x))))
message('Saved: ',out)
