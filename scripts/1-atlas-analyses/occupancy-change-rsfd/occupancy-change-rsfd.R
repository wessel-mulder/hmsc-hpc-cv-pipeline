# Run from the repository root or this analysis subfolder.
root <- normalizePath(if (dir.exists('data-shared')) '.' else '../../..')
out <- file.path(root, 'scripts/1-atlas-analyses/occupancy-change-rsfd')
design <- read.csv(file.path(root, 'data-new/0-pre/raw/design/studyDesign.csv'), check.names=FALSE)
environment <- read.csv(file.path(root, 'data-new/0-pre/derived/environmental/X_environment.csv'), check.names=FALSE)
stopifnot(all(c('site','atlas','survey') %in% names(design)),
          all(c('kvadratkod','atlas','pct_class_0','pct_class_77') %in% names(environment)),
          !anyDuplicated(design$survey), !anyDuplicated(environment[c('kvadratkod','atlas')]),
          setequal(design$survey, environment$survey),
          all(environment$survey == paste0(environment$kvadratkod, '_', environment$atlas)))
environment$water_percent <- environment$pct_class_0 + environment$pct_class_77
stopifnot(!anyNA(environment$water_percent), all(environment$water_percent >= 0 & environment$water_percent <= 100 + 1e-6))
eligible_sites <- lapply(1:3, function(a) intersect(
  as.character(environment$kvadratkod[environment$atlas == a & environment$water_percent <= 75]),
  as.character(design$site[design$atlas == a])))
sites <- sort(Reduce(intersect, eligible_sites))
effort_path <- file.path(root, 'data-new/0-pre/raw/effort/dækning atlas1+2.xlsx')
effort1 <- readxl::read_excel(effort_path, sheet='kvadrater atlas1')
effort2 <- readxl::read_excel(effort_path, sheet='kvadrater atlas2')
good1 <- as.character(effort1$KVADRATKOD[effort1[['UNDS_KODE_']] %in% c(1, 2)])
good2 <- as.character(effort2$KVADRATKOD[effort2[['UNDS_KOD_1']] %in% c(1, 2)])
sites <- sort(Reduce(intersect, list(sites, good1, good2)))
y <- read.csv(file.path(root,'data-new/0-pre/raw/occurrences/Y_occurrences.csv'), row.names=1, check.names=FALSE)
matrices <- lapply(1:3,function(a) {
 ids <- paste(sites,a,sep='_'); stopifnot(all(ids %in% rownames(y)))
 z <- as.matrix(y[ids,,drop=FALSE]); rownames(z)<-sites; z
})
N <- length(sites); S <- ncol(y); total <- N*S
stopifnot(!anyDuplicated(rownames(y)),!anyDuplicated(names(y)),all(vapply(matrices,function(z) !anyNA(z)&&all(z%in%c(0,1)),logical(1))))
atlas_summary <- do.call(rbind,lapply(1:3,function(a) data.frame(atlas=a,cells=N,species=S,occupied=sum(matrices[[a]]),occupied_percent=100*mean(matrices[[a]]),mean_richness=mean(rowSums(matrices[[a]])))))
species_occupancy <- data.frame(species=names(y),sapply(matrices,colSums),check.names=FALSE)
names(species_occupancy)[-1]<-paste0('atlas',1:3)
cell_richness <- data.frame(site=sites,sapply(matrices,rowSums)); names(cell_richness)[-1]<-paste0('atlas',1:3)
ratio <- function(x,d) ifelse(d==0,NA_real_,100*x/d)
summarise_transitions <- function(p,l,g,a) data.frame(persistence=p,losses=l,gains=g,continued_absence=a,earlier=p+l,later=p+g,changed=l+g,changed_percent=ratio(l+g,p+l+g+a),net=g-l,relative_change_percent=ratio(g-l,p+l),persistence_percent=ratio(p,p+l),loss_percent=ratio(l,p+l),gain_percent=ratio(g,g+a),jaccard_change_percent=ratio(l+g,p+l+g))
pairs <- list(c(1,2),c(2,3),c(1,3))
species_transitions<-cell_transitions<-pairwise<-list()
for(i in seq_along(pairs)) {
 ab<-pairs[[i]]; x<-matrices[[ab[1]]]; z<-matrices[[ab[2]]]; label<-paste(ab,collapse='-')
 transitions<-list(x==1 & z==1,x==1 & z==0,x==0 & z==1,x==0 & z==0)
 counts<-lapply(transitions,sum)
 pairwise[[i]]<-cbind(pair=label,do.call(summarise_transitions,counts))
 species_transitions[[i]]<-cbind(pair=label,species=names(y),do.call(summarise_transitions,lapply(transitions,colSums)))
 cell_transitions[[i]]<-cbind(pair=label,site=sites,do.call(summarise_transitions,lapply(transitions,rowSums)))
 stopifnot(sum(unlist(counts))==total,sum(z)-sum(x)==counts[[3]]-counts[[2]],sum(x!=z)==counts[[2]]+counts[[3]])
 stopifnot(all(colSums(z)-colSums(x)==species_transitions[[i]]$net),all(rowSums(z)-rowSums(x)==cell_transitions[[i]]$net))
}
pairwise<-do.call(rbind,pairwise); species_transitions<-do.call(rbind,species_transitions);cell_transitions<-do.call(rbind,cell_transitions)
stopifnot(all(colSums(species_occupancy[-1])==atlas_summary$occupied),all(colSums(cell_richness[-1])==atlas_summary$occupied))
library(ggplot2)
plots<-list(occupancy=ggplot(atlas_summary,aes(factor(atlas),occupied))+geom_col(fill='#287D8E')+labs(x='Atlas',y='Occupied species–cell entries'), gains_losses=ggplot(species_transitions,aes(losses,gains))+geom_abline(slope=1,intercept=0,colour='grey60')+geom_point(alpha=.5)+facet_wrap(~pair)+coord_equal()+labs(x='Cells with recorded losses',y='Cells with recorded gains'))

for(n in names(plots)) ggsave(file.path(out,paste0('occupancy-change-',n,'.png')),plots[[n]]+theme_bw(base_size=12),width=9,height=5,dpi=180)
message('Validated ',N,' common cells × ',S,' species; ',total,' entries per atlas.')

# One shared figure: occupancy levels and the gains/losses behind their differences.
number_label <- function(x) format(x, big.mark=',', scientific=FALSE, trim=TRUE)
occupancy_panel <- ggplot(atlas_summary,aes(atlas,occupied)) +
 geom_line(colour='#287D8E',linewidth=1) + geom_point(colour='#287D8E',size=3) +
 geom_text(aes(label=number_label(occupied)),vjust=-1,size=4) +
 scale_x_continuous(breaks=1:3,labels=paste('Atlas',1:3)) +
 scale_y_continuous(labels=number_label,limits=c(0,125000),expand=expansion(mult=c(0,.02))) +
 labs(x=NULL,y='Recorded presences',title='Total recorded occupancy') + theme_bw(base_size=12)
flows <- rbind(data.frame(pair=pairwise$pair,type='Gains',entries=pairwise$gains),
               data.frame(pair=pairwise$pair,type='Losses',entries=-pairwise$losses))
flows$pair <- factor(flows$pair,levels=pairwise$pair)
flow_panel <- ggplot(flows,aes(pair,entries,fill=type)) +
 geom_col(width=.55) + geom_hline(yintercept=0,colour='grey35') +
 geom_text(aes(label=number_label(abs(entries)),vjust=ifelse(entries>0,-.5,1.4)),size=3.8) +
 geom_point(data=pairwise,aes(x=pair,y=net),inherit.aes=FALSE,shape=23,fill='white',size=3) +
 scale_fill_manual(values=c(Gains='#238B72',Losses='#C86642')) +
 scale_x_discrete(labels=c('1-2'='Atlas 1 → 2','2-3'='Atlas 2 → 3','1-3'='Atlas 1 → 3')) +
 scale_y_continuous(labels=number_label,expand=expansion(mult=c(.2,.18))) +
 labs(x=NULL,y='Gains (+) and losses (−)',fill=NULL,title='Entries gained and lost',
 subtitle='White diamonds show net change (gains minus losses)') +
 theme_bw(base_size=12) + theme(legend.position='top')
combined <- cowplot::plot_grid(occupancy_panel,flow_panel,ncol=1,labels=c('A','B'),align='v',axis='lr',rel_heights=c(1,1.25))
ggsave(file.path(out,'occupancy-change-combined.png'),combined,width=9,height=9,dpi=200,bg='white')

# Range-size frequency distribution on the identical common-cell sample.
# Zero-record species in an atlas remain in the fixed species denominator.
range_size <- do.call(rbind,lapply(1:3,function(a) data.frame(
 species=colnames(matrices[[a]]),atlas=a,
 occupied_cells=colSums(matrices[[a]]),
 occupied_percent=100*colSums(matrices[[a]])/N)))
rsfd <- do.call(rbind,lapply(1:3,function(a) {
 values <- range_size$occupied_percent[range_size$atlas==a]
 thresholds <- sort(unique(c(0,values,100)))
 data.frame(atlas=a,occupied_percent=thresholds,
 cumulative_species_fraction=vapply(thresholds,function(t) mean(values<=t),numeric(1)))
}))
for(a in 1:3) {
 z<-rsfd[rsfd$atlas==a,]
 stopifnot(sum(range_size$atlas==a)==S,all(diff(z$cumulative_species_fraction)>=0),
 tail(z$cumulative_species_fraction,1)==1,
 all(z$cumulative_species_fraction>=0 & z$cumulative_species_fraction<=1),
 all(range_size$occupied_percent[range_size$atlas==a]>=0 & range_size$occupied_percent[range_size$atlas==a]<=100))
 # Check the zero boundary, including species absent from this atlas.
 stopifnot(z$cumulative_species_fraction[1]==mean(colSums(matrices[[a]])==0))
}
rsfd_plot <- ggplot(rsfd,aes(occupied_percent,cumulative_species_fraction,
 colour=factor(atlas),linetype=factor(atlas))) +
 geom_step(direction='hv',linewidth=.9) +
 scale_colour_manual(values=c('#287D8E','#C86642','#6B4C9A'),labels=paste('Atlas',1:3)) +
 scale_linetype_discrete(labels=paste('Atlas',1:3)) +
 scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20),expand=expansion(mult=0)) +
 scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),expand=expansion(mult=0)) +
 labs(x='Percentage of common grid cells occupied',y='Cumulative species fraction',
 colour=NULL,linetype=NULL,title='Cumulative distribution of species occupancy') +
 theme_bw(base_size=12) + theme(legend.position='top')
ggsave(file.path(out,'occupancy-change-cumulative-occupancy.png'),rsfd_plot,
 width=8,height=6,dpi=200,bg='white')

# Species-matched occupancy comparisons: each point is the same species at two dates.
matched_occupancy <- do.call(rbind,lapply(list(c(1,2),c(2,3)),function(ab) {
 earlier <- range_size[range_size$atlas==ab[1],]
 later <- range_size[range_size$atlas==ab[2],]
 j <- match(earlier$species,later$species)
 stopifnot(!anyNA(j),!anyDuplicated(earlier$species),nrow(earlier)==S)
 data.frame(species=earlier$species,pair=paste(ab,collapse='-'),
 earlier_percent=earlier$occupied_percent,later_percent=later$occupied_percent[j])
}))
matched_plots <- lapply(c('1-2','2-3'),function(pair) {
 ab <- strsplit(pair,'-')[[1]]
 z <- matched_occupancy[matched_occupancy$pair==pair,]
 stopifnot(nrow(z)==S,all(z$earlier_percent>=0 & z$earlier_percent<=100),
 all(z$later_percent>=0 & z$later_percent<=100))
 ggplot(z,aes(earlier_percent,later_percent)) +
 geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40') +
 geom_point(colour='#287D8E',alpha=.55,size=1.8) +
 scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20)) +
 scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20)) +
 coord_equal() + theme_bw(base_size=12) +
 labs(x=paste0('Atlas ',ab[1],' occupancy (% of common cells)'),
 y=paste0('Atlas ',ab[2],' occupancy (% of common cells)'),
 title=paste0('Atlas ',ab[1],' → ',ab[2]),
 subtitle='Each point is one species; dashed line indicates no change')
})
for(i in seq_along(matched_plots)) ggsave(
 file.path(out,paste0('occupancy-change-species-matched-',c('1-2','2-3')[i],'.png')),
 matched_plots[[i]],width=7,height=7,dpi=200,bg='white')
