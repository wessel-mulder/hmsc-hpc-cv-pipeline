# Small scientific checks; no repository inputs or outputs are changed.
suppressPackageStartupMessages({library(terra);library(sf)})
source('scripts/0-pre/environmental/process-hilda.R')
source('scripts/0-pre/environmental/process-chelsacruts.R')
d<-tempfile('environment-tests-');dir.create(d)
near<-function(a,b)isTRUE(all.equal(a,b,tolerance=1e-6))
r<-rast(nrows=1,ncols=4,xmin=10,xmax=10.04,ymin=55,ymax=55.01,crs='EPSG:4326')
values(r)<-c(0,77,41,42)
writeRaster(r,file.path(d,'hilda_plus_states_1971_GLOB-v2_wgs84.tif'))
g<-st_as_sf(as.polygons(ext(r),crs=crs(r)));g$kvadratkod<-'test'
h<-process_hilda_v2(g,list(1971),d,file.path(d,'hilda'))
stopifnot(near(h$annual$hh,.75),near(h$annual$pct_class_41,25),near(h$annual$pct_class_42,25),
          near(h$annual$pct_class_0,25),near(h$annual$pct_class_77,25))
# Annual diversity must be averaged before pooling period class percentages.
a<-h$annual[rep(1,2),];a$year<-1971:1972;a[paste0('pct_class_',hilda_codes)]<-0
a$pct_class_41[1]<-100;a$pct_class_42[2]<-100
h<-summarise_hilda_classes(a,list(1971:1972),file.path(d,'hilda'))
stopifnot(near(h$periods[[1]]$hh,0))
# Breeding-only climate: mean C from raw C*10, sum monthly mm, no-data cells retained.
m<-chelsacruts_manifest(list(1971:1972),'https://example.invalid')
r<-rast(nrows=1,ncols=2,xmin=10,xmax=10.02,ymin=55,ymax=55.01,crs='EPSG:4326')
for(i in seq_len(nrow(m))){
  z<-m[i,];v<-if(z$variable=='tmin')100 else if(z$variable=='tmax')200 else 10
  if(z$year==1972)v<-v+if(z$variable=='prec')10 else 20
  values(r)<-c(v,NA);path<-file.path(d,z$variable,z$filename)
  dir.create(dirname(path),showWarnings=FALSE);writeRaster(r,path,overwrite=TRUE)
}
values(r)<-1:2;g<-st_as_sf(as.polygons(r));g$kvadratkod<-c('A','B')
c<-summarise_chelsacruts(g,list(1971:1972),m,d)
stopifnot(near(c$annual$tmean_breeding[c$annual$kvadratkod=='A'],c(15,17)),
          near(c$annual$prec_breeding[c$annual$kvadratkod=='A'],c(40,80)),
          near(c$periods$tmean_breeding[1],16),is.na(c$periods$tmean_breeding[2]),all(m$month %in% 3:6))
unlink(file.path(d,'prec',m$filename[m$variable=='prec'&m$date=='1971-03']))
stopifnot(inherits(try(summarise_chelsacruts(g,list(1971:1972),m,d),silent=TRUE),'try-error'))
unlink(d,recursive=TRUE)
cat('PASS: original classes, water-inclusive diversity, annual aggregation, breeding units, missing cells and missing-month failure\n')
