hilda_codes <- c(0,11,22,23,24,33,40:45,55,66,77)
hilda_labels <- c('Ocean','Urban','Annual crops','Tree crops','Agroforestry','Pasture/rangeland',
                  'Forest unknown/other','Forest evergreen needleleaf','Forest evergreen broadleaf',
                  'Forest deciduous needleleaf','Forest deciduous broadleaf','Forest mixed',
                  'Unmanaged grass/shrub','Sparse/no vegetation','Inland water')

summarise_hilda_classes <- function(annual, atlas, output_path) {
  cols <- paste0('pct_class_',hilda_codes)
  stopifnot(all(cols %in% names(annual)),!anyDuplicated(annual[c('kvadratkod','year')]))
  annual <- annual[c('kvadratkod','year','valid_km2','land_km2','valid_fraction',cols)]
  p <- as.matrix(annual[cols]); p <- p/rowSums(p)
  annual$hh <- 1-rowSums(p^2)
  annual$hh[!is.finite(annual$hh)] <- NA_real_
  dir.create(output_path,recursive=TRUE,showWarnings=FALSE)
  write.csv(annual,file.path(output_path,'hilda_v2_annual.csv'),row.names=FALSE)
  periods <- lapply(seq_along(atlas),function(k) {
    z <- annual[annual$year %in% atlas[[k]],]
    stopifnot(all(table(z$kvadratkod)==length(atlas[[k]])))
    result <- aggregate(z[setdiff(names(z),c('kvadratkod','year'))],
                        list(kvadratkod=z$kvadratkod),mean,na.action=na.pass)
    result$survey <- paste0(result$kvadratkod,'_',k)
    result$atlas <- k; result$start_year <- min(atlas[[k]]); result$end_year <- max(atlas[[k]])
    write.csv(result,file.path(output_path,paste0('hilda_atlas',k,'.csv')),row.names=FALSE)
    result
  })
  write.csv(data.frame(code=hilda_codes,label=hilda_labels),file.path(output_path,'class_dictionary.csv'),row.names=FALSE)
  writeLines(c('HILDA+ v2: https://doi.org/10.1594/PANGAEA.974335',
               'Exact overlap fraction times native pixel area; missing codes 99/NA excluded.',
               'All 15 original classes separate; percentages use all valid area including ocean and water.',
               'hh: mean ANNUAL Simpson diversity over all 15 classes; missing classes are zero, missing source area is not a class.',
               'Period metrics are arithmetic means; missing years are not silently omitted.'),file.path(output_path,'README.txt'))
  list(annual=annual,periods=periods)
}

process_hilda_v2 <- function(grid,atlas,input_path,output_path) {
  stopifnot(inherits(grid,'sf'),!anyDuplicated(grid$kvadratkod),all(sf::st_is_valid(grid)))
  years <- sort(unique(unlist(atlas)))
  files <- file.path(input_path,sprintf('hilda_plus_states_%d_GLOB-v2_wgs84.tif',years))
  if(any(!file.exists(files))) stop('Missing HILDA annual states: ',paste(files[!file.exists(files)],collapse=', '))
  annual <- lapply(seq_along(years),function(k) {
    message('HILDA+ v2: ',years[k])
    r <- terra::rast(files[k]); polygons <- terra::project(terra::vect(grid),terra::crs(r))
    r <- terra::crop(r,polygons,snap='out');names(r)<-'class_code'
    area <- terra::cellSize(r,unit='km');names(area)<-'pixel_km2'
    ex <- terra::extract(c(r,area),polygons,exact=TRUE)
    unknown <- setdiff(unique(ex$class_code[!is.na(ex$class_code)]),c(hilda_codes,99))
    if(length(unknown)) stop('Unexpected HILDA codes: ',paste(unknown,collapse=', '))
    ex$area <- ex$fraction*ex$pixel_km2
    do.call(rbind,lapply(seq_len(nrow(grid)),function(i) {
      z<-ex[ex$ID==i,];valid<-!is.na(z$class_code)&z$class_code!=99
      a<-vapply(hilda_codes,function(code)sum(z$area[valid&z$class_code==code],na.rm=TRUE),numeric(1))
      total<-sum(a);p<-if(total>0)a/total else rep(NA_real_,length(a))
      out<-data.frame(kvadratkod=as.character(grid$kvadratkod[i]),year=years[k],valid_km2=total,
                      land_km2=sum(a[!hilda_codes %in% c(0,77)]),valid_fraction=total/sum(z$area,na.rm=TRUE))
      for(j in seq_along(hilda_codes))out[[paste0('pct_class_',hilda_codes[j])]]<-100*p[j]
      out
    }))
  })
  summarise_hilda_classes(do.call(rbind,annual),atlas,output_path)
}
