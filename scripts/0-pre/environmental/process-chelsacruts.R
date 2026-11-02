chelsacruts_manifest <- function(atlas,base_url='https://os.zhdk.cloud.switch.ch/chelsav1/chelsa_cruts') {
  rows<-list()
  for(year in sort(unique(unlist(atlas))))for(variable in c('tmin','tmax','prec'))
    rows[[length(rows)+1L]]<-data.frame(variable,date=sprintf('%d-%02d',year,3:6))
  m<-do.call(rbind,rows);m<-m[order(m$date,m$variable),]
  m$year<-as.integer(substr(m$date,1,4));m$month<-as.integer(substr(m$date,6,7))
  if(any(m$year<1901|m$year>2016))stop('CHELSAcruts covers 1901-2016 only.')
  m$filename<-sprintf('CHELSAcruts_%s_%d_%d_V.1.0.tif',m$variable,m$month,m$year)
  m$url<-paste0(base_url,'/',m$variable,'/',m$filename);rownames(m)<-NULL;m
}
cache_chelsacruts_month <- function(row,grid,cache_path,offline=FALSE) {
  target<-file.path(cache_path,row$variable,row$filename)
  if(file.exists(target)){
    r<-terra::rast(target);g<-terra::project(terra::vect(grid),terra::crs(r))
    e<-as.vector(terra::ext(g));re<-as.vector(terra::ext(r))
    if(re[1]>e[1]||re[2]<e[2]||re[3]>e[3]||re[4]<e[4])stop('Cache crop does not cover grid: ',target)
    return(target)
  }
  if(offline)stop('Offline: missing monthly crop ',target)
  message('Reading/cropping ',row$filename)
  r<-terra::rast(paste0('/vsicurl/',row$url));terra::NAflag(r)<- -32768
  g<-terra::project(terra::vect(grid),terra::crs(r));r<-terra::crop(r,g,snap='out')
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  tmp<-paste0(target,'.partial.tif');on.exit(unlink(tmp),add=TRUE)
  terra::writeRaster(r,tmp,overwrite=TRUE,gdal='COMPRESS=DEFLATE',datatype='FLT4S')
  if(!file.rename(tmp,target))stop('Cannot commit crop: ',target)
  target
}
summarise_chelsacruts <- function(grid,atlas,manifest,cache_path) {
  monthly<-function(variable,dates){
    z<-manifest[manifest$variable==variable&manifest$date %in% dates,];z<-z[match(dates,z$date),]
    if(anyNA(z$date))stop('Incomplete monthly manifest')
    paths<-file.path(cache_path,z$variable,z$filename)
    if(any(!file.exists(paths)))stop('Missing monthly crops')
    r<-terra::rast(paths);if(variable!='prec')r<-r/10;r
  }
  first<-terra::rast(file.path(cache_path,manifest$variable[1],manifest$filename[1]))
  g<-sf::st_transform(grid,terra::crs(first));area<-terra::cellSize(first,unit='km')
  annual<-list()
  for(k in seq_along(atlas))for(year in atlas[[k]]){
    message('Summarising climate: atlas ',k,', ',year)
    dates<-sprintf('%d-%02d',year,3:6)
    temp<-(monthly('tmin',dates)+monthly('tmax',dates))/2
    metrics<-list(tmean_breeding=mean(temp,na.rm=FALSE),prec_breeding=sum(monthly('prec',dates),na.rm=FALSE))
    stack<-terra::rast(metrics);names(stack)<-names(metrics)
    values<-exactextractr::exact_extract(stack,g,'weighted_mean',weights=area,progress=FALSE)
    names(values)<-names(metrics)
    values[]<-lapply(values,function(x){x[!is.finite(x)]<-NA_real_;x})
    annual[[length(annual)+1L]]<-cbind(data.frame(kvadratkod=as.character(grid$kvadratkod),atlas=k,year),values)
  }
  annual<-do.call(rbind,annual)
  periods<-lapply(seq_along(atlas),function(k){
    z<-annual[annual$atlas==k,]
    result<-do.call(rbind,lapply(as.character(grid$kvadratkod),function(site){
      part<-z[z$kvadratkod==site,];stopifnot(nrow(part)==length(atlas[[k]]))
      data.frame(kvadratkod=site,atlas=k,tmean_breeding=mean(part$tmean_breeding),prec_breeding=mean(part$prec_breeding))
    }))
    result$survey<-paste0(result$kvadratkod,'_',k)
    result$climate_start_year<-min(atlas[[k]]);result$climate_end_year<-max(atlas[[k]]);result$climate_n_years<-length(atlas[[k]])
    result
  })
  list(annual=annual,periods=do.call(rbind,periods))
}
process_chelsacruts <- function(grid,atlas,cache_path,output_path,plan_only=FALSE,workers=1L,offline=FALSE) {
  m<-chelsacruts_manifest(atlas);dir.create(output_path,recursive=TRUE,showWarnings=FALSE)
  write.csv(m,file.path(output_path,'monthly_manifest.csv'),row.names=FALSE)
  if(plan_only)return(m)
  signature<-file.path(cache_path,'crop_config.rds')
  config<-list(geometry=sf::st_as_binary(sf::st_geometry(grid)),crs=sf::st_crs(grid)$wkt,
               base_url='https://os.zhdk.cloud.switch.ch/chelsav1/chelsa_cruts')
  if(file.exists(signature)&&!identical(readRDS(signature),config))stop('Climate cache belongs to a different grid/source')
  dir.create(cache_path,recursive=TRUE,showWarnings=FALSE)
  task<-function(i)tryCatch(cache_chelsacruts_month(m[i,],grid,cache_path,offline),error=function(e)structure(conditionMessage(e),class='crop_error'))
  results<-if(workers>1L&&.Platform$OS.type!='windows')parallel::mclapply(seq_len(nrow(m)),task,mc.cores=workers)else lapply(seq_len(nrow(m)),task)
  failed<-vapply(results,function(x)inherits(x,'crop_error')||inherits(x,'try-error')||is.null(x),logical(1))
  if(any(failed))stop('Climate crop failure: ',as.character(results[[which(failed)[1]]]))
  if(!file.exists(signature))saveRDS(config,signature)
  result<-summarise_chelsacruts(grid,atlas,m,cache_path)
  write.csv(result$annual,file.path(output_path,'chelsacruts_annual.csv'),row.names=FALSE)
  write.csv(result$periods,file.path(output_path,'chelsacruts_per_atlas.csv'),row.names=FALSE)
  writeLines(c('CHELSAcruts V1.0: https://doi.org/10.16904/envidat.159',
               'March-June mean monthly (tmin+tmax)/2 in degrees C (source divided by 10); precipitation sum in mm/season (unscaled).',
               'Annual metrics averaged across period years. Native pixels, exact overlaps weighted by pixel area.',
               'Missing source -32768/NA excluded spatially; missing months and years not silently omitted.',
               'Atlas 3 climate 2014-2016; HILDA 2014-2017.'),file.path(output_path,'README.txt'))
  result
}
write_x_environment <- function(grid,atlas,climate,hilda_path,output_file) {
  cols<-c('survey','kvadratkod','atlas','tmean_breeding','prec_breeding','hh',paste0('pct_class_',hilda_codes))
  pieces<-lapply(seq_along(atlas),function(k){
    h<-read.csv(file.path(hilda_path,paste0('hilda_atlas',k,'.csv')));c<-climate[climate$atlas==k,]
    sites<-as.character(grid$kvadratkod)
    stopifnot(!anyDuplicated(h$kvadratkod),!anyDuplicated(c$kvadratkod),setequal(sites,h$kvadratkod),setequal(sites,c$kvadratkod))
    h<-h[match(sites,h$kvadratkod),];c<-c[match(sites,c$kvadratkod),]
    stopifnot(all(h$start_year==min(atlas[[k]])),all(h$end_year==max(atlas[[k]])))
    z<-cbind(h,c[c('tmean_breeding','prec_breeding')]);z[cols]
  })
  x<-do.call(rbind,pieces)
  stopifnot(nrow(x)==nrow(grid)*length(atlas),!anyDuplicated(x$survey),
            all(abs(rowSums(x[paste0('pct_class_',hilda_codes)])-100)<1e-7))
  dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
  write.csv(x,output_file,row.names=FALSE);x
}
