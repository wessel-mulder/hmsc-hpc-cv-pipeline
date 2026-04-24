rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra)

#### LOAD DATA ####
ref <- rast('~/box/PhD/logistics/data/databases/Rahbek-range-maps/Images/ref.tif')
load("~/box/PhD/logistics/data/databases/Rahbek-range-maps/dis.RData")


#example1: show the range of a species 
sp="Accipiter nisus"
range <- init(ref, NA)
range[dis[[sp]]]=1
plot(trim(range))

#example 2 make a richness map
range <- init(ref, 0)
for (x in dis) range[x] <- range[x] + 1
range[range == 0] <- NA
plot(trim(range))

lapply(dis,length)
ref
max(dis[[sp]])
plot(ref)
max(ref)
min(ref)
