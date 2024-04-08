library(terra)
# library(FedData)
# library(plotKML)

setwd('/projectnb/modislc/users/sjstone/planet/')

# data("worldgrids_pal")
# igbp.colors <- as.vector(worldgrids_pal$IGBP)

ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))

modis.path <- '/projectnb/modislc/users/sjstone/phenology_model_project/data/phenology/V061/'
modis.lc.path <- '/projectnb/modislc/users/sjstone/phenology_model_project/data/LC/V061/'
# HLS.path <- '/projectnb/modislc/users/mkmoon/MuSLI/V1_0/From_AWS/product_v011/'
HLS.path <- '/projectnb/modislc/users/sjstone/planet/data/HLS/'
main.path <- '/projectnb/modislc/users/sjstone/planet/'

# cut.off <- 0.8

planet.path <- list.files('/projectnb/planet/PLSP/Archive/Product_010/',pattern='03_2017_50PCGI*',full.names=T,recursive=T)

nasum <- function(x){sum(is.na(x))}
nsamp <- function(x){length(x)}
QAperc <- function(x){(1 - length(which(is.na(x) == T | x > 2)) / length(x))}

# pull in the MODIS grid
sinu.grid <- vect('/projectnb/modislc/users/sjstone/phenology_model_project/data/modis_grid/modis_sinusoidal_grid_world.shp',
                  crs='+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')

# this funciton will compute the mean HLS value of the pixels that fall within a modis pixel, must specify the band number for HLS and modis 
# HLS bands: 2 - 15% increase, 3 - 50% increase, 4 - 90% increase, 6 - 10% decrease, 7 - 50% decrease, 8 - 85% decrease, 10 - EVI amp, 11 - EVI area
# modis bands: 2 - 15% increase, 3 - 50% increase, 4 - 90% increase, 6 - 10% decrease, 7 - 50% decrease, 8 - 85% decrease, 9 - EVI min, 10 - EVI amp, 11 - EVI area
# when running EVI metrics comment out lines 92 and 143

HLS.MODIS <- function(planet.file, year, HLS.band, mod.band){
  # set.seed(107644)
  # finding the difference between 1970 and the selected year for the MODIS phenology data
  date1 <- as.Date("1970-01-01")
  date2 <- as.Date(paste(year,"-01-01",sep=''))
  difference <- as.numeric(difftime(date2,date1,units='days'))
  
  # creating the planet raster
  plan.raster <- rast(planet.file)
  plan.crs <- crs(plan.raster, proj=TRUE)
  plan.site <- strsplit(planet.file,split='/')[[1]][8]
  
  # pulling which UTM zone the raster is in to then get which grid it falls in
  plan_utm <- as.numeric(strsplit(strsplit(plan.crs,split=' ')[[1]][2],'=')[[1]][2])
  plan_utm <- formatC(plan_utm, width = 2, flag = "0")
  
  if(plan.site == 'US-NC3__NC_Clearcut#3'){
    HLS.full.file <- list.files(paste(HLS.path, 'US-NC3__NC_Clearcut_3', sep = ''), glob2rx(paste('*', year, '.nc', sep = '')), full.names = TRUE)
  }
  else{
    HLS.full.file <- list.files(paste(HLS.path, plan.site, sep = ''), glob2rx(paste('*', year, '.nc', sep = '')), full.names = TRUE)
  }
  
  if(isTRUE(length(HLS.full.file)=='0')==TRUE){
    return(list(location=planet.file, data=NA))
  }
  
  if(length(HLS.full.file) > 1){
    HLS.full.file <- HLS.full.file[grep(paste('_', plan_utm, sep = ''), HLS.full.file)]
  }
  
  if(HLS.band == 2 | HLS.band == 3 | HLS.band == 4){
    QAband <- 23
  }
  else{
    QAband <- 24
  }
  
  HLS.GUP <- rast(HLS.full.file[1])[[c(HLS.band, QAband)]]
  
  # have select sites where the HLS data needs to be reprojected to intersect with the planet site
  if(plan.site == 'US-ALQ__Allequash_Creek_Site' | plan.site == 'US-xST__NEON_Steigerwaldt_Land_Services' | 
     plan.site == 'US-xTL__NEON_Toolik' | plan.site == 'US-xTR__NEON_Treehaven' | plan.site == 'US-xUN__NEON_University_of_Notre_Dame_Environmental_Research_Center'){
    HLS.GUP <- terra::project(HLS.GUP, crs(plan.raster), res = 30, method = 'near')
  }
  
  HLS.crop <- terra::crop(HLS.GUP, ext(plan.raster)) # clipping to the extent of the planet site
  HLS.proj <- terra::project(HLS.crop,'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs',method='near',res=30) # reprojecting to the projection of the MODIS data
  
  # locating the modis tiles that intersect with the extent of the clipped HLS files
  modis.tiles <- terra::intersect(sinu.grid,ext(HLS.proj))
  modis.hv <- cbind(modis.tiles$h,modis.tiles$v) # putting the horizontal and vertical tile into a table
  modis.hv <- formatC(modis.hv,width=2,flag="0") # formating the numbers below 10 to have a leading 0
  
  modis.file <- lapply(c(1:dim(modis.hv)[1]), function(x) list.files(modis.path, pattern = glob2rx(paste('*', year, '001.h', modis.hv[x, 1], 'v', modis.hv[x,2],'*', sep = '')), full.names = T))
  modis <- lapply(modis.file, function(x) sds(x))
  modis.rast <- mosaic(sprc(lapply(modis, function(x) x[[mod.band]][[1]])))
  modis.QArast <- mosaic(sprc(lapply(modis, function(x) x[[12]][[1]])))
  
  modis.lc.file <- lapply(c(1:dim(modis.hv)[1]), function(x) list.files(modis.lc.path, pattern = glob2rx(paste('*', year, '001.h', modis.hv[x, 1], 'v', modis.hv[x,2],'*', sep = '')), full.names = T))
  modis.lc.rast <- mosaic(sprc(lapply(modis.lc.file, function(x) sds(x)[1])))
  
  modis.crop <- terra::crop(modis.rast,ext(HLS.proj)) # cropping the MODIS image to the extent of the planet raster
  modis.df <- terra::as.data.frame(modis.crop,xy=T,na.rm=F) # turning the MODIS raster into a data frame
  
  if(mod.band == 10 | mod.band == 11){
    modis.df <- modis.df
  }
  else{
    modis.df[,3] <- modis.df[,3]-difference
  }
  # GUP.na.loc <- which(is.na(modis.df[,3])==T,arr.ind=T)
  
  modis.qa.crop <- terra::crop(modis.QArast, ext(HLS.proj))
  modis.qa.df <- terra::as.data.frame(modis.qa.crop, xy=T, na.rm = F)
  
  modis.lc.crop <- terra::crop(modis.lc.rast,ext(HLS.proj))
  modis.lc.df <- terra::as.data.frame(modis.lc.crop,xy=T)
  LC <- modis.lc.df[,3] # creating a vector of just the LC values
  QA <- modis.qa.df[,3] # creating a vector of just the QA band value
  sample.polygons <- as.polygons(modis.crop,dissolve=F, na.rm = F) # creating polygons from the modis raster
  
  modis.combine <- cbind(modis.df,LC, QA)
  
  # if(length(GUP.na.loc) > 0){
  #   modis.combine <- modis.combine[-GUP.na.loc,]
  # }
  # else{
  #   modis.combine <- modis.combine
  # }
  
  if(dim(modis.combine)[1]==0){
    return(list(location=plan.site,data=NA))
  }
  
  extract.median <- terra::extract(HLS.proj[[1]], sample.polygons, fun = median, na.rm = T)[,2]
  na.count <- terra::extract(HLS.proj[[1]], sample.polygons, fun = nasum)[,2]
  n.HLS <- terra::extract(HLS.proj[[1]], sample.polygons, fun = nsamp)[,2]
  hlsQAperc <- terra::extract(HLS.proj[[2]], sample.polygons, fun = QAperc)[,2]
  
  plan.df <- as.data.frame(cbind(modis.combine,extract.median,na.count,n.HLS, hlsQAperc))
  colnames(plan.df) <- c('x','y','MODIS','LC','modisQA', 'median','na.count','HLS.count', 'hlsQAperc')
  return(list(location=plan.site,data=plan.df))
  
}

# function that is used to identify where there is insufficient HLS data within a MODIS pixel 
# loc <- function(planetdf,threshold){
#   if(is.na(planetdf)[1]==T){
#     insufficient.loc <- NA
#   }
#   else{
#     insufficient.loc <- which(planetdf[,6]/planetdf[,7] > threshold)
#   }
#   return(insufficient.loc)
# }
# 
# # function that calls the above function and then removes and updates the list with only the valid comparison pixels
# RemoveLoc <- function(median.list){
#   # finding the locations that do not meet the minimum HLS data density for the pixel, can adjust 'cut.off' variable above
#   insufficient.loc <- loc(median.list[[2]],cut.off)
#   
#   if(length(insufficient.loc)=='0'){
#     updated.df <- median.list[[2]]
#   }
#   else if(is.na(insufficient.loc)[1] ==T){
#     updated.df <- NA
#   }
#   else{
#     updated.df <- median.list[[2]][-insufficient.loc,]
#   }
#   final.list <- list(location = median.list[[1]], updated_df = updated.df)
#   return(final.list)
# }

# old plotting function, all of the plotting is done in seperate scripts now 
# PlotFun <- function(updated.median,year){
#   df <- updated.median[[2]]
#   site.name <- updated.median[[1]]
#   if(isTRUE(dim(df)[1] == '0')==TRUE){
#     plot(0,1,main=paste('null plot',site.name,year,sep=': '))
#     hist(0,main = paste('null plot',site.name,year,sep=': '))
#   }
#   else if(is.na(df)[1]==TRUE){
#     plot(0,1,main=paste('null plot',site.name,year,sep=': '))
#     hist(0,main = paste('null plot',site.name,year,sep=': '))
#   }
#   else{
#     lm.compare <- summary(lm(df[,'median']~df[,'MODIS']))
#     rsquare <- round((lm.compare$r.squared),2)
#     plot(df[,'MODIS'],df[,'median'],
#          xlab='modis',
#          ylab='HLS median',
#          main=paste(site.name,year,sep=': '),
#          cex.main=0.8,
#          col=igbp.colors[df[,'LC']+1],
#          pch=16)
#     #points(modis[,3],planet[,3], col='red')
#     mtext(paste('n=',dim(df)[1],', R2=',rsquare,sep=" "),cex=0.75)
#     abline(0,1)
#     hist(df[,'LC'],
#          breaks=0:15,
#          main='MCD12Q1',
#          xlab='LC Type',
#          col=igbp.colors[2:16])
#   }
# }

# test <- HLS.MODIS(planet.path[[11]], '2017', 2, 2)

bands <- c(2, 3, 4, 6, 7, 8, 10, 11)
metrics <- c('15GUP', '50GUP', '90GUP', '10GD', '50GD', '85GD', 'EVIamp', 'EVIarea')

selected.band <- bands[ID]
selected.metric <- metrics[ID]

for(yr in 2017:2019){
  HLSmod.comp <- lapply(planet.path, function(x) HLS.MODIS(x, yr, selected.band, selected.band))
  if(selected.metric == '50GUP'){
    saveRDS(HLSmod.comp, paste('/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_', yr, '_wQA.RDS', sep = ''))
  }
  else{
    saveRDS(HLSmod.comp, paste('/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_', selected.metric, '_', yr, '_wQA.RDS', sep = ''))
  }
  print(paste(selected.metric, yr, 'saved'))
}

# 15% GUP comparison
# HLSmod.comp.2017 <- lapply(planet.path, function(x) HLS.MODIS(x,'2017',2,2))
# saveRDS(HLSmod.comp.2017, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_15GUP_2017_wQA.RDS')
# print('saved')

# HLSmod.comp.2018 <- lapply(planet.path, function(x) HLS.MODIS(x,'2018',2,2))
# saveRDS(HLSmod.comp.2018, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_15GUP_2018_wQA.RDS')
# print('saved')
# 
# HLSmod.comp.2019 <- lapply(planet.path, function(x) HLS.MODIS(x,'2019',2,2))
# saveRDS(HLSmod.comp.2019, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_15GUP_2019_wQA.RDS')
# print('saved')

# 90% GUP comparison
# HLSmod.comp.2017 <- lapply(planet.path, function(x) HLS.MODIS(x,'2017',4,4))
# saveRDS(HLSmod.comp.2017, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_90GUP_2017_wQA.RDS')
# print('saved')
# HLSmod.comp.2018 <- lapply(planet.path, function(x) HLS.MODIS(x,'2018',4,4))
# saveRDS(HLSmod.comp.2018, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_90GUP_2018_wQA.RDS')
# print('saved')
# HLSmod.comp.2019 <- lapply(planet.path, function(x) HLS.MODIS(x,'2019',4,4))
# saveRDS(HLSmod.comp.2019, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_90GUP_2019_wQA.RDS')
# print('saved')

# 10% GD comparison
# HLSmod.comp.2017 <- lapply(planet.path, function(x) HLS.MODIS(x,'2017',6,6))
# saveRDS(HLSmod.comp.2017, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_10GD_2017_wQA.RDS')
# HLSmod.comp.2018 <- lapply(planet.path, function(x) HLS.MODIS(x,'2018',6,6))
# saveRDS(HLSmod.comp.2018, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_10GD_2018_wQA.RDS')
# HLSmod.comp.2019 <- lapply(planet.path, function(x) HLS.MODIS(x,'2019',6,6))
# saveRDS(HLSmod.comp.2019, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_10GD_2019_wQA.RDS')

# 50% GD comparison
# HLSmod.comp.2017 <- lapply(planet.path, function(x) HLS.MODIS(x,'2017',7,7))
# saveRDS(HLSmod.comp.2017, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_50GD_2017_wQA.RDS')
# HLSmod.comp.2018 <- lapply(planet.path, function(x) HLS.MODIS(x,'2018',7,7))
# saveRDS(HLSmod.comp.2018, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_50GD_2018_wQA.RDS')
# HLSmod.comp.2019 <- lapply(planet.path, function(x) HLS.MODIS(x,'2019',7,7))
# saveRDS(HLSmod.comp.2019, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_50GD_2019_wQA.RDS')

# 85% GD comparison
# HLSmod.comp.2017 <- lapply(planet.path, function(x) HLS.MODIS(x,'2017',8,8))
# saveRDS(HLSmod.comp.2017, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_85GD_2017_wQA.RDS')
# HLSmod.comp.2018 <- lapply(planet.path, function(x) HLS.MODIS(x,'2018',8,8))
# saveRDS(HLSmod.comp.2018, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_85GD_2018_wQA.RDS')
# HLSmod.comp.2019 <- lapply(planet.path, function(x) HLS.MODIS(x,'2019',8,8))
# saveRDS(HLSmod.comp.2019, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_85GD_2019_wQA.RDS')

# EVI amp comparison
# HLSmod.comp.2017 <- lapply(planet.path, function(x) HLS.MODIS(x,'2017',10,10))
# saveRDS(HLSmod.comp.2017, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_EVIamp_2017.RDS')
# HLSmod.comp.2018 <- lapply(planet.path, function(x) HLS.MODIS(x,'2018',10,10))
# saveRDS(HLSmod.comp.2018, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_EVIamp_2018.RDS')
# HLSmod.comp.2019 <- lapply(planet.path, function(x) HLS.MODIS(x,'2019',10,10))
# saveRDS(HLSmod.comp.2019, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_EVIamp_2019.RDS')

# EVI area comparison
# HLSmod.comp.2017 <- lapply(planet.path, function(x) HLS.MODIS(x,'2017',11,11))
# saveRDS(HLSmod.comp.2017, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_EVIarea_2017.RDS')
# HLSmod.comp.2018 <- lapply(planet.path, function(x) HLS.MODIS(x,'2018',11,11))
# saveRDS(HLSmod.comp.2018, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_EVIarea_2018.RDS')
# HLSmod.comp.2019 <- lapply(planet.path, function(x) HLS.MODIS(x,'2019',11,11))
# saveRDS(HLSmod.comp.2019, '/projectnb/modislc/users/sjstone/planet/data/upscale/MODIS_HLS_median_EVIarea_2019.RDS')


