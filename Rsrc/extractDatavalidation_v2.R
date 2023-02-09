library(data.table)
library(raster)
library(sf)

valData <- fread("/scratch/project_2000994/PREBASruns/assessCarbon/data/valData.txt",header = T)
valData_shp <- shapefile("/scratch/project_2000994/PREBASruns/assessCarbon/data/valShp/AC_2019_accuracy_plots.shp")

tiles <- c("35VLJ", "34VEQ", "35WMN")



dataAll <- list()
for(i in 1:length(tiles)){
  tileX <- tiles[i]
  pathX <- paste0("/scratch/project_2000994/PREBASruns/assessCarbon/rasters/Finland/AC_training_FI_",tileX,"/outRast/init2016/")
  rastX <- raster(paste0(pathX,"Bda2019.tif"))
  dataAll[[i]] <- crop(valData_shp,extent(rastX))
###Ba 
  rastX <- raster(paste0(pathX,"Bda2019.tif"))
  dataAll[[i]]$Bda2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  rastX <- raster(paste0(pathX,"Bm2019.tif"))
  dataAll[[i]]$Bm2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  rastX <- raster(paste0(pathX,"Bs2019.tif"))
  dataAll[[i]]$Bs2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  
  plot(dataAll[[i]]$Bs2019,dataAll[[i]]$G,ylab="measuremens",xlab="models",main="B")
  points(dataAll[[i]]$Bm2019,dataAll[[i]]$G,col=2,pch=20)
  points(dataAll[[i]]$Bda2019,dataAll[[i]]$G,col=3,pch=20)
  
  #H 
  rastX <- raster(paste0(pathX,"Hda2019.tif"))
  dataAll[[i]]$Hda2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  rastX <- raster(paste0(pathX,"Hm2019.tif"))
  dataAll[[i]]$Hm2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  rastX <- raster(paste0(pathX,"Hs2019.tif"))
  dataAll[[i]]$Hs2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  
  plot(dataAll[[i]]$Hs2019,dataAll[[i]]$H,ylab="measuremens",xlab="models",main="H")
  points(dataAll[[i]]$Hm2019,dataAll[[i]]$H,col=2,pch=20)
  points(dataAll[[i]]$Hda2019,dataAll[[i]]$H,col=3,pch=20)
  
  #D 
  rastX <- raster(paste0(pathX,"Dda2019.tif"))
  dataAll[[i]]$Dda2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  rastX <- raster(paste0(pathX,"Dm2019.tif"))
  dataAll[[i]]$Dm2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  rastX <- raster(paste0(pathX,"Ds2019.tif"))
  dataAll[[i]]$Ds2019 <-raster::extract(rastX, dataAll[[i]], fun = mean, na.rm = TRUE)
  
  plot(dataAll[[i]]$Ds2019,dataAll[[i]]$D,ylab="measuremens",xlab="models",main="D")
  points(dataAll[[i]]$Dm2019,dataAll[[i]]$D,col=2,pch=20)
  points(dataAll[[i]]$Dda2019,dataAll[[i]]$D,col=3,pch=20)
  
}


dataRes2019 <- as.data.table(dataAll[[1]])
dataRes2019$Tile <- tiles[1]
dataX <- as.data.table(dataAll[[2]])
dataX$Tile <- tiles[2]
dataRes2019 <- rbind(dataRes2019,dataX)
dataX <- as.data.table(dataAll[[3]])
dataX$Tile <- tiles[3]
dataRes2019 <- rbind(dataRes2019,dataX)

save(dataRes2019,file="/scratch/project_2000994/PREBASruns/assessCarbon/dataRes2019.rdata")