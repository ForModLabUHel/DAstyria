source("localSettings.r")
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")

setwd(generalPath)

shapeRastOr <- raster("2015/FCM_STY_2015_G_10M_8BITS-Styria.tif")
if(climData == "cliPick"){
  shapeRast <- aggregate(shapeRastOr,fact=500)
  projection(shapeRast)
  
  myRaster <- projectRaster(shapeRast,crs =  "+proj=longlat +datum=WGS84")
  plot(myRaster)
  
  # extract the Altri weather data from clipick from 1.1.2000-2.1.2000 at 10 km grid
  wDs <- getWD("Styria", myRaster, 10, startingYear, 1, 1, yearEnd, 12, 31)
  
  # Find the unique data tables
  # weather data for the site is here 
  weather_data <- unique(wDs)
  
  length(wDs)
  length(weather_data)
  
  
  # raster of climID:s is here
  climID_raster <- climRaster(weather_data, wDs, "Styria", myRaster, 10000)
  
  ###Converts in the orginal coordinate system of the input raster or shape File
  climID_raster <- projectRaster(climID_raster,
                                 crs = crs(shapeRastOr),
                                 res = res(shapeRastOr),
                                 method = "ngb") 
  
  climID_raster <- crop(climID_raster,extent(shapeRastOr))
  climID_raster <- raster(vals=values(climID_raster),ext=extent(shapeRastOr),crs=crs(shapeRastOr),
                          nrows=dim(shapeRastOr)[1],ncols=dim(shapeRastOr)[2])
  
  climID_raster <- mask(climID_raster,shapeRastOr)
  
  writeRaster(climID_raster,filename = "weatherInputs/climID_cliPick.tif",overwrite=T)
  
  # # weather data as one data.table, climID as "id"
  # weather_data_table <- rbindlist(weather_data, idcol="id")
  
  prebasWeather(weather_data)
  
}

if(climData == "eObs"){
  weatherDataBase <- fread("weatherInputs/prebas_styria_1992-2021_eObs.csv")
  weatherDataBase[,year:=year(time)]
  weatherDataBase[,climID:=climID+1]
  climIDs <- unique(weatherDataBase[,.(lon,lat,climID)])
  
  # raster of climID:s is here
  climID_raster <- rasterFromXYZ(climIDs,crs =  "+proj=longlat +datum=WGS84")
  
  # shapeRast <- aggregate(shapeRastOr,fact=500)
  # projection(shapeRast)
  # 
  
  ###Converts in the orginal coordinate system of the input raster or shape File
  climID_raster <- projectRaster(climID_raster,
                                 crs = crs(shapeRastOr),
                                 res = res(shapeRastOr),
                                 method = "ngb") 
  
  climID_raster <- crop(climID_raster,extent(shapeRastOr))
  climID_raster <- raster(vals=values(climID_raster),ext=extent(shapeRastOr),crs=crs(shapeRastOr),
                          nrows=dim(shapeRastOr)[1],ncols=dim(shapeRastOr)[2])
  
  climID_raster <- mask(climID_raster,shapeRastOr)
  climIDs_Styria <- unique(getValues(climID_raster))
  climIDs_Styria <- climIDs_Styria[!is.na(climIDs_Styria)]
  climIDs_StyriaOr <- sort(climIDs_Styria)
  climIDs_StyriaNew <- 1:length(climIDs_StyriaOr)
  
  df <- data.frame(idOr=climIDs_StyriaOr, idNew=climIDs_StyriaNew)
  climID_raster <- subs(climID_raster, df)
  
  weatherStyria <- weatherDataBase[climID %in% climIDs_StyriaOr]
  setnames(weatherStyria,"climID","climID_or")
  
  weatherStyria$climID <- match(weatherStyria$climID_or,df$idOr)
  
  writeRaster(climID_raster,filename = "weatherInputs/climID_eObs.tif",overwrite=T)
  
  # # weather data as one data.table, climID as "id"
  # weather_data_table <- rbindlist(weather_data, idcol="id")
  setkey(weatherStyria,climID,year,DOY)
  weatherStyria <- weatherStyria[DOY!=366]
  ndays <-  nrow(weatherStyria[climID==1])
  nClimIDs <- length(climIDs_StyriaNew)
  PAR <- matrix(weatherStyria$PAR,ncol = ndays,nrow = nClimIDs,byrow = T)
  Precip <- matrix(weatherStyria$Precip,ncol = ndays,nrow = nClimIDs,byrow = T)
  TAir <- matrix(weatherStyria$TAir,ncol = ndays,nrow = nClimIDs,byrow = T)
  VPD <- matrix(weatherStyria$VPD,ncol = ndays,nrow = nClimIDs,byrow = T)
  CO2 <- matrix(380,ncol = ndays,nrow = nClimIDs,byrow = T)
  
  save(PAR,Precip,VPD,TAir,CO2,file=paste0("weatherInputs/weather_",climData,".rdata"))
  
}