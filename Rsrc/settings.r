###choose PREBAS version
vPREBAS <- "master"   #### choose PREBAS version to run the model  "master"


#####Settings####
if(!exists("testRun")) testRun = T ####set to TRUE to test the code on a small raster proportion
if(!exists("extNew")) extNew <- c(550000,550500,5200000,5200500)
if(!exists("CSCrun")){
  CSCrun = F ### set to TRUE if you are running on CSC
}
fracTest <- 0.2 ###fraction of test area
maxSitesRun <- 20000
maxSitesRunTest <- 20000
saveVars <- c(1,11:13,17,30,43,44) ####select variables to save
varHD <- FALSE #### if true will vary H and D of pine and spruce using siteType

###library path in CSC project_2000994
if(CSCrun){
  .libPaths(c("/projappl/project_2000994/project_rpackages", .libPaths()))
  libpath <- .libPaths()[1]
}

##load libraries
library(lfda)
library(mvtnorm)
library(reshape2)
library(plyr)
library(data.table)
require(sm)
require(rgdal)
library(raster)
library(parallel)
library(MASS)
library(readxl)
library(minpack.lm)
library(sf)
library(fasterize)
library(abind)


###check prebas version and install if needed
if(!CSCrun){
  devtools::install_github("ForModLabUHel/Rprebasso", ref=vPREBAS)
}
require(Rprebasso)

####indicate rasterPath and climID path
if(CSCrun){
  generalPath <- "/scratch/project_2000994/PREBASruns/FCMaustria/"
}else{
  generalPath <- "C:/Users/minunno/Documents/research/ForestCarbonMonitoring/FCM_CNN/data/Forest_Structure_layers/"
}

rasterPath <- generalPath

rasterPath2015 <- paste0(generalPath,"2015/")
rasterPath2018 <- paste0(generalPath,"2018/")
rasterPath2021 <- paste0(generalPath,"2021/")
procDataPath <- paste0(generalPath,"procData/")
outPath <- paste0(generalPath,"output/")
initPrebasPath <- paste0(generalPath,"initPrebas/")


if(!exists("climData")) climData = "eObs"
if(CSCrun){
  if(!exists("climIDpath")) climIDpath <- paste0("weatherInputs/climID_",climData,".tif")
  if(!exists("climatepath")) climatepath = paste0("weatherInputs/weather_",climData,".rdata") #### local fm
}else{
  if(!exists("climatepath")) climatepath = "C:/Users/minunno/Documents/research/extarctWeather/inputs/" #### local fm
  if(!exists("climIDpath")) climIDpath <- "C:/Users/minunno/Documents/research/FinSeg/some stuff/climID10km.tif"
}

startYearWeather <- 1971 ###1971 for Finnish weather dataBase
startingYear <- 2015  #2019
year2 <- 2018 ###year of the second measurement
yearEnd <- 2021     #2024


####indicate raster files
# tileX = "34VEQ"
# areaID <- "FI"
year <- 2015
baRast2015 <-  paste0(rasterPath,year,"/FCM_STY_2015_G_10M_8BITS-Styria.tif")
blPerRast2015 <- paste0(rasterPath,year,"/FCM_STY_2015_BLP_10M_8BITS-Styria.tif")
dbhRast2015 <- paste0(rasterPath,year,"/FCM_STY_2015_D_10M_8BITS-Styria.tif")
vRast2015 <- paste0(rasterPath,year,"/FCM_STY_2015_V_10M_16BITS-Styria.tif")
hRast2015 <- paste0(rasterPath,year,"/FCM_STY_2015_H_10M_16BITS-Styria.tif")
conifPerRast2015 <- paste0(rasterPath,year,"/FCM_STY_2015_CP_10M_8BITS-Styria.tif")
siteTypeRast2015 <- paste0(rasterPath,year,"/FCM_STY_2015_SITE_10M_8BITS-Styria.tif")
# pinePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_pine_10M_1CHS_8BITS.tif")
# sprucePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_spruce_10M_1CHS_8BITS.tif")
year <- 2018
baRast2018 <-  paste0(rasterPath,year,"/FCM_STY_2018_G_10M_8BITS-Styria.tif")
blPerRast2018 <- paste0(rasterPath,year,"/FCM_STY_2018_BLP_10M_8BITS-Styria.tif")
dbhRast2018 <- paste0(rasterPath,year,"/FCM_STY_2018_D_10M_8BITS-Styria.tif")
vRast2018 <- paste0(rasterPath,year,"/FCM_STY_2018_V_10M_16BITS-Styria.tif")
hRast2018 <- paste0(rasterPath,year,"/FCM_STY_2018_H_10M_16BITS-Styria.tif")
conifPerRast2018 <- paste0(rasterPath,year,"/FCM_STY_2018_CP_10M_8BITS-Styria.tif")
siteTypeRast2018 <- paste0(rasterPath,year,"/FCM_STY_2018_SITE_10M_8BITS-Styria.tif")
# pinePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_pine_10M_1CHS_8BITS.tif")
# sprucePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_spruce_10M_1CHS_8BITS.tif")
year <- 2021
baRast2021 <-  paste0(rasterPath,year,"/FCM_STY_2021_G_10M_8BITS-Styria.tif")
blPerRast2021 <- paste0(rasterPath,year,"/FCM_STY_2021_BLP_10M_8BITS-Styria.tif")
dbhRast2021 <- paste0(rasterPath,year,"/FCM_STY_2021_D_10M_8BITS-Styria.tif")
vRast2021 <- paste0(rasterPath,year,"/FCM_STY_2021_V_10M_16BITS-Styria.tif")
hRast2021 <- paste0(rasterPath,year,"/FCM_STY_2021_H_10M_16BITS-Styria.tif")
conifPerRast2021 <- paste0(rasterPath,year,"/FCM_STY_2021_CP_10M_8BITS-Styria.tif")
siteTypeRast2021 <- paste0(rasterPath,year,"/FCM_STY_2021_SITE_10M_8BITS-Styria.tif")
# pinePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_pine_10M_1CHS_8BITS.tif")
# sprucePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_spruce_10M_1CHS_8BITS.tif")

# Source of tile-specific settings. Defined in batch job script. When set to TRUE will overwrite the tile-specific 
# settings in this script (lines: 41-49, 53-56, 72-81) with settings from filepath in mySettings variable.
# if(exists("tileSettings")){
#   if(tileSettings) {
#     source(mySettings)
#     
#     # Indicate raster files
#     baRast <-  paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_BA_10M_1CHS_8BITS.tif")
#     blPerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_BLP_10M_1CHS_8BITS.tif")
#     dbhRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_DIA_10M_1CHS_8BITS.tif")
#     vRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_GSV_10M_1CHS_16BITS.tif")
#     hRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_HGT_10M_1CHS_16BITS.tif")
#     pinePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_pine_10M_1CHS_8BITS.tif")
#     sprucePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_spruce_10M_1CHS_8BITS.tif")
#     siteTypeRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_SITE_10M_1CHS_8BITS.tif")
#     siteTypeRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_SITE_10M_1CHS_8BITS.tif")
#     vRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_GSV_10M_1CHS_16BITS.tif")
#     baRast2 <-  paste0(rasterPath,areaID,"_",tileX,"-",year2,"_BA_10M_1CHS_8BITS.tif")
#     dbhRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_DIA_10M_1CHS_8BITS.tif")
#     hRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_HGT_10M_1CHS_16BITS.tif")
#     pinePerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_P_pine_10M_1CHS_8BITS.tif")
#     sprucePerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_P_spruce_10M_1CHS_8BITS.tif")
#     blPerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_BLP_10M_1CHS_8BITS.tif")
#     mgmtmaskRast <- paste0(rasterPath, areaID, "_", tileX, "_mgmtmask.tif")
#   }
# }


nYears <-  yearEnd - startingYear ## number of simulation years
domSPrun = 0.
mgmtmask = F # switch for masking of management


resX <- 10 ### pixel resolution in meters

### define weather inputs (CurrClim, or climate models)
weather = "CurrClim"

###set harvests
defaultThin = 0.
ClCut = 0.
harvscen = "NoHarv"

####set values for NAs and convert factor for prebas units
baNA <- c(253:255); baConv<- 1
blPerNA <- c(253:255); blPerConv<- 0.01
dbhNA <- c(253:255); dbhConv <- 1
vNA <- c(65533:65535); vConv <- 1
hNA <- c(65533:65535); hConv <- 0.1
pinePerNA <- c(253:255); pinePerConv <- 0.01
sprucePerNA <- c(253:255); sprucePerConv <- 0.01
siteTypeNA <- c(253:255); siteTypeConv <- 1

####settings for sitetype estimation
stXruns <- TRUE
siteTypeX <- startingYear #startingYear #year2 #startingYear #1:5


# Set TRUE to enable running 1.8_optST, 2_InitPreb and 3_runModel in parallel. Set to FALSE, these scripts run as serial.
parallelRun <- FALSE
coresN <- 20L ### Set number of cores to use in parallel run 

# Set whether to split unique data in 1.1 to smaller parts. If
# TRUE, data is split.
splitRun <- FALSE
nSplit <- 20
if(splitRun){    # Range/number of split parts
  splitRange <- 1:nSplit
}

####thresholds for variables to reset stand from plantation
maxDens <- 10000
initH <- 1.5
initDBH <- 0.5
initN <- 2200
initBA <- pi*(initDBH/200)^2*initN

#####settings for data extraction
varDT <- c(11:13,30)   ####variables to extract in DT
extrFun <- c("mean","mean","sum","sum")
layerDT <- "all" ###layerID to report in data.tables, if layerDT==all the all layers are considered as sum or mean according to extrFun

#####settings for raster creation
varRast <- varDT  #c(44,30)   ####variables to extract in DT
yearOut <- yearEnd#c(2019,2024)

#####filter model output raster
minX <- c(0,0,0,0)
maxX <- c(70,40,60,1000)


#####end Settings####
