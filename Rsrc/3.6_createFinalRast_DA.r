if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists

# Run settings 
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/functions.r")
setwd(generalPath)


baRast2015 <- "finalRast/FCM_STY_2015_G_10M_8BITS-Styria.tif"
blPerRast2015 <- "finalRast/FCM_STY_2015_BLP_10M_8BITS-Styria.tif"
dbhRast2015 <- "finalRast/FCM_STY_2015_D_10M_8BITS-Styria.tif"
hRast2015 <- "finalRast/FCM_STY_2015_H_10M_16BITS-Styria.tif"
conifPerRast2015 <- "finalRast/FCM_STY_2015_CP_10M_8BITS-Styria.tif"

baRast2018 <- "finalRast/FCM_STY_2018_G_10M_8BITS-Styria.tif"
blPerRast2018 <- "finalRast/FCM_STY_2018_BLP_10M_8BITS-Styria.tif"
dbhRast2018 <- "finalRast/FCM_STY_2018_D_10M_8BITS-Styria.tif"
hRast2018 <- "finalRast/FCM_STY_2018_H_10M_16BITS-Styria.tif"
conifPerRast2018 <- "finalRast/FCM_STY_2018_CP_10M_8BITS-Styria.tif"

baRast2021 <- "finalRast/FCM_STY_2021_G_10M_8BITS-Styria.tif"
blPerRast2021 <- "finalRast/FCM_STY_2021_BLP_10M_8BITS-Styria.tif"
dbhRast2021 <- "finalRast/FCM_STY_2021_D_10M_8BITS-Styria.tif"
hRast2021 <- "finalRast/FCM_STY_2021_H_10M_16BITS-Styria.tif"
conifPerRast2021 <- "finalRast/FCM_STY_2021_CP_10M_8BITS-Styria.tif"


if(cal=="austria"){
  load("/scratch/project_2000994/calibrations/calAustria/outCal/pCROBASaustria.rdata")
  pCROBAS <- pCROBASaustria
}else{
  pCROBAS <- pCROB
}
load(paste0("surErrMods/surMod_Step1_cal",cal,".rdata"))
load(paste0("surErrMods/surMod_Step2_cal",cal,".rdata"))

# source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/rmvweisd.r")
# source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/Fweibull_Arithmetic Mean and Variance.R")

# if(modifiedSettings) {
#   source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings.r") # in CSC
# }

# if (splitRun) {
#   print(paste("# of splits", nSplit))
#   # If output is set to be split to smaller parts (splitRun = TRUE), create separate
#   # folder for the split data tables.
#   mkfldr_split <- "procDataFinal/split"
#   if(!dir.exists(file.path(generalPath, mkfldr_split))) {
#     dir.create(file.path(generalPath, mkfldr_split), recursive = TRUE)
#   }
# } 
# 
# mkfldr <- "procDataFinal"
# if(!dir.exists(file.path(generalPath, mkfldr))) {
#   dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
# }

###extract CurrClim IDs
rastRef <- raster(baRast2015)
if(testRun){
  extNew <- extent(rastRef)
  extNew[1:2]   <- (extent(rastRef)[1] + (extent(rastRef)[2]))/2 + c(10000,30000)
  extNew[3:4]   <- (extent(rastRef)[3] + (extent(rastRef)[4]))/2 + c(-10000,10000)
  rastRef <- crop(rastRef,extNew)
  maxSitesRun <- maxSitesRunTest
}
print("loading")
load(file=paste0("procData/splitFinal/finalRasterData_id",i,".rdata"))
print("loaded")
#      
# load(file="procData/finalRasterData.rdata")
# 
# splitRange <- splitIDsRaster
# split_length <- ceiling(nrow(data.all)/length(splitRange))
# data.all <- data.all[, split_id := NA]
# 
# 
# data.all$split_id[1:split_length] <- 1
# for (i in 2:(max(splitRange)-1)) {
#   data.all$split_id[((i-1)*split_length+1):(split_length*i)] <- i
# }
# data.all$split_id[((length(splitRange)-1)*split_length+1):(nrow(data.all))] <- length(splitRange)
# 
print("computing")

system.time({
  ops <- calcVW(dataX,yearRast,pCROBAS,
                step.modelSVIx = step.modelSVIx)
})

print("calculations completed")

###convert biomasses from kgC/ha to tons(DM)/ha
ops$abgWprebas <- ops$abgWprebas *2 / 1000
ops$blgWprebas <- ops$blgWprebas *2 / 1000

dataX <- ops[,.(x,y,SVIpreb,WblgPreb,WabgPreb,Vpreb)]
save(dataX,file="procData/splitFinal/data_id",splitID,".rdata")

SVIpreb <- rasterFromXYZ(ops[,.(x,y,SVIprebas)],crs=crs(rastRef))
WblgPreb <- rasterFromXYZ(ops[,.(x,y,blgWprebas)],crs=crs(rastRef))
WabgPreb <- rasterFromXYZ(ops[,.(x,y,abgWprebas)],crs=crs(rastRef))
Vpreb <- rasterFromXYZ(ops[,.(x,y,Vprebas)],crs=crs(rastRef))
# Save split tables

writeRaster(SVIpreb,filename = paste0("finalRast/SVI_prebas_",yearRast,splitID,".tif"),overwrite=TRUE)
writeRaster(WblgPreb,filename = paste0("finalRast/Wbl_prebas_",yearRast,splitID,".tif"),overwrite=TRUE)
writeRaster(WabgPreb,filename = paste0("finalRast/Wabg_prebas_",yearRast,splitID,".tif"),overwrite=TRUE)
writeRaster(Vpreb,filename = paste0("finalRast/V_prebas_",yearRast,splitID,".tif"),overwrite=TRUE)


