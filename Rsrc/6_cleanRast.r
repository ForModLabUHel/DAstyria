library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file in local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")
setwd(generalPath)

# raster
# load(paste0(procDataPath,startingYear,"/XYsegID.rdata"))
# crsX <- crs(raster(baRast))
# 
clims <- weather
mans <- harvscen
# 
for(ij in yearOut){
  for(varX in 1:length(varRast)){
    rastName <- paste0("outRast/","init",startingYear,"/st",siteTypeX,"/",clims,"_",mans,"_var",varNames[varRast[varX]],
                       "_spec",layerDT,"_yearStart",startingYear,"_yearOut",ij,".tif")
    rastX <- raster(rastName)
    # print(length(which(!is.na(getValues(rastX)))))
    # print(length(which(getValues(rastX)>maxX[varX])))
    rastX[rastX < minX[varX]] <- NA
    rastX[rastX > maxX[varX]] <- NA
    # print(length(which(!is.na(getValues(rastX)))))
    rastNameCleaned <- paste0("outRast/","init",startingYear,"/st",siteTypeX,"/",clims,"_",mans,"_var",varNames[varRast[varX]],
                       "_spec",layerDT,"_yearStart",startingYear,"_yearOut",ij,"_cleaned.tif")
    writeRaster(rastX,filename = rastNameCleaned,overwrite=T)
  }
  
  
  # rastName <- paste0("outRast/",startingYear,"/","startV_","startYear",startingYear,
  #                    "_layer",layerDT,".tif")
  # rastX <- raster(rastName)
  # print(length(which(!is.na(getValues(rastX)))))
  # # print(length(which(getValues(rastX)>200)))
  # # rastX[rastX < minX[varX]] <- NA
  # # rastX[rastX > maxX[varX]] <- NA
  # print(length(which(!is.na(getValues(rastX)))))
  # 
}


