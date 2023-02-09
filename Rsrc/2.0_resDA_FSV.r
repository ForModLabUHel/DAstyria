library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")

###check and create output directories
###check and create output directories
setwd(generalPath)
mkfldr <- paste0("outRast/",paste0("init",startingYear,"/DA",year2))
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

# yearX <- 3
# nSample = 1000 ###number of samples from the error distribution


###join data
for(i in 1:nSplit){
  load(paste0("posterior/pMvn_FSV_split",i,".rdata"))
  pMvNorm <- data.table(pMvNorm)
  # pMvNormAll <- pMvNorm
  # load("pMvn_FSV_split2.rdata")
  # pMvNorm <- data.table(pMvNorm)
  # pMvNormAll <- rbind(pMvNormAll,pMvNorm)
  
  pMvNorm$varNam <- rep(
    c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
      "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
      "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
    times = nrow(pMvNorm)/126)
  
  pMvNorm <- pMvNorm[!(varNam=="varcov1" | varNam=="varcov2" | varNam=="varcov3")]
  
  
  dataX <- data.table(dcast(data = pMvNorm,
                              formula = segID~varNam,value.var = "V1"))
  if(i ==  1) dataAll <- dataX
  if(i>1) dataAll <- rbind(dataAll,dataX)
  
  rm(pMvNorm,dataX); gc()
  print(i) 
}
print(range(dataAll$Hm2019))

load(paste0(procDataPath,"init",startingYear,"/DA",year2,"/XYsegID.rdata"))  

setkey(XYsegID,segID)
setkey(dataAll,segID)
outXY <- merge(XYsegID,dataAll,all = T)


crsX <- crs(raster(baRast))

  ###remove coordinates NAs
  outXY <- outXY[!is.na(x)]
  
  ###create rasters
  vars <- names(outXY)
  vars <- vars[!vars %in% c("x","y","segID")]
  for(varX in vars){
    rastX <- rasterFromXYZ(outXY[,c("x","y",varX),with=F])
    crs(rastX) <- crsX
    
    rastName <- paste0("outRast/","init",startingYear,"/DA",year2,"/",varX,".tif")
    
    writeRaster(rastX,filename = rastName,overwrite=T)
    print(varX)
  }
