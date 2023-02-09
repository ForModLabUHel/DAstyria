library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file in local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")

# Check and create output directories
setwd(generalPath)
mkfldr <- paste0("outRast/","init",startingYear)
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}
 
load(paste0(procDataPath,"init",startingYear,"/","DA",yearOut,"/XYsegID.rdata"))  
crsX <- crs(raster(baRast))

setkey(XYsegID,segID)
  Dda2019 <- Dm2019 <- Ds2019 <- 
  Hda2019 <- Hm2019 <- Hs2019 <- 
  Bda2019 <- Bm2019 <- Bs2019 <- data.table()
for(i in 1:20){
  load(paste0("posterior/pMvn_FSV_split",i,".rdata"))
  pMvNorm$varNam <- rep(
    c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
      "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
      "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
    times = nrow(pMvNorm)/126)
  setkey(pMvNorm,segID)
  Dda2019 <- rbind(Dda2019,merge(XYsegID,pMvNorm[varNam=="DDA2019"]))
  Dm2019 <- rbind(Dm2019,merge(XYsegID,pMvNorm[varNam=="Dm2019"]))
  Ds2019 <- rbind(Ds2019,merge(XYsegID,pMvNorm[varNam=="Ds2019"]))
  Hda2019 <- rbind(Hda2019,merge(XYsegID,pMvNorm[varNam=="HDA2019"]))
  Hm2019 <- rbind(Hm2019,merge(XYsegID,pMvNorm[varNam=="Hm2019"]))
  Hs2019 <- rbind(Hs2019,merge(XYsegID,pMvNorm[varNam=="Hs2019"]))
  Bda2019 <- rbind(Bda2019,merge(XYsegID,pMvNorm[varNam=="BDA2019"]))
  Bm2019 <- rbind(Bm2019,merge(XYsegID,pMvNorm[varNam=="Bm2019"]))
  Bs2019 <- rbind(Bs2019,merge(XYsegID,pMvNorm[varNam=="Bs2019"]))
  print(i)
}

outXs <- c("Dda2019", "Dm2019", "Ds2019",
                 "Hda2019", "Hm2019", "Hs2019", 
                 "Bda2019", "Bm2019", "Bs2019")
for(outX in outXs){
  rastX <- rasterFromXYZ(get(outX)[,c("x","y","V1"),with=F])
  crs(rastX) <- crsX
  rastName <- paste0("outRast/","init",startingYear,"/",outX,".tif")
  writeRaster(rastX,filename = rastName,overwrite=T)
}