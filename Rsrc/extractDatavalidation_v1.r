library(data.table)
library(raster)
setwd("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/batchFiles")
source("runSettings.r")
load("/scratch/project_2000994/PREBASruns/assessCarbon/data/traningSites.rdata")
# load("/scratch/project_2000994/PREBASruns/assessCarbon/data/")

tileX <- "34VEQ"
setwd(paste0("/scratch/project_2000994/PREBASruns/assessCarbon/rasters/Finland/AC_training_FI_",tileX,"/"))
load("procData/init2016/DA2019/XYsegID.rdata") 
rastX <- rasterFromXYZ(XYsegID)
  

# tileX <- "35VLJ"
data2019x <- data2019[S2Tile==tileX]
segIDs <- extract(rastX,data2019x[,.(XCOORD,YCOORD)])
data2019x$segID <- segIDs

nSplit <- 20
data2019res <- data.table()
###join data
for(i in 1:nSplit){
  load(paste0("posterior/pMvn_FSV_split",i,".rdata"))
  pMvNorm <- data.table(pMvNorm)
  pMvNorm$varNam <- rep(
    c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
      "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
      "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
    times = nrow(pMvNorm)/126)

  dataX <- data2019x[segID %in% unique(pMvNorm$segID)]
  # if (i==13) 
    # dataX <- dataX[!duplicated(segID)]
  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "Bm2019"]$V1
  dataX$Bm2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 
  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "Bs2019"]$V1
  dataX$Bs2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 
  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "BDA2019"]$V1
  dataX$BDA2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 

  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "Dm2019"]$V1
  dataX$Dm2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 
  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "Ds2019"]$V1
  dataX$Ds2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 
  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "DDA2019"]$V1
  dataX$DDA2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 

  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "Hm2019"]$V1
  dataX$Hm2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 
  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "Hs2019"]$V1
  dataX$Hs2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 
  xxx <- pMvNorm[segID %in% dataX$segID & pMvNorm$varNam == "HDA2019"]$V1
  dataX$HDA2019 <- xxx[match(dataX$segID,unique(dataX$segID))] 

  data2019res <- rbind(data2019res,dataX)
  rm(pMvNorm,dataX); gc()
  print(i) 
}
save(data2019res,file=paste0("data2019res_",tileX,".rdata"))