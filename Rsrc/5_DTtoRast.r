library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file in local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")

# Check and create output directories
setwd(generalPath)
mkfldr <- paste0("outRast/","init",startingYear,"/st",siteTypeX)
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}
 
load(paste0(procDataPath,"init",startingYear,"/","DA",yearOut,"/XYsegID.rdata"))  
crsX <- crs(raster(baRast))
rastX <- rasterFromXYZ(XYsegID)
segIDx <- extract(rastX,data2019[S2Tile==tileX,.(XCOORD,YCOORD)])

plot(rastX)
points(data2019[S2Tile==tileX,.(XCOORD,YCOORD)],pch=20,col=2)

setwd(generalPath)


load("posterior/pMvn_FSV_split1.rdata")

load("/scratch/project_2000994/PREBASruns/assessCarbon/data/traningSites.rdata")
segIDx <- extract(rastX,data2019[S2Tile==tileX,.(XCOORD,YCOORD)])
data2019[S2Tile==tileX,segID:=segIDx]
setkey(data2019,segID)
setkey(pMvNorm,segID)

pMvNorm$varNam <- rep(
  c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
    "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
    "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
  times = nrow(pMvNorm)/126)

oo <- merge(data2019, pMvNorm[varNam=="Hm2019",1:2], by.x = "segID", 
      by.y = "segID", all.x = TRUE, all.y = FALSE)

setkey(XYsegID,segID)
  Dda2019 <- Dm2019 <- Ds2019 <- 
  Hda2019 <- Hm2019 <- Hs2019 <- 
  Bda2019 <- Bm2019 <- Bs2019 <- data.table()
for(i in 1:20){
  load(paste0("posterior/pMvn_FSV_split",i,".rdata"))
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
}
for(i in 2:20){
  load(paste0("posterior/pMvn_FSV_split",i,".rdata"))
  setkey(oo,segID)
  setkey(pMvNorm,segID)
       
  pMvNorm$varNam <- rep(
         c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
           "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
           "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
         times = nrow(pMvNorm)/126)
       
  assign(paste0("Hm2019_",i), merge(data2019, pMvNorm[varNam=="Hm2019",1:2], by.x = "segID", 
                   by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("HDA2019_",i), merge(data2019, pMvNorm[varNam=="HDA2019",1:2], by.x = "segID", 
                               by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("Hs2019_",i), merge(data2019, pMvNorm[varNam=="Hs2019",1:2], by.x = "segID", 
                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("Bm2019_",i), merge(data2019, pMvNorm[varNam=="Bm2019",1:2], by.x = "segID", 
                                    by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("BDA2019_",i), merge(data2019, pMvNorm[varNam=="BDA2019",1:2], by.x = "segID", 
                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("Bs2019_",i), merge(data2019, pMvNorm[varNam=="Bs2019",1:2], by.x = "segID", 
                                    by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("Dm2019_",i), merge(data2019, pMvNorm[varNam=="Dm2019",1:2], by.x = "segID", 
                                    by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("DDA2019_",i), merge(data2019, pMvNorm[varNam=="DDA2019",1:2], by.x = "segID", 
                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
  assign(paste0("Ds2019_",i), merge(data2019, pMvNorm[varNam=="Ds2019",1:2], by.x = "segID", 
                                    by.y = "segID", all.x = FALSE, all.y = FALSE))
  print(i)
}




data2019res <- rbind(Ds2019_10,Ds2019_11)
data2019res <- rbind(data2019res,Ds2019_12)
data2019res <- rbind(data2019res,Ds2019_13)
data2019res <- rbind(data2019res,Ds2019_14)
data2019res <- rbind(data2019res,Ds2019_15)
data2019res <- rbind(data2019res,Ds2019_16)
data2019res <- rbind(data2019res,Ds2019_17)
data2019res <- rbind(data2019res,Ds2019_18)
data2019res <- rbind(data2019res,Ds2019_19)
data2019res <- rbind(data2019res,Ds2019_20)

setnames(data2019res,"V1","Ds2019")

data2019res$DDA2019 <- c(DDA2019_10$V1,DDA2019_11$V1,DDA2019_12$V1,DDA2019_13$V1,DDA2019_14$V1,DDA2019_15$V1,
                         DDA2019_16$V1,DDA2019_17$V1,DDA2019_18$V1,DDA2019_19$V1,DDA2019_20$V1)

data2019res$Dm2019 <- c(Dm2019_10$V1,Dm2019_11$V1,Dm2019_12$V1,Dm2019_13$V1,Dm2019_14$V1,Dm2019_15$V1,
                         Dm2019_16$V1,Dm2019_17$V1,Dm2019_18$V1,Dm2019_19$V1,Dm2019_20$V1)



data2019res$Hm2019 <- c(Hm2019_10$V1,Hm2019_11$V1,Hm2019_12$V1,Hm2019_13$V1,Hm2019_14$V1,Hm2019_15$V1,
                        Hm2019_16$V1,Hm2019_17$V1,Hm2019_18$V1,Hm2019_19$V1,Hm2019_20$V1)

data2019res$Hs2019 <- c(Hs2019_10$V1,Hs2019_11$V1,Hs2019_12$V1,Hs2019_13$V1,Hs2019_14$V1,Hs2019_15$V1,
                        Hs2019_16$V1,Hs2019_17$V1,Hs2019_18$V1,Hs2019_19$V1,Hs2019_20$V1)

data2019res$HDA2019 <- c(HDA2019_10$V1,HDA2019_11$V1,HDA2019_12$V1,HDA2019_13$V1,HDA2019_14$V1,HDA2019_15$V1,
                        HDA2019_16$V1,HDA2019_17$V1,HDA2019_18$V1,HDA2019_19$V1,HDA2019_20$V1)


data2019res$Bm2019 <- c(Bm2019_10$V1,Bm2019_11$V1,Bm2019_12$V1,Bm2019_13$V1,Bm2019_14$V1,Bm2019_15$V1,
                        Bm2019_16$V1,Bm2019_17$V1,Bm2019_18$V1,Bm2019_19$V1,Bm2019_20$V1)

data2019res$Bs2019 <- c(Bs2019_10$V1,Bs2019_11$V1,Bs2019_12$V1,Bs2019_13$V1,Bs2019_14$V1,Bs2019_15$V1,
                        Bs2019_16$V1,Bs2019_17$V1,Bs2019_18$V1,Bs2019_19$V1,Bs2019_20$V1)

data2019res$BDA2019 <- c(BDA2019_10$V1,BDA2019_11$V1,BDA2019_12$V1,BDA2019_13$V1,BDA2019_14$V1,BDA2019_15$V1,
                         BDA2019_16$V1,BDA2019_17$V1,BDA2019_18$V1,BDA2019_19$V1,BDA2019_20$V1)








pMvNorm$varNam <- rep(
  c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
    "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
    "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
  times = nrow(pMvNorm)/126)


rastX <- rasterFromXYZ(XYsedID)

segIDx <- extract(rastX,data2019[,.(XCOORD.YCOORD)])

clims <- weather
mans <- harvscen

for(ij in yearOut){
  for(varX in varRast){
    createTifFromDT(clims, mans, ij, varX, layerDT, startingYear,XYsegID,crsX = crsX)
    print(varNames[varX])
  }
}



if(startingYear==siteTypeX){
  fileDT=paste0("outDT/","init",startingYear,"/","startV_layerall.rdata")
  load(fileDT)
  setkey(XYsegID,segID)
  setkey(startV,segID)
  outXY <- merge(XYsegID,startV,all = T)
  ###remove coordinates NAs
  outXY <- outXY[!is.na(x)]
  
  ###create raster 
  rastX <- rasterFromXYZ(outXY[,c("x","y","value"),with=F])
  crs(rastX) <- crsX
  
  rastName <- paste0("outRast/","init",startingYear,"/","startV_startYear",startingYear,"_layerall.tif")
  writeRaster(rastX,filename = rastName,overwrite=T)
  
}
