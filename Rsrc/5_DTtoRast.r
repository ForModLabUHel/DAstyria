library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/functions.r")

# Check and create output directories
setwd(generalPath)
mkfldr <- "outRast/"
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}


load(paste0(procDataPath,"/XYsegID.rdata"))  
crsX <- crs(raster(baRast2015))
# rastX <- rasterFromXYZ(XYsegID)
# # segIDx <- extract(rastX,data2019[S2Tile==tileX,.(XCOORD,YCOORD)])
# 
# plot(rastX)
# points(data2019[S2Tile==tileX,.(XCOORD,YCOORD)],pch=20,col=2)



varNam <- c("H","D","B","conif","brleaf")
dataSource <- c("PREBAS","S2","DA")
years <- c(2018,2021)

outNames <- c(paste0(c("H","D","B","conif","brleaf"),"_m_2018"),
              paste0(c("H","D","B","conif","brleaf"),"_s2_2018"),
              paste0(c("H","D","B","conif","brleaf"),"_da_2018"),
              paste0(c("H","D","B","conif","brleaf"),"_m_2021"),
              paste0(c("H","D","B","conif","brleaf"),"_s2_2021"),
              paste0(c("H","D","B","conif","brleaf"),"_da_2021"))

namesX <- c("segID","Hmod","Dmod","Bmod","Conifmod","Blmod")
###combine data
allData <- data.table()
# load("posterior/pMvn_FSV_split1.rdata")
# 
# setkey(pMvNorm,segID)
# pMvNorm$outNam <- rep(outNames, times = length(unique(pMvNorm$segID)))

# oo <- merge(XYsegID, pMvNorm[outNam=="H_da_2021"], by.x = "segID", 
#             by.y = "segID", all.x = TRUE, all.y = FALSE)

# setkey(XYsegID,segID)
# Dda2019 <- Dm2019 <- Ds2019 <- 
#   Hda2019 <- Hm2019 <- Hs2019 <- 
#   Bda2019 <- Bm2019 <- Bs2019 <- data.table()
for(i in 1:nSplit){
  load(paste0("posterior/pMvn_FSV_split",i,".rdata"))
  setkey(pMvNorm,segID)
  pMvNorm$varID <- rep(1:30,nrow(pMvNorm)/30)
  
  procData <- reshape(pMvNorm,                          # Reshape data from long to wide format
                      timevar   = "varID", 
                      idvar     = "segID", 
                      direction = "wide",
  )
  DA1 <- procData[,c(1,12:16)]
  DA2 <- procData[,c(1,27:31)]
  setnames(DA1,namesX)
  setnames(DA2,namesX)
  DA1$year <- 2018
  DA2$year <- 2021
  procData <- rbind(DA1,DA2)
  allData <- rbind(allData,procData)
  # Dda2019 <- rbind(Dda2019,merge(XYsegID,pMvNorm[varNam=="DDA2019"]))
  # Dm2019 <- rbind(Dm2019,merge(XYsegID,pMvNorm[varNam=="Dm2019"]))
  # Ds2019 <- rbind(Ds2019,merge(XYsegID,pMvNorm[varNam=="Ds2019"]))
  # Hda2019 <- rbind(Hda2019,merge(XYsegID,pMvNorm[varNam=="HDA2019"]))
  # Hm2019 <- rbind(Hm2019,merge(XYsegID,pMvNorm[varNam=="Hm2019"]))
  # Hs2019 <- rbind(Hs2019,merge(XYsegID,pMvNorm[varNam=="Hs2019"]))
  # Bda2019 <- rbind(Bda2019,merge(XYsegID,pMvNorm[varNam=="BDA2019"]))
  # Bm2019 <- rbind(Bm2019,merge(XYsegID,pMvNorm[varNam=="Bm2019"]))
  # Bs2019 <- rbind(Bs2019,merge(XYsegID,pMvNorm[varNam=="Bs2019"]))
  print(i)
}


####calculate V and biomasses
allData[,Bhmod:=Hmod*Bmod]
allData[,BAconifmod:=Bmod*Conifmod/(Conifmod+Blmod)]
allData[,BAblmod:= Bmod*Blmod/(Conifmod+Blmod)]
allData[,st:=as.factor(3)]

allData$V <- -999
allData$Wabg <- -999
allData$Wblg <- -999
setkey(allData,segID)
allData[year==2018]$V <- pmax(0.,predict(step.modelV,newdata=allData[year==2018]))
allData[year==2018]$Wabg <- pmax(0.,predict(step.modelWabg,newdata=allData[year==2018]))
allData[year==2018]$Wblg <- pmax(0.,predict(step.modelWblg,newdata=allData[year==2018]))
allData[year==2021]$V <- pmax(0.,predict(step.modelV2,newdata=allData[year==2021]))
allData[year==2021]$Wabg <- pmax(0.,predict(step.modelWabg2,newdata=allData[year==2021]))
allData[year==2021]$Wblg <- pmax(0.,predict(step.modelWblg2,newdata=allData[year==2021]))


oo <- merge(XYsegID, allData[year==2018], by.x = "segID",
            by.y = "segID", all.x = TRUE, all.y = FALSE)

H_rast <- rasterFromXYZ(oo[,.(x,y,Hmod)])
D_rast <- rasterFromXYZ(oo[,.(x,y,Dmod)])
G_rast <- rasterFromXYZ(oo[,.(x,y,Bmod)])
conif_rast <- rasterFromXYZ(oo[,.(x,y,Conifmod)])
brleaf_rast <- rasterFromXYZ(oo[,.(x,y,Blmod)])
V_rast <- rasterFromXYZ(oo[,.(x,y,V)])
Wabg_rast <- rasterFromXYZ(oo[,.(x,y,Wabg)])
Wblg_rast <- rasterFromXYZ(oo[,.(x,y,Wblg)])

writeRaster(H_rast,file="outRast/H_da2018.tif")
writeRaster(D_rast,file="outRast/D_da2018.tif")
writeRaster(G_rast,file="outRast/G_da2018.tif")
writeRaster(conif_rast,file="outRast/conif_da2018.tif")
writeRaster(brleaf_rast,file="outRast/brleaf_da2018.tif")
writeRaster(V_rast,file="outRast/V_da2018.tif")
writeRaster(Wabg_rast,file="outRast/Wabg_da2018.tif")
writeRaster(Wblg_rast,file="outRast/Wblg_da2018.tif")



###rasters 2021
oo <- merge(XYsegID, allData[year==2021], by.x = "segID",
            by.y = "segID", all.x = TRUE, all.y = FALSE)

H_rast <- rasterFromXYZ(oo[,.(x,y,Hmod)])
D_rast <- rasterFromXYZ(oo[,.(x,y,Dmod)])
G_rast <- rasterFromXYZ(oo[,.(x,y,Bmod)])
conif_rast <- rasterFromXYZ(oo[,.(x,y,Conifmod)])
brleaf_rast <- rasterFromXYZ(oo[,.(x,y,Blmod)])
V_rast <- rasterFromXYZ(oo[,.(x,y,V)])
Wabg_rast <- rasterFromXYZ(oo[,.(x,y,Wabg)])
Wblg_rast <- rasterFromXYZ(oo[,.(x,y,Wblg)])

writeRaster(H_rast,file="outRast/H_da2021.tif")
writeRaster(D_rast,file="outRast/D_da2021.tif")
writeRaster(G_rast,file="outRast/G_da2021.tif")
writeRaster(conif_rast,file="outRast/conif_da2021.tif")
writeRaster(brleaf_rast,file="outRast/brleaf_da2021.tif")
writeRaster(V_rast,file="outRast/V_da2021.tif")
writeRaster(Wabg_rast,file="outRast/Wabg_da2021.tif")
writeRaster(Wblg_rast,file="outRast/Wblg_da2021.tif")





# dataX <- as.data.frame(allData[segID==1]$V1)
# rownames(dataX) <- outNames
# dataX <- data.table(t(dataX))
# dataSel <- data.table(Hmod=dataX$H_da_2018,
#                       Dmod=dataX$D_da_2018,
#                       Bhmod=dataX$B_da_2018 * dataX$H_da_2018,
#                       BAconifmod=dataX$B_da_2018 * dataX$conif_da_2018/(dataX$conif_da_2018 + dataX$brleaf_da_2018),
#                       BAblmod=dataX$B_da_2018 * dataX$brleaf_da_2018/(dataX$conif_da_2018 + dataX$brleaf_da_2018),
#                       st=as.factor(3))
# max(0.,predict(step.modelV,newdata=dataSel))
# max(0.,predict(step.modelWabg,newdata=dataSel))
# max(0.,predict(step.modelWblg,newdata=dataSel))
# 
# # for(i in 2:nSplit){
# #   load(paste0("posterior/pMvn_FSV_split",i,".rdata"))
# #   setkey(oo,segID)
# #   setkey(pMvNorm,segID)
# #   
# #   pMvNorm$varNam <- rep(
# #     c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
# #       "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
# #       "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
# #     times = nrow(pMvNorm)/126)
# #   
# #   assign(paste0("Hm2019_",i), merge(data2019, pMvNorm[varNam=="Hm2019",1:2], by.x = "segID", 
# #                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("HDA2019_",i), merge(data2019, pMvNorm[varNam=="HDA2019",1:2], by.x = "segID", 
# #                                      by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("Hs2019_",i), merge(data2019, pMvNorm[varNam=="Hs2019",1:2], by.x = "segID", 
# #                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("Bm2019_",i), merge(data2019, pMvNorm[varNam=="Bm2019",1:2], by.x = "segID", 
# #                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("BDA2019_",i), merge(data2019, pMvNorm[varNam=="BDA2019",1:2], by.x = "segID", 
# #                                      by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("Bs2019_",i), merge(data2019, pMvNorm[varNam=="Bs2019",1:2], by.x = "segID", 
# #                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("Dm2019_",i), merge(data2019, pMvNorm[varNam=="Dm2019",1:2], by.x = "segID", 
# #                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("DDA2019_",i), merge(data2019, pMvNorm[varNam=="DDA2019",1:2], by.x = "segID", 
# #                                      by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   assign(paste0("Ds2019_",i), merge(data2019, pMvNorm[varNam=="Ds2019",1:2], by.x = "segID", 
# #                                     by.y = "segID", all.x = FALSE, all.y = FALSE))
# #   print(i)
# # }
# # 
# 
# 
# 
# data2019res <- rbind(Ds2019_10,Ds2019_11)
# data2019res <- rbind(data2019res,Ds2019_12)
# data2019res <- rbind(data2019res,Ds2019_13)
# data2019res <- rbind(data2019res,Ds2019_14)
# data2019res <- rbind(data2019res,Ds2019_15)
# data2019res <- rbind(data2019res,Ds2019_16)
# data2019res <- rbind(data2019res,Ds2019_17)
# data2019res <- rbind(data2019res,Ds2019_18)
# data2019res <- rbind(data2019res,Ds2019_19)
# data2019res <- rbind(data2019res,Ds2019_20)
# 
# setnames(data2019res,"V1","Ds2019")
# 
# data2019res$DDA2019 <- c(DDA2019_10$V1,DDA2019_11$V1,DDA2019_12$V1,DDA2019_13$V1,DDA2019_14$V1,DDA2019_15$V1,
#                          DDA2019_16$V1,DDA2019_17$V1,DDA2019_18$V1,DDA2019_19$V1,DDA2019_20$V1)
# 
# data2019res$Dm2019 <- c(Dm2019_10$V1,Dm2019_11$V1,Dm2019_12$V1,Dm2019_13$V1,Dm2019_14$V1,Dm2019_15$V1,
#                         Dm2019_16$V1,Dm2019_17$V1,Dm2019_18$V1,Dm2019_19$V1,Dm2019_20$V1)
# 
# 
# 
# data2019res$Hm2019 <- c(Hm2019_10$V1,Hm2019_11$V1,Hm2019_12$V1,Hm2019_13$V1,Hm2019_14$V1,Hm2019_15$V1,
#                         Hm2019_16$V1,Hm2019_17$V1,Hm2019_18$V1,Hm2019_19$V1,Hm2019_20$V1)
# 
# data2019res$Hs2019 <- c(Hs2019_10$V1,Hs2019_11$V1,Hs2019_12$V1,Hs2019_13$V1,Hs2019_14$V1,Hs2019_15$V1,
#                         Hs2019_16$V1,Hs2019_17$V1,Hs2019_18$V1,Hs2019_19$V1,Hs2019_20$V1)
# 
# data2019res$HDA2019 <- c(HDA2019_10$V1,HDA2019_11$V1,HDA2019_12$V1,HDA2019_13$V1,HDA2019_14$V1,HDA2019_15$V1,
#                          HDA2019_16$V1,HDA2019_17$V1,HDA2019_18$V1,HDA2019_19$V1,HDA2019_20$V1)
# 
# 
# data2019res$Bm2019 <- c(Bm2019_10$V1,Bm2019_11$V1,Bm2019_12$V1,Bm2019_13$V1,Bm2019_14$V1,Bm2019_15$V1,
#                         Bm2019_16$V1,Bm2019_17$V1,Bm2019_18$V1,Bm2019_19$V1,Bm2019_20$V1)
# 
# data2019res$Bs2019 <- c(Bs2019_10$V1,Bs2019_11$V1,Bs2019_12$V1,Bs2019_13$V1,Bs2019_14$V1,Bs2019_15$V1,
#                         Bs2019_16$V1,Bs2019_17$V1,Bs2019_18$V1,Bs2019_19$V1,Bs2019_20$V1)
# 
# data2019res$BDA2019 <- c(BDA2019_10$V1,BDA2019_11$V1,BDA2019_12$V1,BDA2019_13$V1,BDA2019_14$V1,BDA2019_15$V1,
#                          BDA2019_16$V1,BDA2019_17$V1,BDA2019_18$V1,BDA2019_19$V1,BDA2019_20$V1)
# 
# 
# 
# 
# 
# 
# 
# 
# pMvNorm$varNam <- rep(
#   c("Hm2019","Dm2019","Bm2019","perPm2019","perSPm2019","perBm2019",rep("varcov1",36),
#     "Hs2019","Ds2019","Bs2019","perPs2019","perSPs2019","perBs2019",rep("varcov2",36),
#     "HDA2019","DDA2019","BDA2019","perPDA2019","perSPDA2019","perBDA2019",rep("varcov3",36)),
#   times = nrow(pMvNorm)/126)
# 
# 
# rastX <- rasterFromXYZ(XYsedID)
# 
# segIDx <- extract(rastX,data2019[,.(XCOORD.YCOORD)])
# 
# clims <- weather
# mans <- harvscen
# 
# for(ij in yearOut){
#   for(varX in varRast){
#     createTifFromDT(clims, mans, ij, varX, layerDT, startingYear,XYsegID,crsX = crsX)
#     print(varNames[varX])
#   }
# }
# 
# 
# 
# if(startingYear==siteTypeX){
#   fileDT=paste0("outDT/","init",startingYear,"/","startV_layerall.rdata")
#   load(fileDT)
#   setkey(XYsegID,segID)
#   setkey(startV,segID)
#   outXY <- merge(XYsegID,startV,all = T)
#   ###remove coordinates NAs
#   outXY <- outXY[!is.na(x)]
#   
#   ###create raster 
#   rastX <- rasterFromXYZ(outXY[,c("x","y","value"),with=F])
#   crs(rastX) <- crsX
#   
#   rastName <- paste0("outRast/","init",startingYear,"/","startV_startYear",startingYear,"_layerall.tif")
#   writeRaster(rastX,filename = rastName,overwrite=T)
#   
# }
