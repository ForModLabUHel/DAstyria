
# Run settings 
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")
setwd(generalPath)
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists

# source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/rmvweisd.r")
# source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/Fweibull_Arithmetic Mean and Variance.R")

# if(modifiedSettings) {
#   source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings.r") # in CSC
# }

# Create folders for outputs.
if (splitRun) {
  # If output is set to be split to smaller parts (splitRun = TRUE), create separate
  # folder for the split data tables.
  mkfldr_split <- paste0("procData/",paste0("init",startingYear,"/DA",year2,"_split"))
  if(!dir.exists(file.path(generalPath, mkfldr_split))) {
    dir.create(file.path(generalPath, mkfldr_split), recursive = TRUE)
  }
} 
  
mkfldr <- paste0("procData/",paste0("init",startingYear,"/DA",year2))
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}



###extract CurrClim IDs
rastX <- raster(baRast2015)
if(testRun){
  extNew <- extent(rastX)
  extNew[2]   <- (extent(rastX)[1] + (extent(rastX)[2] - extent(rastX)[1])*fracTest)
  extNew[4]   <- (extent(rastX)[3] + (extent(rastX)[4] - extent(rastX)[3])*fracTest)
  rastX <- crop(rastX,extNew)
  maxSitesRun <- maxSitesRunTest
}

climID <- raster(climIDpath)

rm(rastX)
gc()


fileNames <- c(baRast2015,baRast2018,baRast2021,
               blPerRast2015,blPerRast2018,blPerRast2021,
               dbhRast2015,dbhRast2018,dbhRast2021,
               vRast2015,vRast2018,vRast2021,
               hRast2015,hRast2018,hRast2021,
               conifPerRast2015,conifPerRast2018,conifPerRast2021,
               # sprucePerRast,
               siteTypeRast2015,siteTypeRast2018,siteTypeRast2021,
               if (mgmtmask==T) mgmtmaskRast)


for(i in 1:length(fileNames)){
  rastX <- raster(fileNames[i])
  if(testRun) rastX <- crop(rastX,extNew)    ####if it is a test run rasters are cropped to a smaller area
  dataX <- data.table(rasterToPoints(rastX))
  if(i==1){
    data.all <- dataX 
  }else{
    data.all <- merge(data.all,dataX)
  }
  print(fileNames[i])
}


###attach weather ID
data.all$climID <- extract(climID,data.all[,.(x,y)])
# dataX <- data.table(rasterToPoints(climIDs))
# data.all <- merge(data.all,dataX)
years <- c(2015,2018,2021)
setnames(data.all,c("x","y",paste0("ba",years),paste0("bl",years),
                    paste0("dbh",years),paste0("v",years),paste0("h",years),
                    paste0("conif",years),paste0("siteType",years),
                    if (mgmtmask==T) "mgmtmask","climID"))

##filter data 
if(mgmtmask==T) data.all <- data.all[mgmtmask == 0]
data.all <- data.all[!ba2015 %in% baNA]
data.all <- data.all[!ba2018 %in% baNA]
data.all <- data.all[!ba2021 %in% baNA]
data.all <- data.all[!bl2015 %in% blPerNA]
data.all <- data.all[!bl2018 %in% blPerNA]
data.all <- data.all[!bl2021 %in% blPerNA]
data.all <- data.all[!dbh2015 %in% dbhNA]
data.all <- data.all[!dbh2018 %in% dbhNA]
data.all <- data.all[!dbh2021 %in% dbhNA]
data.all <- data.all[!v2015 %in% vNA]
data.all <- data.all[!v2018 %in% vNA]
data.all <- data.all[!v2021 %in% vNA]
data.all <- data.all[!h2015 %in% hNA]
data.all <- data.all[!h2018 %in% hNA]
data.all <- data.all[!h2021 %in% hNA]
data.all <- data.all[!conif2015 %in% pinePerNA]
data.all <- data.all[!conif2018 %in% pinePerNA]
data.all <- data.all[!conif2021 %in% pinePerNA]
# data.all <- data.all[!spruceP %in% sprucePerNA]
data.all <- data.all[!siteType2015 %in% siteTypeNA]
data.all <- data.all[!siteType2018 %in% siteTypeNA]
data.all <- data.all[!siteType2021 %in% siteTypeNA]

####convert data to prebas units
data.all <- data.all[, ba2015 := ba2015 * baConv]
data.all <- data.all[, ba2018 := ba2018 * baConv]
data.all <- data.all[, ba2021 := ba2021 * baConv]
data.all <- data.all[, bl2015 := bl2015 * blPerConv]
data.all <- data.all[, bl2018 := bl2018 * blPerConv]
data.all <- data.all[, bl2021 := bl2021 * blPerConv]
data.all <- data.all[, dbh2015 := dbh2015 * dbhConv]
data.all <- data.all[, dbh2018 := dbh2018 * dbhConv]
data.all <- data.all[, dbh2021 := dbh2021 * dbhConv]
data.all <- data.all[, v2015 := v2015 * vConv]
data.all <- data.all[, v2018 := v2018 * vConv]
data.all <- data.all[, v2021 := v2021 * vConv]
data.all <- data.all[, h2015 := h2015 * hConv]
data.all <- data.all[, h2018 := h2018 * hConv]
data.all <- data.all[, h2021 := h2021 * hConv]
data.all <- data.all[, conif2015 := conif2015 * pinePerConv]
data.all <- data.all[, conif2018 := conif2018 * pinePerConv]
data.all <- data.all[, conif2021 := conif2021 * pinePerConv]
data.all <- data.all[, siteType2015 := siteType2015 * siteTypeConv]
data.all <- data.all[, siteType2018 := siteType2018 * siteTypeConv]
data.all <- data.all[, siteType2021 := siteType2021 * siteTypeConv]


###tocheck!
if(siteTypeX==year2){
  data.all[,siteType:=siteType2]  
}else if(siteTypeX==startingYear){
  data.all[,siteType:=siteType1]  
}else{
  data.all[,siteType:=siteTypeX]  
}
data.all[siteType2015>5,siteType2015:=5]
data.all[siteType2018>5,siteType2018:=5]
data.all[siteType2021>5,siteType2021:=5]


#####I'm excluding from the runs the areas that have been clearcutted and have ba=0 
# data.all[h==0. & dbh==0 & ba==0,clCut:=1]
data.all[,clCut2015:=0];data.all[,clCut2018:=0];data.all[,clCut2021:=0]
data.all[ba2015==0,clCut2015:=1]
data.all[ba2018==0,clCut2018:=1]
data.all[ba2021==0,clCut2021:=1]

###calculate tree density
data.all[clCut2015==0,N2015:=ba2015/(pi*(dbh2015/200)^2)]
data.all[clCut2018==0,N2018:=ba2018/(pi*(dbh2018/200)^2)]
data.all[clCut2021==0,N2021:=ba2021/(pi*(dbh2021/200)^2)]

!!!!
####check where H is below minimum initial height and replace
smallH <- intersect(which(data.all$h2015 < initH), which(data.all$clCut2015==0))
data.all[smallH, h2015:=initH]
smallH <- intersect(which(data.all$h2018 < initH), which(data.all$clCut2018==0))
data.all[smallH, h2018:=initH]
smallH <- intersect(which(data.all$h2021 < initH), which(data.all$clCut2021==0))
data.all[smallH, h2021:=initH]

###check where density is too high and replace stand variables with initial conditions
tooDens <- intersect(which(data.all$N> maxDens), which(data.all$clCut==0))
data.all[tooDens,h:=initH]
data.all[tooDens,ba:=initBA]
data.all[tooDens,dbh:=initDBH]
data.all[tooDens,N:=initN]


data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType ==1, blp:=1  ]
data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType <= 3 & siteType > 1, spruceP:=1  ]
data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType >= 4, pineP:=1  ]

###!!!!!!!!!!!!########careful with this part##########!!!!!!!!#########

####calculate dV, dBA, dH, dDBH
# data.all[,dV := v2-v]
data.all[,dVy := (v2-v)/(year2 - startingYear)]
# data.all[,dBA := ba2-ba]
data.all[,dBAy := (ba2-ba)/(year2 - startingYear)]
# data.all[,dH := h2-h]
data.all[,dHy := (h2-h)/(year2 - startingYear)]
# data.all[,dDBH := dbh2-dbh]
data.all[,dDBHy := (dbh2-dbh)/(year2 - startingYear)]

####group pixels by same values
data.all[, segID := .GRP, by = .(ba, blp,dbh, h, pineP, spruceP, 
                                 siteType1,siteType2, climID,dVy,v2,
                                 dBAy,ba2,dHy,h2,dDBHy,dbh2,
                                 pineP2, spruceP2,blp2)]

####Count segID pix
data.all[, npix:=.N, segID]

# uniqueData <- data.table()
####find unique initial conditions
uniqueData <- unique(data.all[clCut==0,.(segID,npix,climID,ba,blp,dbh,h,pineP,spruceP,
                                         siteType1,siteType2,dBAy,ba2,dVy,v2,
                                         dHy,h2,dDBHy,dbh2,pineP2, spruceP2,blp2)])

uniqueData[,uniqueKey:=1:nrow(uniqueData)]
setkey(uniqueData, uniqueKey)
# uniqueData[,N:=ba/(pi*(dbh/200)^2)]
# range(uniqueData$N)

uniqueData[,area:=npix*resX^2/10000]

###assign ID to similar pixels
XYsegID <- data.all[,.(x,y,segID)]

###!!!!!!!!!!!!########end careful with this part##########!!!!!!!!#########

# nSamples <- ceiling(dim(uniqueData)[1]/20000)
# sampleID <- 1
# 
# for(sampleID in sampleIDs){
#   set.seed(1)
#   samplesX <- split(uniqueData, sample(1:nSample, nrow(uniqueData), replace=T))
#   sampleX <- ops[[sampleID]]
#   sampleX[,area := N*resX^2/10000]
#   # sampleX[,id:=climID]
# }


nSamples <- ceiling(dim(uniqueData)[1]/maxSitesRun)
set.seed(1)
sampleset <- sample(1:nSamples, nrow(uniqueData),  replace=T)
samples <- split(uniqueData, sampleset) 

# adding sampleID, sampleRow (= row within sample) 
uniqueData[,sampleID:=sampleset]
uniqueData[,sampleRow:=1:length(h),by=sampleID]

segID <- numeric(0)
for(i in 1:nSamples){
  sampleX <- samples[[i]]
  segID <- c(segID,sampleX$segID)
}

save(data.all,file=paste0(procDataPath,"init",startingYear,"/DA",year2,"/allData.rdata"))         ### All data
save(uniqueData,file=paste0(procDataPath,"init",startingYear,"/DA",year2,"/uniqueData.rdata"))    ### unique pixel combination to run in PREBAS
save(samples,file=paste0(procDataPath,"init",startingYear,"/DA",year2,"/samples.rdata"))    ### unique pixel combination to run in PREBAS
save(XYsegID,segID,file=paste0(procDataPath,"init",startingYear,"/DA",year2,"/XYsegID.rdata"))    ### Coordinates and segID of all pixels

#### If needed (splitRun = TRUE), unique data is split to separate tables here to enable 
#    running further scripts in multiple sections. Number of split parts is defined in splitRange variable (in settings).
#    Running in multiple sections reduces processing time

if (splitRun) {
  
  # Create split_id column which is used in splitting the table. NOTICE that the last section might be of unequal size compared to the others.
  split_length <- ceiling(nrow(uniqueData)/length(splitRange))
  uniqueData <- uniqueData[, split_id := NA]


  uniqueData$split_id[1:split_length] <- 1
  for (i in 2:(max(splitRange)-1)) {
    uniqueData$split_id[((i-1)*split_length+1):(split_length*i)] <- i
  }
  uniqueData$split_id[((length(splitRange)-1)*split_length+1):(nrow(uniqueData))] <- length(splitRange)

  # Split the table to list of elements. Splitting is done based on the split_id.
  split_list <- split(uniqueData,uniqueData$split_id)

  for (i in 1:max(splitRange)) {
  
    # Convert the split results to separate data tables
    uniqueDataSplit <- as.data.table(split_list[[i]])
  
    # Remove split_id column
    uniqueDataSplit <- uniqueDataSplit[, split_id:=NULL]
  
    # Save split tables
    save(uniqueDataSplit,file=paste0(procDataPath,"init",startingYear,"/DA",year2,"_split/uniqueData",i,".rdata"))  
  
    rm(uniqueDataSplit)
  }
}
