# Run settings 
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists


if(startingYear != siteTypeX){
  siteTypeX = startingYear
  warning("siteTypeX changed to startingYear")
} 

setwd(generalPath)

if (splitRun) {
  # If output is set to be split to smaller parts (splitRun = TRUE), create separate
  # folder for the split data tables.
  mkfldr_split <- paste0("procData/",paste0("init",startingYear,"/ForUn",year2,"_split"))
  if(!dir.exists(file.path(generalPath, mkfldr_split))) {
    dir.create(file.path(generalPath, mkfldr_split), recursive = TRUE)
  }
}

mkfldr <- paste0("procData/",paste0("init",startingYear,"/ForUn",year2))
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}


###extract CurrClim IDs
rastRef <- raster(baRast2015)
if(testRun){
  extNew <- extent(rastRef)
  extNew[1:2]   <- (extent(rastRef)[1] + (extent(rastRef)[2]))/2 + c(-10000,10000)
  extNew[3:4]   <- (extent(rastRef)[3] + (extent(rastRef)[4]))/2 + c(-10000,10000)
  rastRef <- crop(rastRef,extNew)
  maxSitesRun <- maxSitesRunTest
}


climID <- raster(climIDpath)
climID <- crop(climID,rastRef)
climID <- raster(vals=values(climID),ext=extent(rastRef),crs=crs(rastRef),
                 nrows=dim(rastRef)[1],ncols=dim(rastRef)[2])

fileNames <- c(baRast2015,baRast2018,baRast2021,
               blPerRast2015,blPerRast2018,blPerRast2021,
               dbhRast2015,dbhRast2018,dbhRast2021,
               vRast2015,vRast2018,vRast2021,
               hRast2015,hRast2018,hRast2021,
               conifPerRast2015,conifPerRast2018,conifPerRast2021,
               # sprucePerRast,
               siteTypeRast2015,siteTypeRast2018,siteTypeRast2021,
               if (mgmtmask==T) mgmtmaskRast)

years <- c(2015,2018,2021)
varNames <- c(paste0("ba",years),paste0("bl",years),
              paste0("dbh",years),paste0("v",years),paste0("h",years),
              paste0("conif",years),paste0("siteType",years),
              if (mgmtmask==T) "mgmtmask","climID")

for(i in 1:length(fileNames)){
  rastX <- raster(fileNames[i])
  if(testRun) rastX <- crop(rastX,rastRef)    ####if it is a test run rasters are cropped to a smaller area
  rastX <- raster(vals=values(rastX),ext=extent(rastRef),crs=crs(rastRef),
                  nrows=dim(rastRef)[1],ncols=dim(rastRef)[2])
  
  dataX <- data.table(rasterToPoints(rastX))
  setnames(dataX,c("x","y",varNames[i]))
  setkey(dataX,x,y)
  if(i==1){
    data.all <- dataX 
    setkey(data.all,x,y)
  }else{
    data.all <- merge(data.all,dataX,all=T)
    setkey(data.all,x,y)
  }
  print(fileNames[i])
}

###attach weather ID
data.all$climID <- extract(climID,data.all[,.(x,y)])

##filter data 
data.all <- data.all[!ba %in% baNA]
data.all <- data.all[!blp %in% blPerNA]
data.all <- data.all[!dbh %in% dbhNA]
data.all <- data.all[!v %in% vNA]
data.all <- data.all[!h %in% hNA]
data.all <- data.all[!pineP %in% pinePerNA]
data.all <- data.all[!spruceP %in% sprucePerNA]
data.all <- data.all[!siteType %in% siteTypeNA]

####convert data to prebas units
data.all <- data.all[, ba := ba * baConv]
data.all <- data.all[, blp := blp * blPerConv]
data.all <- data.all[, dbh := dbh * dbhConv]
data.all <- data.all[, v := v * vConv]
data.all <- data.all[, h := h * hConv]
data.all <- data.all[, pineP := pineP * pinePerConv]
data.all <- data.all[, spruceP := spruceP * sprucePerConv]
data.all <- data.all[, siteType := siteType * siteTypeConv]

####group pixels by same values
data.all[, segID := .GRP, by = .(ba, blp,dbh, h, pineP, spruceP, siteType, climID)]
data.all[,clCut:=0]

#####I'm excluding from the runs the areas that have been clearcutted and have ba=0 
# data.all[h==0. & dbh==0 & ba==0,clCut:=1]
data.all[ba==0,clCut:=1]

###calculate tree density
data.all[clCut==0,N:=ba/(pi*(dbh/200)^2)]

####check where H is below minimum initial height and replace
smallH <- intersect(which(data.all$h < initH), which(data.all$clCut==0))
data.all[smallH, h:=initH]

###check where density is too high and replase stand variables with initial conditions
tooDens <- intersect(which(data.all$N> maxDens), which(data.all$clCut==0))
data.all[tooDens,h:=initH]
data.all[tooDens,ba:=initBA]
data.all[tooDens,dbh:=initDBH]
data.all[tooDens,N:=initN]


data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType ==1, blp:=1  ]
data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType <= 3 & siteType > 1, spruceP:=1  ]
data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType >= 4, pineP:=1  ]

####Count segID pix
data.all[, npix:=.N, segID]

# uniqueData <- data.table()
####find unique initial conditions
uniqueData <- unique(data.all[clCut2015==0,.(segID,npix,ba2015,blPerRast2015,dbh2015,h2015,conif2015,siteType2015,N2015,climID)])

uniqueData[,uniqueKey:=1:nrow(uniqueData)]
setkey(uniqueData, uniqueKey)
# range(uniqueData$N)
# uniqueData[,area:=npix*resX^2/10000]

###assign ID to similar pixels
XYsegID <- data.all[,.(x,y,segID)]


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
uniqueData[,sampleRow:=1:length(h2015),by=sampleID]

segID <- numeric(0)
for(i in 1:nSamples){
  sampleX <- samples[[i]]
  segID <- c(segID,sampleX$segID)
}


save(data.all,file=paste0(procDataPath,"init",startingYear,"/","ForUn",yearEnd,"/allData.rdata"))         ### All data
save(uniqueData,file=paste0(procDataPath,"init",startingYear,"/","ForUn",yearEnd,"/uniqueData.rdata"))    ### unique pixel combination to run in PREBAS
save(samples,file=paste0(procDataPath,"init",startingYear,"/","ForUn",yearEnd,"/samples.rdata"))    ### unique pixel combination to run in PREBAS
save(XYsegID,segID,file=paste0(procDataPath,"init",startingYear,"/","ForUn",yearEnd,"/XYsegID.rdata"))    ### Coordinates and segID of all pixels


#### If needed (splitRun = TRUE), unique data is split to separate tables here to enable 
#    running further scripts in multiple sections. Number of split parts is defined in splitRange variable (in settings).
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
    save(uniqueDataSplit,file=paste0(procDataPath,"init",startingYear,"/","ForUn",yearEnd,"_split/uniqueData",i,".rdata"))  
    
    rm(uniqueDataSplit)
  }
}



