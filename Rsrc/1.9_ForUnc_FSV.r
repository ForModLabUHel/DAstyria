library(MASS)
library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file in local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")

###check and create output directories
###check and create output directories
setwd(generalPath)
mkfldr <- "posterior/"
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

yearX <- 3
nSample = 1000 ###number of samples from the error distribution

# Load unique data.
# If data is processed in split parts, define to variable split_id which split part to process (in batch job script).
# If splitRun is not needed, the unique data dataset for the whole tile is loaded.
if (splitRun) {
  uniqueData_file <- load(paste0("procData/init",startingYear,"/ForUn",yearEnd,"_split/uniqueData", split_id, ".rdata"))
  uniqueData <- get(uniqueData_file)
  rm(list = uniqueData_file)
  rm(uniqueData_file)
} else{
  load(paste0("procData/init",startingYear,"/ForUn",yearEnd,"/uniqueData.rdata"))  
}

####load error models
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/inputUncer.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/logisticPureF.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/step.probit.rdata"))
###load surrMods
###this needs to be changed loading the 1.5 surr mod calibrated with the right uniqueData
load("surErrMods/surMod.rdata")

uniqueData[,BAp:= (ba * pineP/(pineP+spruceP+blp))]
uniqueData[,BAsp:= (ba * spruceP/(pineP+spruceP+blp))]
uniqueData[,BAb:= (ba * blp/(pineP+spruceP+blp))]
# uniqueData[,BAp2:= (ba2 * pineP2/(pineP2+spruceP2+blp2))]
# uniqueData[,BAsp2:= (ba2 * spruceP2/(pineP2+spruceP2+blp2))]
# uniqueData[,BAb2:= (ba2 * blp2/(pineP2+spruceP2+blp2))]

dataSurMod <- uniqueData[,.(h,dbh,BAp,BAsp,BAb,siteType,segID)] 
setnames(dataSurMod,c("H","D","BAp","BAsp","BAb","st","segID"))


dataSurMod[,BApPer:=.(BAp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAspPer:=.(BAsp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAbPer:=.(BAb/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAtot:=.(sum(BAp,BAsp,BAb)),by=segID]

# nSeg <- nrow(dataSurMod)  ##200
# load("stProbMod.rdata")
# colnames(stProb) <- paste0("pST",1:5)
# # dataSurMod <- cbind(dataSurMod[1:nSeg],stProb[1:nSeg,])
# pMvn <- matrix(NA,nSeg,43)
# system.time({
#   for(i in 1:nSeg){
#     pMvn[i,] <- prForUnc(dataSurMod[i],nSample = nSample,yearUnc=startingYear,
#                          tileX=tileX)
#     if (i %% 100 == 0) { print(i) }
#   }
# })
# nSeg=100
# test <- dataSurMod[1:nSeg]

pMvn <- data.table()

if(parallelRun){
  system.time({ # PARALLEL PROCESSING
    # Number of cores used for processing is defined with mc.cores argument (in settings). mc.cores = 1 disables parallel processing.
    pMvn <- mclapply(1, function(i){
      pMvn <- dataSurMod[, prForUnc(.SD,nSample = nSample,yearUnc=startingYear,
                                    tileX=tileX), by =segID]
    },mc.cores = coresN)
  })
 
  pMvn <- as.data.table(pMvn) # convert from list to data table
  
} else {
  system.time({ # SERIAL PROCESSING
    pMvn <- dataSurMod[, prForUnc(.SD,nSample = nSample,yearUnc=startingYear,
                                  tileX=tileX), by =segID]
  })
}

if(splitRun) { 
  save(pMvn, file = paste0("posterior/pMvn_ForUnc_split",split_id,".rdata")) # save processed split data, combine split parts after processing all parts
  
} else {
  save(pMvn,file="posterior/pMvn_ForUnc.rdata") # pMvn finished for the whole dataset
}
