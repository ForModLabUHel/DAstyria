library(MASS)

# Run settings 
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists

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
  uniqueData_file <- load(paste0("procData/init",startingYear,"/DA",year2,"_split/uniqueData", split_id, ".rdata"))
  uniqueData <- get(uniqueData_file)
  rm(list = uniqueData_file)
  rm(uniqueData_file)
  gc()
  load(paste0("procData/init",startingYear,"/calST_split/stProbMod",split_id,".rdata"))
  stProb <- data.table(stProb)
} else{
  load(paste0("procData/init",startingYear,"/DA",year2,"/uniqueData.rdata"))  
  load("stProbMod.rdata")
  stProb <- data.table(stProb)
}

####load error models
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/inputUncer.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/logisticPureF.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/step.probit.rdata"))
###load surrMods !!!change name
load("surErrMods/surMod.rdata")

uniqueData[,BAp:= (ba * pineP/(pineP+spruceP+blp))]
uniqueData[,BAsp:= (ba * spruceP/(pineP+spruceP+blp))]
uniqueData[,BAb:= (ba * blp/(pineP+spruceP+blp))]
uniqueData[,BAp2:= (ba2 * pineP2/(pineP2+spruceP2+blp2))]
uniqueData[,BAsp2:= (ba2 * spruceP2/(pineP2+spruceP2+blp2))]
uniqueData[,BAb2:= (ba2 * blp2/(pineP2+spruceP2+blp2))]

dataSurMod <- uniqueData[,.(segID,h,dbh,BAp,BAsp,BAb,siteType1,
                            siteType2,v2,ba2,h2,dbh2,
                            BAp2,BAsp2,BAb2)] 
setnames(dataSurMod,c("segID","H","D","BAp","BAsp","BAb","st1",
                      "st2","V2","ba2","h2","dbh2",
                      "BAp2","BAsp2","BAb2"))


dataSurMod[,BApPer:=.(BAp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAspPer:=.(BAsp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAbPer:=.(BAb/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAtot:=.(sum(BAp,BAsp,BAb)),by=segID]
dataSurMod[,BApPer2:=.(BAp2/sum(BAp2,BAsp2,BAb2)*100),by=segID]
dataSurMod[,BAspPer2:=.(BAsp2/sum(BAp2,BAsp2,BAb2)*100),by=segID]
dataSurMod[,BAbPer2:=.(BAb2/sum(BAp2,BAsp2,BAb2)*100),by=segID]
dataSurMod[,BAtot2:=.(sum(BAp2,BAsp2,BAb2)),by=segID]

nSeg <- nrow(dataSurMod)  ##200

dataSurMod <- merge(dataSurMod,stProb)

#pMvNormBASE <- matrix(NA,127,nSeg)
#pMvNorm <- matrix(NA,127,nSeg)
pMvNorm <- data.table()

if(parallelRun){

  system.time({ # PARALLEL PROCESSING
    # Number of cores used for processing is defined with mc.cores argument (in settings). mc.cores = 1 disables parallel processing.
    #pMvNorm <- mclapply(1:ncol(pMvNorm), function(i){
    pMvNorm <- mclapply(1, function(i){
      # pMvNorm[,i] <- pSVDA(dataSurMod[i],nSample,year1=startingYear,
      #                     year2=year2,tileX=tileX)
      pMvNorm <- dataSurMod[, pSVDA(.SD,nSample = nSample,year1=startingYear,
                                    year2=year2,tileX=tileX), by = segID]
    
      },mc.cores = coresN)
  })

  # Modify the output to correct form 
  #pMvNormDF <- as.data.frame(pMvNorm)
  #colnames(pMvNormDF) <- colnames(pMvNormBASE)
  #pMvNorm <- data.matrix(pMvNormDF)
  pMvNorm <- as.data.table(pMvNorm)

} else {

  # system.time({ # SERIAL PROCESSING
    # for(i in 1:nSeg){
    #   pMvNorm <- pSVDA(dataSurMod[i],nSample,year1=startingYear,
    #                        year2=year2,tileX=tileX)
    #   if (i %% 100 == 0) { print(i) }
    # }
  # })

  system.time({
   pMvNorm <- dataSurMod[, pSVDA(.SD,nSample = nSample,year1=startingYear,
                                year2=year2,tileX=tileX), by = segID]
  })
}

if(splitRun) {  ##  If run in split parts, output produced with each split part is saved temporarily (as pMvn_FSV_split*split_id*.rdata).
                ##  Split outputs will be combined and saved to pMvn_FSV.rdata file when the last split part is processed
  
  save(pMvNorm, file = paste0("posterior/pMvn_FSV_split",split_id,".rdata"))
  
# if(split_id==max(splitRange)){
#   
#   List <- list()
#   
#   # Iterate through all split parts of pMvn. 
#   # Bind split parts to a single matrix
#   for (i in 1:max(splitRange)) {
#     pMvn_file <- load(paste0("pMvn_FSV_split",i,".rdata"))
#     pMvn_split <- get(pMvn_file)
#     rm(pMvn_file)
#     
#     List[[i]] <- pMvn_split
#     
#     rm(pMvn_split)
#   }
#   
#   pMvNorm <- do.call(cbind,List)
#   save(pMvNorm,file="pMvn_FSV.rdata") # pMvNorm finished for the whole dataset
#   
# } 
  
} else {
save(pMvNorm,file="posterior/pMvn_FSV.rdata") # pMvNorm finished for the whole dataset
}
