library(MASS)
library(devtools)
library(mvtnorm)

source("localSettings.r")
# Run settings 
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/functions.r")

###check and create output directories
setwd(generalPath)
mkfldr <- paste0("posterior/")
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

yearX <- 3
nSample = 1000 ###number of samples from the error distribution

# Load unique data.
# If data is processed in split parts, define to variable split_id which split part to process (in batch job script).
# If splitRun is not needed, the unique data dataset for the whole tile is loaded.
if (splitRun) {
  uniqueData_file <- load(paste0("procData/split/uniqueData", split_id, ".rdata"))
  uniqueData <- get(uniqueData_file)
  rm(list = uniqueData_file)
  rm(uniqueData_file)
  gc()
  # load(paste0("procData/init",startingYear,"/calST_split/stProbMod",split_id,".rdata"))
  # stProb <- data.table(stProb)
} else{
  load("procData/uniqueData.rdata")
  # load("stProbMod.rdata")
  # stProb <- data.table(stProb)
}

####load error models
load(url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/data/inputUncer.rdata"))
# load(url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/data/logisticPureF.rdata"))
# load(url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/data/step.probit.rdata"))
###load surrMods !!!change name
load("surErrMods/surMod_Step1.rdata")
load("surErrMods/surMod_Step2.rdata")

uniqueData[,BAconif2015:= (ba2015 * conif2015/(conif2015+bl2015))]
uniqueData[,BAbl2015:= (ba2015 * bl2015/(conif2015+bl2015))]
uniqueData[,BAconif2018 := (ba2018 * conif2018/(conif2018+bl2018))]
uniqueData[,BAbl2018:= (ba2018 * bl2018/(conif2018+bl2018))]
uniqueData[,BAconif2021 := (ba2021 * conif2021/(conif2021+bl2021))]
uniqueData[,BAbl2021:= (ba2021 * bl2021/(conif2021+bl2021))]

dataSurMod <- uniqueData[,.(segID,h2015,dbh2015,BAconif2015,BAbl2015,siteType2015,
                            siteType2018,v2018,ba2018,h2018,dbh2018,
                            BAconif2018,BAbl2018,
                            v2021,ba2021,h2021,dbh2021,
                            BAconif2021,BAbl2021,siteType2021)] 
setnames(dataSurMod,c("segID","H","D","BAconif","BAbl","st1",
                      "st2","V2","ba2","h2","dbh2",
                      "BAconif2","BAbl2", "v3","ba3",
                      "h3","dbh3","BAconif3","BAbl3","st3"))


dataSurMod[,BAconifPer:=.(BAconif/sum(BAconif,BAbl)*100),by=segID]
dataSurMod[,BAblPer:=.(BAbl/sum(BAconif,BAbl)*100),by=segID]
dataSurMod[,BAtot:=.(sum(BAconif,BAbl)),by=segID]
dataSurMod[,BAconifPer2:=.(BAconif2/sum(BAconif2,BAbl2)*100),by=segID]
dataSurMod[,BAblPer2:=.(BAbl2/sum(BAconif2,BAbl2)*100),by=segID]
dataSurMod[,BAtot2:=.(sum(BAconif2,BAbl2)),by=segID]
dataSurMod[,BAconifPer3:=.(BAconif3/sum(BAconif3,BAbl3)*100),by=segID]
dataSurMod[,BAblPer3:=.(BAbl3/sum(BAconif3,BAbl3)*100),by=segID]
dataSurMod[,BAtot3:=.(sum(BAconif3,BAbl3)),by=segID]

nSeg <- nrow(dataSurMod)  ##200

# dataSurMod <- merge(dataSurMod,stProb)

#pMvNormBASE <- matrix(NA,127,nSeg)
#pMvNorm <- matrix(NA,127,nSeg)
pMvNorm <- data.table()
errDataX <- errData$y2015$all
postErr <- errDataX$muFSVda

print("running DA")

if(parallelRun){
  
  system.time({ # PARALLEL PROCESSING
    # Number of cores used for processing is defined with mc.cores argument (in settings). mc.cores = 1 disables parallel processing.
    #pMvNorm <- mclapply(1:ncol(pMvNorm), function(i){
    pMvNorm <- mclapply(1, function(i){
      # pMvNorm[,i] <- pSVDA(dataSurMod[i],nSample,year1=startingYear,
      #                     year2=year2,tileX=tileX)
      pMvNorm <- dataSurMod[1:10000, pSVDA_2steps(.SD,
                                                  nSample = nSample,
                                                  errData1 = errDataX,
                                                  errData2 = errDataX,
                                                  errData3 = errDataX,
                                                  errorPost = postErr,
                                                  step.modelH=step.modelH,
                                                  step.modelD=step.modelD,
                                                  step.modelB=step.modelB,
                                                  step.modelBconif=step.modelBconif,
                                                  step.modelBbl=step.modelBbl,
                                                  step.modelH2=step.modelH2,
                                                  step.modelD2=step.modelD2,
                                                  step.modelB2=step.modelB2,
                                                  step.modelBconif2=step.modelBconif2,
                                                  step.modelBbl2=step.modelBbl2,
                                                  dist="mvnorm"),
                            by = segID]
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
  pMvNorm <- data.table()
  system.time({ # SERIAL PROCESSING
    for(i in 1:10000){ #nSeg
      pMvNorm <- rbind(pMvNorm,
                       data.table(segID=i,
                                  V1=pSVDA_2steps(segIDx=dataSurMod[i],nSample,
                                                  errDataX,errDataX,errDataX,
                                                  postErr,
                                                  step.modelH=step.modelH,step.modelD=step.modelD,
                                                  step.modelB=step.modelB,
                                                  step.modelBconif=step.modelBconif,
                                                  step.modelBbl=step.modelBbl,
                                                  step.modelH2=step.modelH2,step.modelD2=step.modelD2,
                                                  step.modelB2=step.modelB2,
                                                  step.modelBconif2=step.modelBconif2,
                                                  step.modelBbl2=step.modelBbl2,
                                                  dist="mvnorm")[2:61]))
    }
  })
}

if(splitRun) {  ##  If run in split parts, output produced with each split part is saved temporarily (as pMvn_FSV_split*split_id*.rdata).
  ##  Split outputs will be combined and saved to pMvn_FSV.rdata file when the last split part is processed
  
  save(pMvNorm, file = paste0("posterior/pMvn_FSV_split",split_id,".rdata"))
  print(paste("ID",split_id,"DA completed"))
  
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
  print("DA completed")
}
