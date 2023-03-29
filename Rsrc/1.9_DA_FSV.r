library(MASS)

# Run settings 
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/functions.r")

###check and create output directories
###check and create output directories
setwd(generalPath)
mkfldr <- paste0("posterior",year2,"/")
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

  system.time({ # SERIAL PROCESSING
   for(i in 1:1000){
  
   errDataX <- errData$y2015$all
   postErr <- errData$muFSVda
   # errData2 <- errData$y2015$all
   nX <- i
   
   kk <- pSVDA_2steps(segIDx=dataSurMod[nX],nSample,
                      errDataX,errDataX,errDataX,
                     postErr,
                     step.modelH=step.modelH,step.modelD=step.modelD,
                     step.modelB=step.modelB,
                     step.modelBconif=step.modelBconif,
                     step.modelBbl=step.modelBbl)
   }
  })
   xx <- dataSurMod[nX]
   xx$H <- kk[61]
   xx$D <- kk[62]
   xx$BAtot <- kk[63]
   xx$BAconif <- kk[63]*kk[64]/100
   xx$BAbl <- kk[63]*kk[65]/100
   xx$st2 <- xx$st3
   xx$ba2 <- xx$ba3
   xx$h2 <- xx$h3
   xx$dbh2 <- xx$dbh3
   xx$BAconif2 <- xx$BAconif3
   xx$BAbl2 <- xx$BAbl3
   xx$BAconifPer2 <- xx$BAconifPer3
   xx$BAblPer2 <- xx$BAblPer3
   xx$BAtot2 <- xx$BAtot3
   errDataX <- list(); errDataX$muFSVda <- errData2$muFSVda
   errDataX$sigmaFSVda <- matrix(kk[66:90],5,5)
   
   kk2 <- pSVDA(xx,nSample,
               errDataX,errData2,
               step.modelHx=step.modelH,step.modelDx=step.modelD,
               step.modelBx=step.modelB,
               step.modelBconifx=step.modelBconif,
               step.modelBblx=step.modelBbl)
   
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
