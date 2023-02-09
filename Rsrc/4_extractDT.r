
library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file in local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")

# Check and create output directories
setwd(generalPath)
mkfldr <- paste0("outDT/","init",startingYear,"/st",siteTypeX)
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

# for(clims in weather){
  clims <- weather
  mans <- harvscen
  # r_no <- 12
  # layerDT <- "all" #1,2,3
  # startingYear <- c(2014,2019)
  for(stYear in startingYear){
    print(stYear)
    createDT(clims,"NoHarv",varDT,layerDT,stYear,extrFun,siteTypeX)
  }
  
  
########Need to clean data  
  # ###cleanData
  # sitesX <- numeric(0)
  # for(stYear in startingYear){
  #   for(vars in c(1:3,7,8,12:15)){
  #     varX <- varNames[varDT[vars]]
  #     fileX=paste0("outDT/",varX,"_",mans,"_",clims,"StartYear",stYear,"Spall.rdata")
  #     load(fileX)
  #     siteNeg <- unique(which(get(varX)<0,arr.ind=T)[,1])
  #     siteNA <- unique(which(is.na(get(varX)),arr.ind=T)[,1])
  #     siteX <- union(siteNA,siteNeg)
  #     sitesX <- union(sitesX,siteX[!is.na(siteX)])
  #   }
  #   
  #   # get(varX)[sitesX[[r_no]]]
  #   sitesX <- as.integer(sitesX)
  #   for(vars in 1:length(varDT)){
  #     varX <- varNames[varDT[vars]]
  #     fileX=paste0("outputDT/",varX,"_",mans,"_",clims,"StartYear",stYear,"Spall.rdata")
  #     load(fileX)
  #     set(get(varX),i=sitesX,j=1:3,value = NA)
  #     save(list=varX,file = fileX)
  #   }
  #   # get(varX)[sitesX[[r_no]]]
  #   # range(get(varX),na.rm=T)
  #   print(paste("cleaned year",stYear))
  # }
  # 
  
