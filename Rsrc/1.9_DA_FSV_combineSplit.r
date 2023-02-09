# Run settings 
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists


###set working directory
setwd(generalPath)


# Create an empty data table where the combined data is added
pMvNormBind <- data.table()

# Iterate through the split pMvNorm data tables and bind split parts to a single data table 
for (i in splitRange) {
  pMvNorm_file <- load(paste0("pMvn_FSV_split",i,".rdata"))
  pMvNorm_split <- get(pMvNorm_file)
  
  rm(list = pMvNorm_file)
  rm(pMvNorm_file)
  
  pMvNormBind <- rbindlist(list(pMvNormBind, pMvNorm_split))
  
  rm(pMvNorm_split)
}

pMvNorm <- pMvNormBind
save(pMvNorm,file="pMvn_FSV.rdata") # pMvNorm finished for the whole dataset
