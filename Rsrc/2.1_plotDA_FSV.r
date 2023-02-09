library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")

###check and create output directories
# setwd("~/research/assessCarbon/data/Finland/AC_training_FI_35VLJ")
setwd(generalPath)
mkfldr <- paste0("plots/",paste0("init",startingYear,"/DA",year2))
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

# yearX <- 3
# nSample = 1000 ###number of samples from the error distribution


###join data
# for(i in 1:nSplit){
dataAll <- data.table()
split_id=3
# load("pMvn_FSV_split3.rdata")
load(paste0("posterior/pMvn_FSV_split",split_id,".rdata"))
pMvNorm <- data.table(pMvNorm)[1:126000]
# pMvNorm <- data.table(pMvNorm)
setnames(pMvNorm,"V1","value")
# pMvNormAll <- pMvNorm
# load("pMvn_FSV_split2.rdata")
# pMvNorm <- data.table(pMvNorm)
# pMvNormAll <- rbind(pMvNormAll,pMvNorm)

npix <- nrow(pMvNorm)/126
pMvNorm$run <- rep(c(rep("m2019",42),rep("s2019",42),rep("DA2019",42)),times = npix)

pMvNorm$varNam <- rep(c("H","D","B","pP","pS","pB",paste0("vc_",1:36)),times = npix*3)

vxind <-  c(1,8,15,22,29,36)
covX <- 1:36
covX <- covX[-vxind]
dataAll <- pMvNorm[!varNam %in% paste0("vc_",covX)]
rm(pMvNorm);gc()
library(stringr)
# vrns$varNam <- str_replace(vrns$varNam, "vc", "v")
dataAll$varNam <- str_replace_all(dataAll$varNam, "vc_15", "vB")
dataAll$varNam <- str_replace_all(dataAll$varNam, "vc_1", "vH")
dataAll$varNam <- str_replace_all(dataAll$varNam, "vc_8", "vD")
dataAll$varNam <- str_replace_all(dataAll$varNam, "vc_22", "vpP")
dataAll$varNam <- str_replace_all(dataAll$varNam, "vc_29", "vpS")
dataAll$varNam <- str_replace_all(dataAll$varNam, "vc_36", "vpB")

# 
# dataX <- data.table(dcast(data = pMvNorm,
#                           formula = segID~varNam,value.var = "value"))
# if(i ==  1) dataAll <- dataX
# if(i>1) dataAll <- rbind(dataAll,dataX)
# varAll <- data.table(dcast(data = vrns,
#                            formula = segID~varNam,value.var = "value"))
# rm(vrns); gc()
# print(i) 
# }


####load error models
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/inputUncer.rdata"))
# load("C:/Users/minunno/GitHub/satRuns/data/inputUncer.rdata")
errData$t35VLJ$sigmaFSVda
vs2016 <- t(diag(errData$t35VLJ$sigmaFSVda))
colnames(vs2016) <- c("B","D","H" ,"pP","pS","pB")
vs2016  <- as.data.table(vs2016)

ops <- dataAll[run=="s2019" & varNam %in% c("vB","vH","vD","vpP","vpS","vpB")]
ops[,run:="s2016"]
ops[,value:=NA]
ops[varNam=="vB",value:=vs2016$B]
ops[varNam=="vH",value:=vs2016$H]
ops[varNam=="vD",value:=vs2016$D]
ops[varNam=="vpP",value:=vs2016$pP]
ops[varNam=="vpS",value:=vs2016$pS]
ops[varNam=="vpB",value:=vs2016$pB]

# load("uniqueData3.rdata")
load(paste0("procData/init",startingYear,"/DA",year2,"_split/uniqueData", split_id, ".rdata"))

data2016 <- uniqueDataSplit[segID %in% unique(ops$segID),.(segID,ba,h,dbh,pineP,spruceP,blp)]
rm(uniqueDataSplit); gc()
setnames(data2016,c("segID", "B","H","D" ,"pP","pS","pB"))
setkey(ops,segID)
setkey(data2016,segID)

ops <- rbind(ops,data.table(segID=data2016$segID,value=data2016$B,
                            run="s2016",varNam="B"))
ops <- rbind(ops,data.table(segID=data2016$segID,value=data2016$H,
                            run="s2016",varNam="H"))
ops <- rbind(ops,data.table(segID=data2016$segID,value=data2016$D,
                            run="s2016",varNam="D"))
ops <- rbind(ops,data.table(segID=data2016$segID,value=data2016$pP,
                            run="s2016",varNam="pP"))
ops <- rbind(ops,data.table(segID=data2016$segID,value=data2016$pS,
                            run="s2016",varNam="pS"))
ops <- rbind(ops,data.table(segID=data2016$segID,value=data2016$pB,
                            run="s2016",varNam="pB"))


dataAll <- rbind(dataAll,ops)
rm(ops);gc()

# load("XYsegID.rdata")
load(paste0(procDataPath,"init",startingYear,"/DA",year2,"/XYsegID.rdata"))
setkey(XYsegID,segID)
setkey(dataAll,segID)
dataAll <- merge(XYsegID,dataAll)
dataAll$run <- factor(dataAll$run,levels = c("s2016","m2019","s2019","DA2019"))
save(dataAll,file = paste0("plots/","init",startingYear,"/DA",year2,"/dataForPlots.rdata"))


library(data.table)

library(ggplot2)
library(ggridges)
library(ggpubr)
load("plots/init2016/DA2019/dataForPlots.rdata")

vars <- unique(dataAll$varNam)
pX <- list()
units <- c("H (m)","D (cm)","BA (m2)", "cover (%)", "cover (%)", "cover (%)",
           "variance (m)","variance (cm)","variance (m2)",
           "variance (%)", "variance (%)", "variance (%)")
for(ij in 1:length(vars)){
  varX <- vars[ij]
  pX[[ij]] <- ggplot(dataAll[varNam==varX],
                     aes(x = value, y = run, fill = stat(x))
  ) + 
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "", option = "C") +
    labs(title = varX)+ xlab(units[ij])+ ylab(NULL)
}

# png("plots/init2016/DA2019/fsvPlot.png")
  pFSV <- ggarrange(pX[[1]],pX[[2]],pX[[3]],pX[[4]],pX[[5]],pX[[6]],
          nrow = 2,ncol=3,legend = "none")
  # dev.off()

  # png("plots/init2016/DA2019/varPlot.png")
  pVar = ggarrange(pX[[7]],pX[[8]],pX[[9]],pX[[10]],pX[[11]],pX[[12]],
                   nrow = 2,ncol=3,legend = "none")
  # dev.off()
  p1 <-  pFSV <- ggarrange(pX[[1]],pX[[2]],pX[[3]],pX[[7]],pX[[8]],pX[[9]],
                           nrow = 2,ncol=3,legend = "none")
  p2 <-  pFSV <- ggarrange(pX[[4]],pX[[5]],pX[[6]],pX[[10]],pX[[11]],pX[[12]],
                           nrow = 2,ncol=3,legend = "none")
  ggsave("plots/init2016/DA2019/fsvPlot.pdf",plot = pFSV)
  ggsave("plots/init2016/DA2019/VarfsvPlot.pdf",plot = pVar)
  ggsave("plots/HDBPlot.png",plot = p1)
  ggsave("plots/coverPlot.png",plot = p2)
  # ggsave("plots/init2016/DA2019/fsvPlot.eps",plot = pFSV,device = "eps")
  # ggsave("plots/init2016/DA2019/VarfsvPlot.eps",plot = pVar,device = "eps")
  
