library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings in local directory if one exists


### Check and create output directories
mkfldr <- paste0("surErrMods/")
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}


if(CSCrun){
  load("/scratch/project_2000994/PREBASruns/assessCarbon/data/traningSites.rdata") # in CSC
}else{
  load("data/traningSites.rdata")
}

data2016$year <- 2016
data2019$year <- 2019
dataAll <- rbind(data2016,data2019)

###function to calculate residuals and MVN distribution parameters
calError <- function(dataX){
  resids <- data.table()
  resids$G <- dataX$G.mea- dataX$G.est
  resids$D <- dataX$D.mea- dataX$D.est
  resids$H <- (dataX$H.mea- dataX$H.est)/10
  resids$BAp <- (dataX$PINE.mea- dataX$PINE.est)
  resids$BAsp <- (dataX$SPRUCE.mea- dataX$SPRUCE.est)
  resids$BAb <- (dataX$BL.mea- dataX$BL.est)
  sigmaFSVda <- cov(resids)
  muFSVda <- colMeans(resids)
  ###include V in the multivariate normal distribution for ite type DA
  resids2 <- resids[,.(G,D,H)]
  resids2$V <- (dataX$V.mea- dataX$V.est)
  muSTda <- colMeans(resids2)
  sigmaSTda <- cov(resids2)
  corMatFSVda <- cor(resids)
  corMatSTda <- cor(resids2)
  return(list(residsFSVda = resids,muFSVda=muFSVda,sigmaFSVda=sigmaFSVda,
              corMatFSVda=corMatFSVda,
              residsSTda = resids,muSTda=muSTda,sigmaSTda=sigmaSTda,
              corMatSTda=corMatSTda))
}
errData <- list()
errData$all <- calError(dataAll)
errData$y2016$all <- calError(data2016)
errData$y2019$all <- calError(data2019)
errData$y2016$t35VLJ <- calError(data2016[S2Tile == "35VLJ"])
errData$y2019$t35VLJ <- calError(data2019[S2Tile == "35VLJ"])
errData$y2016$t35VNL <- calError(data2016[S2Tile == "35VNL"])
errData$y2019$t35VNL <- calError(data2019[S2Tile == "35VNL"])
errData$y2016$t34VEQ <- calError(data2016[S2Tile == "34VEQ"])
errData$y2019$t34VEQ <- calError(data2019[S2Tile == "34VEQ"])
errData$y2016$t35WMN <- calError(data2016[S2Tile == "35WMN"])
errData$y2019$t35WMN <- calError(data2019[S2Tile == "35VLJ"])
# 
errData$t35VLJ <- calError(dataAll[S2Tile == "35VLJ"])
errData$t35VNL <- calError(dataAll[S2Tile == "35VNL"])
errData$t34VEQ <- calError(dataAll[S2Tile == "34VEQ"])
errData$t35WMN <- calError(dataAll[S2Tile == "35WMN"])

if(CSCrun){
  save(errData,file="/scratch/project_2000994/PREBASruns/assessCarbon/data/inputUncer.rdata") # in CSC
}else{
  save(errData,file="data/inputUncer.rdata")
}



#### Probit regression between sitetype.ref and sitetype.est 
summary(dataAll)
dataAll[Class.est>5,Class.est:=5]
dataAll[Class.ref>5,Class.ref:=5]
dataAll$alpha<-NA
dataAll$alpha[dataAll$Class.est==1]<- 0.2833333 #mean(pCROB['alfar1',1:3])
dataAll$alpha[dataAll$Class.est==2]<- 0.4066667 #mean(pCROB['alfar2',1:3])
dataAll$alpha[dataAll$Class.est==3]<- 0.4966667 #mean(pCROB['alfar3',1:3])
dataAll$alpha[dataAll$Class.est==4]<- 0.6233333 #mean(pCROB['alfar4',1:3])
dataAll$alpha[dataAll$Class.est==5]<- 0.7866667 #mean(pCROB['alfar5',1:3])

dataAll$Class.est.f<-factor(dataAll$Class.est,levels = 1:5)
dataAll$Class.ref.f<-factor(dataAll$Class.ref,levels = 1:5)


calProbit <- function(dataX, yearX="all", tileX="all"){
  if(yearX=="all" & tileX=="all"){
    dataX <- dataX[,.(Class.ref.f, Class.est,H.est,D.est,G.est,PINE.est, SPRUCE.est, BL.est)]
  }else if(yearX=="all" & tileX!="all"){
    dataX <- dataX[S2Tile==tileX,.(Class.ref.f, Class.est,H.est,D.est,G.est,PINE.est, SPRUCE.est, BL.est)]
  }else if(yearX!="all" & tileX=="all"){
    dataX <- dataX[year==yearX,.(Class.ref.f, Class.est,H.est,D.est,G.est,PINE.est, SPRUCE.est, BL.est)]
  }else{
    dataX <- dataX[year==yearX & S2Tile==tileX,.(Class.ref.f, Class.est,H.est,D.est,G.est,PINE.est, SPRUCE.est, BL.est)]
  }
  setnames(dataX,c("st.mea.f","st", "H", "D", "BAtot","BApPer","BAspPer","BAbPer"))
  full.probit<-polr(st.mea.f ~ st+H+D+BAtot+BApPer+BAspPer+BAbPer,data=dataX ,method='probit')
  step.probit <- stepAIC(full.probit, direction = "both",
                         trace = FALSE)
  return(step.probit)
}

step.probit <- list()
step.probit$all <- calProbit(dataAll)
step.probit$y2016$all <- calProbit(dataAll,yearX=2016)
step.probit$y2019$all <- calProbit(dataAll,yearX=2019)
step.probit$y2016$t35VLJ <- calProbit(dataAll,yearX=2016,tileX="35VLJ")
step.probit$y2019$t35VLJ <- calProbit(dataAll,yearX=2019,tileX="35VLJ")
# step.probit$y2016$t35VNL <- calProbit(dataAll,yearX=2016,tileX="35VNL")
step.probit$y2019$t35VNL <- calProbit(dataAll,yearX=2019,tileX="35VNL")
step.probit$y2016$t34VEQ <- calProbit(dataAll,yearX=2016,tileX="34VEQ")
step.probit$y2019$t34VEQ <- calProbit(dataAll,yearX=2019,tileX="34VEQ")
step.probit$y2019$t35WMN <- calProbit(dataAll,yearX="all",tileX="35WMN")

# save(step.probit,file = paste0(generalPath,'surErrMods/stProbit.rdata'))
if(CSCrun){
  save(step.probit,file="/scratch/project_2000994/PREBASruns/assessCarbon/data/step.probit.rdata") # in CSC
}else{
  save(step.probit,file="data/step.probit.rdata")
}

#### logistic function to predict the pure forests

dataAll$max.pro.ref<-apply(dataAll[, c('PINE.mea','SPRUCE.mea','BL.mea')], 1, max)
dataAll$max.pro.est<-apply(dataAll[, c('PINE.est','SPRUCE.est','BL.est')], 1, max)
dataAll$pure.ref<-F
dataAll$pure.ref[dataAll$max.pro.ref>=100]<-T
# dataAll$pure.ref<-as.factor(dataAll$pure.ref)

logistic.model <- list()

logistic.model$all <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                    data = dataAll)

yearX=2016; tileX = "35VLJ"
dataX <- dataAll[year==yearX & S2Tile == tileX]
logistic.model[[paste0("y",yearX)]][[paste0("t",tileX)]] <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                                   data = dataX)
yearX=2019; tileX = "35VLJ"
dataX <- dataAll[year==yearX & S2Tile == tileX]
logistic.model[[paste0("y",yearX)]][[paste0("t",tileX)]] <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                                                                data = dataX)
yearX=2019; tileX = "35VNL"
dataX <- dataAll[year==yearX & S2Tile == tileX]
logistic.model[[paste0("y",yearX)]][[paste0("t",tileX)]] <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                                                                data = dataX)
yearX=2016; tileX = "34VEQ"
dataX <- dataAll[year==yearX & S2Tile == tileX]
logistic.model[[paste0("y",yearX)]][[paste0("t",tileX)]] <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                                                                data = dataX)
yearX=2019; tileX = "34VEQ"
dataX <- dataAll[year==yearX & S2Tile == tileX]
logistic.model[[paste0("y",yearX)]][[paste0("t",tileX)]] <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                                                                data = dataX)
yearX=2016; tileX = "35WMN"
dataX <- dataAll[year==yearX & S2Tile == tileX]
logistic.model[[paste0("y",yearX)]][[paste0("t",tileX)]] <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                                                                data = dataX)
yearX=2019; tileX = "35WMN"
dataX <- dataAll[year==yearX & S2Tile == tileX]
logistic.model[[paste0("y",yearX)]][[paste0("t",tileX)]] <- glm(formula = pure.ref ~ max.pro.est, family = binomial(link = "logit"), 
                                                                data = dataX)
# plot(dataAll$max.pro.est,predict(logistic.model,type="response"))
logisticPureF <- logistic.model
if(CSCrun){
  save(logisticPureF,file="/scratch/project_2000994/PREBASruns/assessCarbon/data/logisticPureF.rdata") # in CSC
}else{
  save(logisticPureF,file="data/logisticPureF.rdata")
}

