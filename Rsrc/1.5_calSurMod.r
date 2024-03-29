setwd("/scratch/project_2000994/PREBASruns/FCMaustria/Rsrc/")
library(MASS)
library(minpack.lm)
library(devtools)
source("localSettings.r") # use settings in local directory if one exists

# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/functions.r")

###check and create output directories
setwd(generalPath)
mkfldr <- paste0("surErrMods/")
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

yearX <- 3

####needs to be changed for ForUnc runs
load(paste0(procDataPath,"/uniqueData.rdata"))
set.seed(123)
data.sample <- uniqueData[sample(1:nrow(uniqueData),nSampleSurrMod)]

sampleID <- 10
rcpfile="CurrClim"

# resample siteType using a uniform distribution 
set.seed(1234)
data.sample$siteType <- sample(1:5,nSampleSurrMod,replace = T)

###load weather inputs
load(climatepath)

# if(rcpfile=="CurrClim"){
#   load(paste(climatepath, rcpfile,".rdata", sep=""))  
#   setnames(dat,"id","climID")
# } else{
#   load(paste(climatepath, rcpfile, sep=""))  
# }

## Prepare the same initial state for all harvest scenarios that are simulated in a loop below

totAreaSample <- sum(data.sample$area)

###check if climID matches
allclIDs <- 1:nrow(PAR)#unique(dat$climID)
samClIds <- unique(data.sample$climID)
newClimIDs <- match(data.sample$climID,samClIds)
data.sample$climID <- newClimIDs
if(!all(samClIds %in% allclIDs)){
  opsClim <- samClIds[which(!samClIds %in% allclIDs)]
  dt = data.table(allclIDs, val = allclIDs) # you'll see why val is needed in a sec
  setnames(dt,c("x","val"))
  # setattr(dt, "sorted", "x")  # let data.table know that w is sorted
  setkey(dt, x) # sorts the data
  # binary search and "roll" to the nearest neighbour
  replX <- dt[J(opsClim), roll = "nearest"]
  data.sample$climID <- mapvalues(data.sample$climID,replX[[1]],replX[[2]])
}


##process weather inputs    
# clim = prep.climate.f(dat, data.sample, startingYear, nYears,startYearWeather)
clim <- list()
clim$PAR <- PAR[samClIds,]
clim$VPD <- VPD[samClIds,]
clim$TAir <- TAir[samClIds,]
clim$Precip <- Precip[samClIds,]
clim$CO2 <- CO2[samClIds,]

# px=pCROB
# px[17,] = pCROB[17,]*0.03
# change nYears
nYears1=year2-startingYear

####
if(cal=="austria"){
  load("/scratch/project_2000994/calibrations/calAustria/outCal/pCROBASaustria.rdata")
  pCROBAS <- pCROBASaustria
}else{
  pCROBAS <- pCROB
}


# Region = nfiareas[ID==r_no, Region]
initPrebas = create_prebas_input.f(clim, pX = pCROBAS,
                                   data.sample=data.sample, nYears = nYears1, 
                                   startSim = startingYear,domSPrun=domSPrun)

###reset names
setnames(data.sample,c("ba","bl","dbh","h","conif"),paste0(c("ba","bl","dbh","h","conif"),startingYear))
print("model initialized")

###run Model
out <- multiPrebas(initPrebas)$multiOut
Vx <- rowSums(out[,yearX,30,,1])
Wabgx <- rowSums(out[,yearX,24,,1]) + rowSums(out[,yearX,31,,1]) + rowSums(out[,yearX,33,,1])
Wblgx <- rowSums(out[,yearX,25,,1]) + rowSums(out[,yearX,32,,1])
Bx <- rowSums(out[,yearX,13,,1])
Bconx <- out[,yearX,13,1,1]
Bblx <- out[,yearX,13,2,1]
Hx <- out[,yearX,11,1,1] * out[,yearX,13,1,1]/Bx +
  out[,yearX,11,2,1] * out[,yearX,13,2,1]/Bx
Dx <- out[,yearX,12,1,1] * out[,yearX,13,1,1]/Bx +
  out[,yearX,12,2,1] * out[,yearX,13,2,1]/Bx
SVI1 <- rowSums(out[,1,43,,1])
Wabg1 <- rowSums(out[,1,24,,1]) + rowSums(out[,1,31,,1]) + rowSums(out[,1,33,,1])
Wblg1 <- rowSums(out[,1,25,,1]) + rowSums(out[,1,32,,1])

PREBx <- data.table(segID=initPrebas$siteInfo[,1],V3 = Vx, B3 = Bx,
                    H3 = Hx, D3 = Dx,Wabg3=Wabgx,Wblg3=Wblgx,
                    Bcon3=Bconx,Bbl3=Bblx,Bh3=Bx*Hx,SVI1=SVI1,
                    Wblg1=Wblg1, Wabg1 = Wabg1)
print("runs completed")


####build surrogate model
### Run settings & functions

# load("C:/Users/minunno/GitHub/satRuns/data/inputUncer.rdata")
load(url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/data/inputUncer.rdata"))

dataX <- data.table(cbind(initPrebas$multiInitVar[,3:5,1],initPrebas$multiInitVar[,5,2],
                          initPrebas$siteInfo[,3],
                          PREBx$V3,PREBx$B3,PREBx$H3,PREBx$D3,
                          PREBx$Bcon3,PREBx$Bbl3,
                          PREBx$Wabg3,PREBx$Wblg3,
                          PREBx$Bh3,PREBx$SVI1,PREBx$Wblg1,PREBx$Wabg1))
setnames(dataX,c("H","D","BAconif","BAbl","st","Vmod","Bmod",
                 "Hmod","Dmod","BAconifmod","BAblmod",
                 "Wabgmod","Wblgmod","Bhmod","SVI1mod",
                 "Wblg1mod","Wabg1mod"))
# if(!all(unique(dataX$st) %in% unique(uniqueData$siteType))) stop("not all siteTypes of the tile are in the sample")

#### Here we use stepwise regression to construct an emulator for stand variables prediction
# dataX$lnVmod<-log(dataX$Vmod)
# dataX$alpha<-NA
dataX$st <- factor(dataX$st)
dataX[,BAtot:=(BAconif+BAbl)]
dataX[,BAh:=BAtot*H]
dataX[,N:=BAtot/(pi*(D/200)^2)]
b = -1.605 ###coefficient of Reineke
dataX[,SDI:=N *(D/10)^b]
dataX[,rootBAconif:=BAconif^0.5]
dataX[,BAconif2:=BAconif^(2)]
full.modelB <-lm(Bmod~H+D+BAh+BAconif+BAbl+st,data=dataX)
step.modelB <- stepAIC(full.modelB, direction = "both",
                       trace = FALSE)
full.modelH <-lm(Hmod~H+D+BAconif+BAbl+st,data=dataX)
step.modelH <- stepAIC(full.modelH, direction = "both",
                       trace = FALSE)
full.modelD <-lm(Dmod~H+D+BAh+BAconif+BAbl+st,data=dataX)
step.modelD <- stepAIC(full.modelD, direction = "both",
                       trace = FALSE)
full.modelBconif <-lm(BAconifmod~H+D+BAh+BAconif+BAbl+st+rootBAconif,data=dataX)
step.modelBconif <- stepAIC(full.modelBconif, direction = "both",
                            trace = FALSE)
full.modelBbl <-lm(BAblmod~H+D+BAh+BAconif+BAbl+st,data=dataX)
step.modelBbl <- stepAIC(full.modelBbl, direction = "both",
                         trace = FALSE)
full.modelV <-lm(Vmod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelV <- stepAIC(full.modelV, direction = "both",
                       trace = FALSE)
full.modelWabg <-lm(Wabgmod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWabg <- stepAIC(full.modelWabg, direction = "both",
                          trace = FALSE)
full.modelWblg <-lm(Wblgmod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWblg <- stepAIC(full.modelWblg, direction = "both",
                          trace = FALSE)
full.modelSVI1 <-lm(SVI1mod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelSVI1 <- stepAIC(full.modelSVI1, direction = "both",
                          trace = FALSE)
full.modelWblg1 <-lm(Wblg1mod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWblg1 <- stepAIC(full.modelWblg1, direction = "both",
                           trace = FALSE)
full.modelWabg1 <-lm(Wabg1mod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWabg1 <- stepAIC(full.modelWabg1, direction = "both",
                           trace = FALSE)
#     
plot(step.modelV$fitted.values,dataX$Vmod,pch=".",col=2)
abline(0,1)
plot(step.modelB$fitted.values,dataX$Bmod,pch=".",col=2)
abline(0,1)
plot(step.modelH$fitted.values,dataX$Hmod,pch=".",col=2)
abline(0,1)
plot(step.modelD$fitted.values,dataX$Dmod,pch=".",col=2)
abline(0,1)
plot(step.modelBconif$fitted.values,dataX$BAconifmod,pch=".",col=2)
abline(0,1)
plot(step.modelBbl$fitted.values,dataX$BAblmod,pch=".",col=2)
abline(0,1)
plot(step.modelWabg$fitted.values,dataX$Wabgmod,pch=".",col=2)
abline(0,1)
plot(step.modelWblg$fitted.values,dataX$Wblgmod,pch=".",col=2)
abline(0,1)
plot(step.modelSVI1$fitted.values,dataX$SVI1mod,pch=".",col=2)
abline(0,1)
plot(step.modelWblg1$fitted.values,dataX$Wblg1mod,pch=".",col=2)
abline(0,1)
plot(step.modelWabg1$fitted.values,dataX$Wabg1mod,pch=".",col=2)
abline(0,1)


# summary(nonlinear)
# summary(step.model)
save(step.modelV,step.modelB,step.modelD,step.modelH,
     step.modelBconif,step.modelBbl,
     step.modelWabg,step.modelWblg,step.modelSVI1,
     step.modelWblg1,step.modelWabg1,
     file=paste0("surErrMods/surMod_Step1_cal",cal,".rdata")) ###needs to be changed update name


#####Run for second step 2018-2021

nYears2 = yearEnd-year2
# Region = nfiareas[ID==r_no, Region]
initPrebas = create_prebas_input.f(clim, pX = pCROBAS,
                                   data.sample, nYears = nYears2, 
                                   startSim = year2,domSPrun=domSPrun)

print("model initialized")


###run Model
out <- multiPrebas(initPrebas)$multiOut
Vx <- rowSums(out[,yearX,30,,1])
Wabgx <- rowSums(out[,yearX,24,,1]) + rowSums(out[,yearX,31,,1]) + rowSums(out[,yearX,33,,1])
Wblgx <- rowSums(out[,yearX,25,,1]) + rowSums(out[,yearX,32,,1])
Bx <- rowSums(out[,yearX,13,,1])
Bconx <- out[,yearX,13,1,1]
Bblx <- out[,yearX,13,2,1]
Hx <- out[,yearX,11,1,1] * out[,yearX,13,1,1]/Bx +
  out[,yearX,11,2,1] * out[,yearX,13,2,1]/Bx
Dx <- out[,yearX,12,1,1] * out[,yearX,13,1,1]/Bx +
  out[,yearX,12,2,1] * out[,yearX,13,2,1]/Bx
SVI1 <- rowSums(out[,1,43,,1])
Wabg1 <- rowSums(out[,1,24,,1]) + rowSums(out[,1,31,,1]) + rowSums(out[,1,33,,1])
Wblg1 <- rowSums(out[,1,25,,1]) + rowSums(out[,1,32,,1])

PREBx <- data.table(segID=initPrebas$siteInfo[,1],V3 = Vx, B3 = Bx,
                    H3 = Hx, D3 = Dx,Wabg3=Wabgx,Wblg3=Wblgx,
                    Bcon3=Bconx,Bbl3=Bblx,Bh3=Bx*Hx,SVI1=SVI1,
                    Wblg1=Wblg1, Wabg1 = Wabg1)

print("runs completed")

####build surrogate model
### Run settings & functions

# load("C:/Users/minunno/GitHub/satRuns/data/inputUncer.rdata")
load(url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/data/inputUncer.rdata"))
dataX <- data.table(cbind(initPrebas$multiInitVar[,3:5,1],initPrebas$multiInitVar[,5,2],
                          initPrebas$siteInfo[,3],
                          PREBx$V3,PREBx$B3,PREBx$H3,PREBx$D3,
                          PREBx$Bcon3,PREBx$Bbl3,
                          PREBx$Wabg3,PREBx$Wblg3,
                          PREBx$Bh3,PREBx$SVI1,PREBx$Wblg1,PREBx$Wabg1))
setnames(dataX,c("H","D","BAconif","BAbl","st","Vmod","Bmod",
                 "Hmod","Dmod","BAconifmod","BAblmod",
                 "Wabgmod","Wblgmod","Bhmod","SVI1mod",
                 "Wblg1mod","Wabg1mod"))

#### Here we use stepwise regression to construct an emulator for stand variables prediction
# dataX$lnVmod<-log(dataX$Vmod)
# dataX$alpha<-NA
dataX$st <- factor(dataX$st)
dataX[,BAtot:=(BAconif+BAbl)]
dataX[,BAh:=BAtot*H]
dataX[,N:=BAtot/(pi*(D/200)^2)]
b = -1.605 ###coefficient of Reineke
dataX[,SDI:=N *(D/10)^b]
dataX[,rootBAconif:=BAconif^0.5]
dataX[,BAconif2:=BAconif^(2)]
full.modelB <-lm(Bmod~H+D+BAh+BAconif+BAbl+st,data=dataX)
step.modelB2 <- stepAIC(full.modelB, direction = "both",
                        trace = FALSE)
full.modelH <-lm(Hmod~H+D+BAconif+BAbl+st,data=dataX)
step.modelH2 <- stepAIC(full.modelH, direction = "both",
                        trace = FALSE)
full.modelD <-lm(Dmod~H+D+BAh+BAconif+BAbl+st,data=dataX)
step.modelD2 <- stepAIC(full.modelD, direction = "both",
                        trace = FALSE)
full.modelBconif <-lm(BAconifmod~H+D+BAh+BAconif+BAbl+st+rootBAconif,data=dataX)
step.modelBconif2 <- stepAIC(full.modelBconif, direction = "both",
                             trace = FALSE)
full.modelBbl <-lm(BAblmod~H+D+BAh+BAconif+BAbl+st,data=dataX)
step.modelBbl2 <- stepAIC(full.modelBbl, direction = "both",
                          trace = FALSE)
full.modelV <-lm(Vmod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelV2 <- stepAIC(full.modelV, direction = "both",
                        trace = FALSE)
full.modelWabg <-lm(Wabgmod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWabg2 <- stepAIC(full.modelWabg, direction = "both",
                           trace = FALSE)
full.modelWblg <-lm(Wblgmod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWblg2 <- stepAIC(full.modelWblg, direction = "both",
                           trace = FALSE)
full.modelSVI1 <-lm(SVI1mod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelSVI1_2 <- stepAIC(full.modelSVI1, direction = "both",
                            trace = FALSE)
full.modelWblg1 <-lm(Wblg1mod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWblg1_2 <- stepAIC(full.modelWblg1, direction = "both",
                             trace = FALSE)
full.modelWabg1 <-lm(Wabg1mod~Hmod+Dmod+Bhmod+BAconifmod+BAblmod+st,data=dataX)
step.modelWabg1_2 <- stepAIC(full.modelWabg1, direction = "both",
                             trace = FALSE)

plot(step.modelV2$fitted.values,dataX$Vmod,pch=".",col=2)
abline(0,1)
plot(step.modelB2$fitted.values,dataX$Bmod,pch=".",col=2)
abline(0,1)
plot(step.modelH2$fitted.values,dataX$Hmod,pch=".",col=2)
abline(0,1)
plot(step.modelD2$fitted.values,dataX$Dmod,pch=".",col=2)
abline(0,1)
plot(step.modelBconif2$fitted.values,dataX$BAconifmod,pch=".",col=2)
abline(0,1)
plot(step.modelBbl2$fitted.values,dataX$BAblmod,pch=".",col=2)
abline(0,1)
plot(step.modelWabg2$fitted.values,dataX$Wabgmod,pch=".",col=2)
abline(0,1)
plot(step.modelWblg2$fitted.values,dataX$Wblgmod,pch=".",col=2)
abline(0,1)
plot(step.modelSVI1_2$fitted.values,dataX$SVI1mod,pch=".",col=2)
abline(0,1)
plot(step.modelWblg1_2$fitted.values,dataX$Wblg1mod,pch=".",col=2)
abline(0,1)
plot(step.modelWabg1_2$fitted.values,dataX$Wabg1mod,pch=".",col=2)
abline(0,1)


# summary(nonlinear)
# summary(step.model)
# summary(nonlinear)
# summary(step.model)
save(step.modelV2,step.modelB2,step.modelD2,step.modelH2,
     step.modelBconif2,step.modelBbl2,
     step.modelWabg2,step.modelWblg2,step.modelSVI1_2,
     step.modelWblg1_2,step.modelWabg1_2,
     file=paste0("surErrMods/surMod_Step2_cal",cal,".rdata")) ###needs to be changed update name
