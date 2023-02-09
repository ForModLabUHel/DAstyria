library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings in local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")


###check and create output directories
setwd(generalPath)
mkfldr <- paste0("surErrMods/")
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

yearX <- 3

####needs to be changed for ForUnc runs
load(paste0(procDataPath,"init",startingYear,"/DA",year2,"/samples.rdata"))  
  
  sampleID <- 10
  rcpfile="CurrClim"

  # resample siteType using a uniform distribution 
  samples[[sampleID]]$siteType <- sample(1:5,nrow(samples[[sampleID]]),replace = T)

  if(rcpfile=="CurrClim"){
    load(paste(climatepath, rcpfile,".rdata", sep=""))  
    setnames(dat,"id","climID")
  } else{
    load(paste(climatepath, rcpfile, sep=""))  
  }
  
  ## Prepare the same initial state for all harvest scenarios that are simulated in a loop below
    data.sample <- samples[[sampleID]]
    totAreaSample <- sum(data.sample$area)
    
    ###check if climID matches
    allclIDs <- unique(dat$climID)
    samClIds <- unique(data.sample$climID)
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
    
    clim = prep.climate.f(dat, data.sample, startingYear, nYears,startYearWeather)
    
    # Region = nfiareas[ID==r_no, Region]
    initPrebas = create_prebas_input.f(r_no, clim, data.sample, nYears = nYears,
                                       startingYear = startingYear,domSPrun=domSPrun)
    
    print("model initialized")

    
    ###run Model
    out <- multiPrebas(initPrebas)$multiOut
    Vx <- rowSums(out[,yearX,30,,1])
    Bx <- rowSums(out[,yearX,13,,1])
    Bpx <- out[,yearX,13,1,1]
    Bspx <- out[,yearX,13,2,1]
    Bdx <- out[,yearX,13,3,1]
    Hx <- out[,yearX,11,1,1] * out[,yearX,13,1,1]/Bx +
      out[,yearX,11,2,1] * out[,yearX,13,2,1]/Bx +
      out[,yearX,11,3,1] * out[,yearX,13,3,1]/Bx
    Dx <- out[,yearX,12,1,1] * out[,yearX,13,1,1]/Bx +
      out[,yearX,12,2,1] * out[,yearX,13,2,1]/Bx +
      out[,yearX,12,3,1] * out[,yearX,13,3,1]/Bx
    
    PREBx <- data.table(segID=initPrebas$siteInfo[,1],V3 = Vx, B3 = Bx,
                        H3 = Hx, D3 = Dx,
                        Bp3=Bpx,Bsp3=Bspx,Bd3=Bdx)
    print("runs completed")
    
    
    
    ####build surrogate model
    library(MASS)
    library(minpack.lm)
    ### Run settings & functions

    # load("C:/Users/minunno/GitHub/satRuns/data/inputUncer.rdata")
    load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/inputUncer.rdata"))
    # load(paste0(procDataPath,"init",startingYear,"/","st",siteTypeX,"/XYsegID.rdata"))  
    # load(paste0("output/init",startingYear,"/","st",siteTypeX,"/CurrClim_sample1.rdata"))  
    # load(paste0("procData/init",startingYear,"/","st",siteTypeX,"/uniqueData.rdata"))  
    # load(paste0("outDT/init",startingYear,"/","st",siteTypeX,"/V_NoHarv_CurrClimlayerall.rdata"))  
    # Vmod2019 <- rowSums(out[,yearX,6,])
    # load(paste0("initPrebas/init",startingYear,"/","st",siteTypeX,"/CurrClim_sample1.rdata"))  
    dataX <- data.table(cbind(initPrebas$multiInitVar[,3:5,1],initPrebas$multiInitVar[,5,2],
                              initPrebas$multiInitVar[,5,3],initPrebas$siteInfo[,3],
                              PREBx$V3,PREBx$B3,PREBx$H3,PREBx$D3,
                              PREBx$Bp3,PREBx$Bsp3,PREBx$Bd3))
    setnames(dataX,c("H","D","BAp","BAsp","BAb","st","Vmod","Bmod",
                     "Hmod","Dmod","BApmod","BAspmod","BAdmod"))
    # if(!all(unique(dataX$st) %in% unique(uniqueData$siteType))) stop("not all siteTypes of the tile are in the sample")
    
    #### Here we use stepwise regression to construct an emulator for stand variables prediction
    # dataX$lnVmod<-log(dataX$Vmod)
    # dataX$alpha<-NA
    dataX$st <- factor(dataX$st)
    dataX[,BAtot:=(BAp+BAsp+BAb)]
    dataX[,BAh:=BAtot*H]
    dataX[,N:=BAtot/(pi*(D/200)^2)]
    b = -1.605 ###coefficient of Reineke
    dataX[,SDI:=N *(D/10)^b]
    dataX[,rootBAp:=BAp^0.5]
    dataX[,BAp2:=BAp^(2)]
    full.modelV <-lm(Vmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
    step.modelV <- stepAIC(full.modelV, direction = "both",
                          trace = FALSE)
    full.modelB <-lm(Bmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
    step.modelB <- stepAIC(full.modelB, direction = "both",
                           trace = FALSE)
    full.modelH <-lm(Hmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
    step.modelH <- stepAIC(full.modelH, direction = "both",
                           trace = FALSE)
    full.modelD <-lm(Dmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
    step.modelD <- stepAIC(full.modelD, direction = "both",
                           trace = FALSE)
    full.modelBp <-lm(BApmod~H+D+SDI+BAh+BAp+BAsp+BAb+st+rootBAp,data=dataX)
    step.modelBp <- stepAIC(full.modelBp, direction = "both",
                           trace = FALSE)
    full.modelBsp <-lm(BAspmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
    step.modelBsp <- stepAIC(full.modelBsp, direction = "both",
                             trace = FALSE)
    full.modelBd <-lm(BAdmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
    step.modelBd <- stepAIC(full.modelBd, direction = "both",
                             trace = FALSE)
    # start<-as.vector(full.model$coefficients)
### Anonther option is to use nonlinear regression, which differed in error assumption. 
    # full.model0<-lm(lnVmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
    # start<-as.vector(full.model0$coefficients)
    # nonlinear<-nlsLM(Vmod~exp(a+b*H+c*D+d*SDI+e*BAh+f*BAp+g*BAsp+h*BAb+
    #                           i2*as.numeric(st==2)+
    #                           i3*as.numeric(st==3)+
    #                           i4*as.numeric(st==4)+
    #                           i5*as.numeric(st==5)
    #                         ),
    #                data=dataX,start = list(a=start[1],
    #                                  b=start[2],
    #                                  c=start[3],
    #                                  d=start[4],
    #                                  e=start[5],
    #                                  f=start[6],
    #                                  g=start[7],
    #                                  h=start[8],
    #                                  i2=start[9],
    #                                  i3=start[10],
    #                                  i4=start[11],
    #                                  i5=start[12]
    #                                  ))
    #     
    plot(step.modelV$fitted.values,dataX$Vmod,pch=".")
    abline(0,1)
    plot(step.modelB$fitted.values,dataX$Bmod,pch=".")
    abline(0,1)
    plot(step.modelH$fitted.values,dataX$Hmod,pch=".")
    abline(0,1)
    plot(step.modelD$fitted.values,dataX$Dmod,pch=".")
    abline(0,1)
    plot(step.modelBp$fitted.values,dataX$BApmod,pch=".")
    abline(0,1)
    plot(step.modelBsp$fitted.values,dataX$BAspmod,pch=".")
    abline(0,1)
    plot(step.modelBd$fitted.values,dataX$BAdmod,pch=".")
    abline(0,1)
    
    # summary(nonlinear)
    # summary(step.model)
    save(step.modelV,step.modelB,step.modelD,step.modelH,
         step.modelBp,step.modelBsp,step.modelBd,
         file="surErrMods/surMod.rdata") ###needs to be changed update name
    
