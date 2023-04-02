
####function to initialize the model
create_prebas_input.f = function(clim, data.sample, nYears, pX=pCROB,
                                 startSim=0,domSPrun=0,nLayers=2) { # dat = climscendataset
  #domSPrun=0 initialize model for mixed forests according to data inputs 
  #domSPrun=1 initialize model only for dominant species 
  nSites <- nrow(data.sample)
  
  ###site Info matrix. nrow = nSites, cols: 1 = siteID; 2 = climID; 3=site type;
  ###4 = nLayers; 5 = nSpecies;
  ###6=SWinit;   7 = CWinit; 8 = SOGinit; 9 = Sinit
  
  siteInfo <- matrix(c(NA,NA,NA,160,0,0,20,nLayers,nLayers,413,0.45,0.118),nSites,12,byrow = T)
  #siteInfo <- matrix(c(NA,NA,NA,3,3,160,0,0,20),nSites,9,byrow = T)
  siteInfo[,1] <- data.sample$segID
  siteInfo[,2] <- data.sample$climID
  siteInfo[,3] <- data.sample$siteType
  
  # litterSize <- matrix(0,3,3)
  # litterSize[1,1:2] <- 30
  # litterSize[1,3] <- 10
  # litterSize[2,] <- 2
  
  ###Initialise model
  # initVardension nSites,variables, nLayers
  # variables: 1 = species; 2 = Age; 3 = H; 4=dbh; 5 = ba; 6 = Hc
  setnames(data.sample,c(paste0("ba",startSim),paste0("h",startSim),
                         paste0("dbh",startSim),paste0("bl",startSim),
                         paste0("conif",startSim)),
           c("ba","h","dbh","bl","conif"))
  initVar <- array(NA, dim=c(nSites,7,nLayers))
  data.sample[,baConif:= (ba * conif/(conif+bl))]
  data.sample[,baBl:= (ba * bl/(conif+bl))]
  data.sample[,dbhConif:= dbh]
  data.sample[,dbhBl:= dbh]
  data.sample[,hConif:= h]
  data.sample[,hBl:= h]
  
  data.sample[,N:=ba/(pi*(dbh/2)^2/10000)]
  
  areas <- data.sample$area
  
  initVar[,1,] <- as.numeric(rep(c(1,3),each=nSites))
  initVar[,2,] <- as.numeric(data.sample[,h])*3.3  # round(as.numeric(data.sample[,age]))  ##### set to 1 because we do not know age
  initVar[,3,] <- as.numeric(data.sample[,h])
  # initVar[,3,][which(initVar[,3,]<1.5)] <- 1.5  ####if H < 1.5 set to 1.5
  
  initVar[,4,] <- as.numeric(data.sample[,dbh])
  if(domSPrun==1){
    ##initialize model only for dominant species##
    initVar[,5,] = 0.
    ix = unlist(data.sample[, which.max(c(pineP, spruceP, blp)), by=1:nrow(data.sample)] [, 2])
    for(jx in 1:nSites) initVar[jx,5,ix[jx]] = as.numeric(data.sample[, ba])[jx]
  } else{
    ###initialize model for mixed forest runs
    initVar[,5,1] <- as.numeric(data.sample[,(ba * conif/(conif+bl))])
    initVar[,5,2] <- as.numeric(data.sample[,(ba * bl/(conif+bl))])
    
    # !!To check start
    # if(varHD){ #### if true will vary H and D of pine and spruce using siteType
    #   ###increase spruceP dbh 10% for spruceP sitetype 1:2
    #   data.sample[conif>0. & spruceP >0. & siteType<2.5 & baSP > baP,X:=(ba-1.1*baSP-baB)/baP]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,dbhP:=X*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,dbhSP:=1.1*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5  & baSP > baP & dbhP<0.5,dbhSP:=((ba-(0.5/dbh)*baP-baB)/baSP)*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5  & baSP > baP & dbhP<0.5,dbhP:=0.5]
    #   
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,dbhSP:=dbh * (ba - 0.9*baP - baB)/baSP]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,dbhP:=pmax(0.9*dbh,0.3)]
    #   
    #   ####increase spruceP h 10% for spruceP sitetype 1:2
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,X:=(ba-1.1*baSP-baB)/baP]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,hP:=X*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,hSP:=1.1*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP & hP<1.5,hSP:=((ba-(1.5/h)*baP-baB)/baSP)*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP & hP<1.5,hP:=1.5]
    #   
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,hSP:=h * (ba - 0.9*baP - baB)/baSP]
    #   data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,hP:=pmax(0.9*h,1.3)]
    #   
    #   ####increase spruceP dbh 5% for spruceP sitetype 3
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP,X:=(ba-1.05*baSP-baB)/baP]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP,dbhP:=X*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP,dbhSP:=1.05*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP & dbhP<0.5,dbhSP:=((ba-(0.5/dbh)*baP-baB)/baSP)*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP & dbhP<0.5,dbhP:=0.5]
    #   
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,dbhSP:=dbh * (ba - 0.95*baP - baB)/baSP]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,dbhP:=pmax(0.95*dbh,0.3)]
    #   
    #   ####increase spruceP h 5% for spruceP sitetype 3
    #   data.sample[pineP>0. & spruceP >0. & siteType==3,X:=(ba-1.05*baSP-baB)/baP]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3,hP:=X*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3,hSP:=1.05*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & hP<1.5,hSP:=((ba-(1.5/h)*baP-baB)/baSP)*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & hP<1.5,hP:=1.5]
    #   
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,hSP:=h * (ba - 0.95*baP - baB)/baSP]
    #   data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,hP:=pmax(0.95*h,1.3)]
    #   
    #   ####increase pineP dbh 10% for sitetype >= 4
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,X:=(ba-1.1*baP-baB)/baSP]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhSP:=X*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhP:=1.1*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & dbhSP<0.5,dbhP:=((ba-(0.5/dbh)*baSP-baB)/baP)*dbh]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & dbhSP<0.5,dbhSP:=0.5]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhP:=dbh * (ba - 0.9*baSP - baB)/baP]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhSP:=pmax(0.9*dbh,0.3)]
    #   ####increase pineP h 10% for sitetype >= 4
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,X:=(ba-1.1*baP-baB)/baSP]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hSP:=X*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hP:=1.1*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & hSP<1.5,hP:=((ba-(1.5/h)*baSP-baB)/baP)*h]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & hSP<1.5,hSP:=1.5]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hP:=h * (ba - 0.9*baSP - baB)/baP]
    #   data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hSP:=pmax(0.9*h,1.3)]
    #   initVar[,3,1] <- as.numeric(data.sample[,hP])
    #   initVar[,3,2] <- as.numeric(data.sample[,hSP])
    #   initVar[,4,1] <- as.numeric(data.sample[,dbhP])
    #   initVar[,4,2] <- as.numeric(data.sample[,dbhSP])
      
    # }
    # !!To check end
  }
  
  # initVar[,6,] <- as.numeric(data.sample[,hc])
  
  ###check which BA ==0. and set to 0 the rest of the variables
  NoConif <- which(initVar[,5,1]==0.)
  NoDecid <- which(initVar[,5,2]==0.)
  
  siteInfo[NoConif,8:9] <- siteInfo[NoConif,8:9] - 1
  siteInfo[NoDecid,8:9] <- siteInfo[NoDecid,8:9] - 1
  
  #siteInfo[NoPine,4] <- siteInfo[NoPine,4] - 1
  #siteInfo[NoSpruce,4] <- siteInfo[NoSpruce,4] - 1
  #siteInfo[NoDecid,4] <- siteInfo[NoDecid,4] - 1
  initVar[NoConif,3:6,1] <- 0.
  initVar[NoDecid,3:6,2] <- 0.
  initVar[NoConif,,1] <- initVar[NoConif,,2]
  # initVar[NoPine,,1:2] <- initVar[NoPine,,2:3]
  
  nLay1 <- which(siteInfo[,8]==1)
  # nLay2 <- which(siteInfo[,8]==2)
  initVar[nLay1,3:7,2] <- 0
  initVar[,7,] <- 0
  # initVar[nLay2,3:7,3] <- 0
  # initVar[which(initVar[,5,1]==0.),,1] <- initVar[which(initVar[,5,1]==0.),,2]
  # initVar[which(initVar[,5,1]==0.),,2] <- initVar[which(initVar[,5,1]==0.),,3]
  # initVar[which(initVar[,5,1]==0.),1,3] <- 1
  # initVar[which(initVar[,5,1]==0.),3:6,3] <- 0
  
  # if (FALSE) {
  #   dat = dat[climID %in% data.sample[, unique(climID)]]
  #   
  #   if(weather!= "CurrClim.rdata"){
  #     # dat[, pvm:= as.Date('1980-01-01') - 1 + rday ]
  #     # dat[, DOY:= as.numeric(format(pvm, "%j"))]
  #     dat[, Year:= as.numeric(floor(rday/366)+1971)]
  #     dat = dat[Year >= startSim]
  #     dat[DOY==366, DOY:=365]
  #   }
  #   PARtran = t( dcast(dat[, list(climID, rday, PAR)], rday ~ climID,
  #                      value.var="PAR")[, -1])
  #   TAirtran = t( dcast(dat[, list(climID, rday, TAir)], rday ~ climID,
  #                       value.var="TAir")[, -1])
  #   VPDtran = t( dcast(dat[, list(climID, rday, VPD)], rday ~ climID,
  #                      value.var="VPD")[, -1])
  #   Preciptran = t( dcast(dat[, list(climID, rday, Precip)], rday ~ climID,
  #                         value.var="Precip")[, -1])
  #   CO2tran = t( dcast(dat[, list(climID, rday, CO2)], rday ~ climID,
  #                      value.var="CO2")[, -1])
  # }
  # siteInfo[, 2]  = match(as.numeric(siteInfo[, 2]), as.numeric(rownames(clim[[1]])))
  # siteInfo[, 2]  = match(siteInfo[,2], unique(dat$climID))
  
  # defaultThin=as.numeric(1-data.sample[, cons])
  # energyCut <- ClCut <- as.numeric(1-data.sample[, cons])
  ## Set to match climate data years
  initPrebas <- InitMultiSite(nYearsMS = rep(nYears,nSites),siteInfo=siteInfo,
                              # litterSize = litterSize,#pAWEN = parsAWEN,
                              defaultThin = defaultThin,
                              ClCut = ClCut,
                              areas =areas,
                              pCROBAS = pX,
                              # energyCut = energyCut, 
                              multiInitVar = as.array(initVar),
                              PAR = clim$PAR[, 1:(nYears*365)],
                              TAir=clim$TAir[, 1:(nYears*365)],
                              VPD=clim$VPD[, 1:(nYears*365)],
                              Precip=clim$Precip[, 1:(nYears*365)],
                              CO2=clim$CO2[, 1:(nYears*365)],
                              yassoRun = 1)
  initPrebas
}




yasso.mean.climate.f = function(dat, data.sample, startingYear, nYears){
  dat = dat[climID %in% data.sample[, unique(climID)]]
  dat[, DOY:=rep(1:365, len=dim(dat)[1])]
  dat[, Year:=rep(1980:2099, each=365)]
  #dat[, Year:= as.numeric(format(pvm, "%Y"))]
  dat = dat[Year >= startingYear & Year <= startingYear+nYears]
  dat[, pvm:= as.Date(paste(Year, '-01-01', sep="")) - 1 + DOY ]
  #dat[, DOY:= as.numeric(format(pvm, "%j"))]
  dat[, Mon:= as.numeric(format(pvm, "%m"))]
  #dat[DOY==366, DOY:=365]
  Tmean = dat[, mean(TAir), by = Year]
  Tsum = dat[, sum(ifelse(TAir>5, TAir-5, 0)), by=.(climID, Year)][, mean(V1), by=Year]
  PAR = dat[, mean(PAR), by = Year]
  VPD = dat[, mean(VPD), by = Year]
  CO2 = dat[, mean(CO2), by = Year]
  Precip = dat[, sum(Precip), by = .(climID, Year)][, mean(V1), by=Year]
  Tampl = dat[, .(mean(TAir)), by = .(climID, Year, Mon)][, (max(V1)-min(V1))/2, by=Year]
  
  out = cbind(Tmean, Precip[, -1], Tampl[, -1], CO2[, -1], PAR[, -1], VPD[, -1], Tsum[, -1])
  colnames(out) = c('Year','Tmean','Precip','Tampl', 'CO2', "PAR", "VPD", "Tsum5")
  out
}


prep.climate.f = function(dat, data.sample, startingYear, nYears,startYearWeather){
  dat = dat[climID %in% data.sample[, unique(climID)]]
  if(weather == "CurrClim"){
    dat[, Year:= as.numeric(floor((rday-1)/365)+startYearWeather)]
    dat1 = dat[Year >= startingYear]
    if(nYears>length(unique(dat1$Year))){
      nSampleYear <- nYears - length(unique(dat1$Year))
      set.seed(123)
      yearX <- sample(startYearWeather:min(startingYear,max(dat$Year)),nSampleYear,replace = F)
      lastYear <- max(dat$Year,startingYear)
      newYears <- lastYear + 1:length(yearX)
      dat2 <- dat[Year %in% yearX,]
      dat2$Year <- newYears[match(dat2$Year,yearX)]
      dat <- rbind(dat1,dat2)
      dat[,DOY:=rep(1:365,nYears),by=climID]
      setorder(dat,climID,Year,DOY)
      dat[,rday:=1:(365*nYears),by=climID]
    }
    
  }else{
    dat[, pvm:= as.Date('1980-01-01') - 1 + rday ]
    dat[, DOY:= as.numeric(format(pvm, "%j"))]
    dat[, Year:= as.numeric(format(pvm, "%Y"))]
    dat = dat[Year >= startingYear]
    dat[DOY==366, DOY:=365]
  }
  climID = dat[,unique(climID)]
  PARtran = t( dcast(dat[, list(climID, rday, PAR)], rday ~ climID,
                     value.var="PAR")[, -1])
  TAirtran = t( dcast(dat[, list(climID, rday, TAir)], rday ~ climID,
                      value.var="TAir")[, -1])
  VPDtran = t( dcast(dat[, list(climID, rday, VPD)], rday ~ climID,
                     value.var="VPD")[, -1])
  Preciptran = t( dcast(dat[, list(climID, rday, Precip)], rday ~ climID,
                        value.var="Precip")[, -1])
  CO2tran = t( dcast(dat[, list(climID, rday, CO2)], rday ~ climID,
                     value.var="CO2")[, -1])
  list(PAR=PARtran, TAir=TAirtran, VPD=VPDtran, 
       Precip=Preciptran, CO2=CO2tran,climID=climID)
}


##function to compile all data and create data.table 
createDT <- function(climate, management,variable, layer,startingYear,funX,siteTypeX){
  
  loadFolder <- paste0(outPath,"init",startingYear,"/",
                       "st",siteTypeX,"/")
  files <- list.files(path= loadFolder)#,pattern = paste0("year",startingYear,"_"))
  startV <- data.table()
  for (ij in variable) assign(varNames[ij],data.table())
  
  for(i in 1:length(files)){
    sampleID <- paste0("sample",i,".")
    
    fileX <- files[grep(sampleID,files,fixed = T)]
    
    load(paste0(loadFolder,fileX))
    
    if(exists("v0")) startV <- rbind(startV,v0)
    margin= 1:2#(length(dim(out$annual[,,variable,]))-1)
    if(layer=="all"){
      ijX <- 0
      for (ij in variable){ 
        ijX = ijX+1
        if(varNames[ij] %in% dimnames(out)$varX){
          varX <- which(dimnames(out)$varX==varNames[ij])
          if(funX[ijX]=="sum") assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                                                    apply(out[,,varX,],margin,sum,na.rm=T))))
          
          if(funX[ijX]=="mean"){
            BAind <- which(saveVars==13)
            if(length(BAind)<1){
              print(paste("BA not saved!! for",varNames[saveVars[varX]],"aritmetic mean was used."))
              assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                                   apply(out[,,varX,],margin,mean,na.rm=T))))
            }else{
              sumBA <- apply(out[,,BAind,],1:2,sum,na.rm=T)
              sumBAs <- array(sumBA,dim = dim(out[,,BAind,]))
              fracBA <- out[,,BAind,]/sumBAs
              wgdVar <- out[,,varX,] * fracBA
              assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                                   apply(wgdVar,margin,sum,na.rm=T))))
            } }
        }else{
          print(paste(varNames[ij],"not saved"))
        }
      }
    }else{
      for(ij in variable){
        if(varNames[ij] %in% dimnames(out)$varX){
          varX <- which(dimnames(out)$varX==varNames[ij])
          assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                               out[,,varX,layer])))
        }else{
          print(paste(varNames[ij],"not saved"))
        }
      } 
    }
    print(i)
  }
  
  for(ij in variable){
    file2save <- paste0("outDT/init",startingYear,"/","st",siteTypeX,"/",varNames[ij],"_",management,"_",climate,
                        "layer",layer,".rdata")
    save(list=varNames[ij],file=file2save)
  }
  if((dim(startV)[1]) > 0 & siteTypeX==startingYear) save(startV,file=paste0("outDT/","init",startingYear,"/","startV_layerall.rdata"))
}


# this function create raster in tif format from data.tables selecting one year or the average of a time priod if yearOut is a vector of years
createTifFromDT <- function(yearOut, varX, stYear,XYsegID,crsX=NA){
  simYear <- yearOut - stYear
  fileDT=paste0("outDT/","init",startingYear,"/st",siteTypeX,"/",varNames[varX],"_",management,"_",climate,
                "layer",layerDT,".rdata")
  load(fileDT)
  
  outX <- t(get(varNames[varX]))
  if (length(simYear)==1) outX <- outX[simYear,]
  if (length(simYear)>1) outX <- colMeans(outX[simYear,],na.rm = T)
  # outX <- data.table(cbind(segID,areas,outX))
  outX <- data.table(cbind(segID,outX))
  
  setnames(outX,c("segID",varNames[varX]))
  
  setkey(XYsegID,segID)
  setkey(outX,segID)
  outXY <- merge(XYsegID,outX,all = T)
  ###remove coordinates NAs
  outXY <- outXY[!is.na(x)]
  
  ###create raster 
  rastX <- rasterFromXYZ(outXY[,c("x","y",varNames[varX]),with=F])
  crs(rastX) <- crsX
  
  rastName <- paste0("outRast/","init",startingYear,"/st",siteTypeX,"/",climate,"_",management,"_var",varNames[varX],
                     "_spec",layerDT,"_yearStart",stYear,"_yearOut",yearOut,".tif")
  writeRaster(rastX,filename = rastName,overwrite=T)
}




####functions used in the siteType calculations
fixBAper <- function(BApers){
  minBA <- min(BApers)
  if(minBA<0) BApers <- BApers - minBA
  BApers <- BApers/sum(BApers)*100
  return(BApers)
}

###function for site type data assimilation
pSTx <- function(segIDx,nSample,year1,year2,tileX){
  mu1 <- errData[[paste0("y",year1)]][[paste0("t",tileX)]]$muFSVda
  if(is.null(mu1)) mu1 <- errData[[paste0("y",year1)]]$all$muFSVda
  sigma1 <- errData[[paste0("y",year1)]][[paste0("t",tileX)]]$sigmaFSVda
  if(is.null(sigma1)) sigma1 <- errData[[paste0("y",year1)]]$all$sigmaFSVda
  mu2 <- errData[[paste0("y",year2)]][[paste0("t",tileX)]]$muSTda
  if(is.null(mu2)) mu2 <- errData[[paste0("y",year2)]]$all$muSTda
  sigma2 <- errData[[paste0("y",year2)]][[paste0("t",tileX)]]$sigmaSTda
  if(is.null(sigma2)) sigma2 <- errData[[paste0("y",year2)]]$all$sigmaSTda
  set.seed(1234)
  sampleError <- data.table(mvrnorm(nSample,mu=mu1,Sigma=sigma1))
  # segIDx <- dataSurV[segID==2]
  sampleX <- data.table()
  sampleX$H <- segIDx$H + sampleError$H
  sampleX$D <- segIDx$D + sampleError$D
  sampleX$BAtot <- segIDx$BAtot + sampleError$G
  sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
  sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
  sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
  sampleX$BAp <- segIDx$BApPer * sampleX$BAtot/100
  sampleX$BAsp <- segIDx$BAspPer * sampleX$BAtot/100
  sampleX$BAb <- segIDx$BAbPer * sampleX$BAtot/100
  
  ###filter data
  minH <- 1.5; minD <- 0.5; minB <- pi*(minD/2)^2/10000*2200
  xx <- unique(c(which(sampleX$H< minH),which(sampleX$D< minD),which(sampleX$BAtot< minB)))
 
  sampleX[, c("BApPer", "BAspPer", "BAbPer"):=
            as.list(fixBAper(unlist(.(BApPer,BAspPer,BAbPer)))), 
          by = seq_len(nrow(sampleX))]
  
  max.pro.est<-apply(segIDx[, c('BApPer','BAspPer','BAbPer')], 1, which.max)
  segIDx[,max.pro.est:=max.pro.est]
  
  if(year1=="all" & tileX=="all"){
    logistic.model <- logisticPureF$all
  }else{
    logistic.model <- logisticPureF[[paste0("y",year1)]][[paste0("t",tileX)]]
    # if(is.null(logistic.model)) logistic.model<- logisticPureF[[paste0("y",year1)]]$all
    if(is.null(logistic.model)) logistic.model<- logisticPureF$all
  }

  # tryCatch({
    set.seed(1234)
    sampleX$pureF <- runif(min(nSample,nrow(sampleX)),0,1)<predict(logistic.model,type="response",newdata = segIDx)
  # }, error=function(e){cat("ERROR logisticPureF:",print(segIDx), "\n")})

  if(max.pro.est==1) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(100,0,0)]
  if(max.pro.est==2) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,100,0)]
  if(max.pro.est==3) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,0,100)]
  
  sampleX[,BAp:=BApPer*BAtot/100]
  sampleX[,BAsp:=BAspPer*BAtot/100]
  sampleX[,BAb:=BAbPer*BAtot/100]
  # sampleX[,st:=segIDx$st]
  sampleX[,V2:=segIDx$V2]
  sampleX[,ba2:=segIDx$ba2]
  sampleX[,h2:=segIDx$h2]
  sampleX[,dbh2:=segIDx$dbh2]
  sampleX[,segID:=segIDx$segID]
  
  sampleX[,BAtot:=(BAp+BAsp+BAb)]
  sampleX[,BAh:=BAtot*H]
  sampleX[,N:=BAtot/(pi*(D/200)^2)]
  b = -1.605 ###coefficient of Reineke
  sampleX[,SDI:=N *(D/10)^b]
  
  # full.model<-lm(lnVmod~H+D+lnBAp+lnBAsp+lnBAb+st,data=dataX)
  # tryCatch(expr = {
    sampleX$st <- factor(1)
    sampleX[,VsurST1 := pmax(0.,predict(step.modelV,newdata=sampleX))]
    sampleX[,BsurST1 := pmax(0.,predict(step.modelB,newdata=sampleX))]
    sampleX[,DsurST1 := pmax(0.,predict(step.modelD,newdata=sampleX))]
    sampleX[,HsurST1 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  # },
  # error=function(e){cat("ERROR ST1:",print(segIDx), "\n")})
  
  ###estimates for negative values based on average growth  
  dH <- mean(sampleX[-xx]$HsurST1 - sampleX[-xx]$H,na.rm=T)
  dV <- mean(sampleX[-xx]$VsurST1 - segIDx$V,na.rm=T)
  dD <- mean(sampleX[-xx]$DsurST1 - sampleX[-xx]$D,na.rm=T)
  dB <- mean(sampleX[-xx]$BsurST1 - sampleX[-xx]$BAtot,na.rm=T)
  sampleX[xx]$HsurST1 <- sampleX[xx]$H + dH
  sampleX[xx]$DsurST1 <- sampleX[xx]$D + dD
  sampleX[xx]$VsurST1 <- segIDx$V + dV
  sampleX[xx]$BsurST1 <- sampleX[xx]$BAtot + dB
  
  sampleX$st <- factor(2)
  sampleX[,VsurST2 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST2 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST2 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST2 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  ###estimates for negative values based on average growth  
  dH <- mean(sampleX[-xx]$HsurST2 - sampleX[-xx]$H,na.rm=T)
  dV <- mean(sampleX[-xx]$VsurST2 - segIDx$V,na.rm=T)
  dD <- mean(sampleX[-xx]$DsurST2 - sampleX[-xx]$D,na.rm=T)
  dB <- mean(sampleX[-xx]$BsurST2 - sampleX[-xx]$BAtot,na.rm=T)
  sampleX[xx]$HsurST2 <- sampleX[xx]$H + dH
  sampleX[xx]$DsurST2 <- sampleX[xx]$D + dD
  sampleX[xx]$VsurST2 <- segIDx$V + dV
  sampleX[xx]$BsurST2 <- sampleX[xx]$BAtot + dB
  
  sampleX$st <- factor(3)
  sampleX[,VsurST3 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST3 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST3 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST3 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  ###estimates for negative values based on average growth  
  dH <- mean(sampleX[-xx]$HsurST3 - sampleX[-xx]$H,na.rm=T)
  dV <- mean(sampleX[-xx]$VsurST3 - segIDx$V,na.rm=T)
  dD <- mean(sampleX[-xx]$DsurST3 - sampleX[-xx]$D,na.rm=T)
  dB <- mean(sampleX[-xx]$BsurST3 - sampleX[-xx]$BAtot,na.rm=T)
  sampleX[xx]$HsurST3 <- sampleX[xx]$H + dH
  sampleX[xx]$DsurST3 <- sampleX[xx]$D + dD
  sampleX[xx]$VsurST3 <- segIDx$V + dV
  sampleX[xx]$BsurST3 <- sampleX[xx]$BAtot + dB
  
  sampleX$st <- factor(4)
  sampleX[,VsurST4 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST4 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST4 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST4 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  ###estimates for negative values based on average growth  
  dH <- mean(sampleX[-xx]$HsurST4 - sampleX[-xx]$H,na.rm=T)
  dV <- mean(sampleX[-xx]$VsurST4 - segIDx$V,na.rm=T)
  dD <- mean(sampleX[-xx]$DsurST4 - sampleX[-xx]$D,na.rm=T)
  dB <- mean(sampleX[-xx]$BsurST4 - sampleX[-xx]$BAtot,na.rm=T)
  sampleX[xx]$HsurST4 <- sampleX[xx]$H + dH
  sampleX[xx]$DsurST4 <- sampleX[xx]$D + dD
  sampleX[xx]$VsurST4 <- segIDx$V + dV
  sampleX[xx]$BsurST4 <- sampleX[xx]$BAtot + dB
  
  sampleX$st <- factor(5)
  sampleX[,VsurST5 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST5 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST5 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST5 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  ###estimates for negative values based on average growth  
  dH <- mean(sampleX[-xx]$HsurST5 - sampleX[-xx]$H,na.rm=T)
  dV <- mean(sampleX[-xx]$VsurST5 - segIDx$V,na.rm=T)
  dD <- mean(sampleX[-xx]$DsurST5 - sampleX[-xx]$D,na.rm=T)
  dB <- mean(sampleX[-xx]$BsurST5 - sampleX[-xx]$BAtot,na.rm=T)
  sampleX[xx]$HsurST5 <- sampleX[xx]$H + dH
  sampleX[xx]$DsurST5 <- sampleX[xx]$D + dD
  sampleX[xx]$VsurST5 <- segIDx$V + dV
  sampleX[xx]$BsurST5 <- sampleX[xx]$BAtot + dB
  
  dx1 <- cbind(sampleX$BsurST1 - segIDx$ba2,sampleX$DsurST1 - segIDx$dbh2,
               sampleX$HsurST1 - segIDx$h2,sampleX$VsurST1 - segIDx$V2)
  dx2 <- cbind(sampleX$BsurST2 - segIDx$ba2,sampleX$DsurST2 - segIDx$dbh2,
               sampleX$HsurST2 - segIDx$h2,sampleX$VsurST2 - segIDx$V2)
  dx3 <- cbind(sampleX$BsurST3 - segIDx$ba2,sampleX$DsurST3 - segIDx$dbh2,
               sampleX$HsurST3 - segIDx$h2,sampleX$VsurST3 - segIDx$V2)
  dx4 <- cbind(sampleX$BsurST4 - segIDx$ba2,sampleX$DsurST4 - segIDx$dbh2,
               sampleX$HsurST4 - segIDx$h2,sampleX$VsurST4 - segIDx$V2)
  dx5 <- cbind(sampleX$BsurST5 - segIDx$ba2,sampleX$DsurST5 - segIDx$dbh2,
               sampleX$HsurST5 - segIDx$h2,sampleX$VsurST5 - segIDx$V2)
  
  pst1 <- mean(dmvnorm(x=dx1, mean=mu2, sigma=sigma2, log=FALSE))
  pst2 <- mean(dmvnorm(x=dx2, mean=mu2, sigma=sigma2, log=FALSE))
  pst3 <- mean(dmvnorm(x=dx3, mean=mu2, sigma=sigma2, log=FALSE))
  pst4 <- mean(dmvnorm(x=dx4, mean=mu2, sigma=sigma2, log=FALSE))
  pst5 <- mean(dmvnorm(x=dx5, mean=mu2, sigma=sigma2, log=FALSE))
  psum <- pst1 +pst2+pst3 +pst4+pst5
  pst1 <- pst1/psum
  pst2 <- pst2/psum
  pst3 <- pst3/psum
  pst4 <- pst4/psum
  pst5 <- pst5/psum
  return(pST=c(segIDx$segID,pst1,pst2,pst3,pst4,pst5)) 
}


###function for structural variables data assimilation 
pSVDA <- function(segIDx,nSample,
                  errData1,errData2,
                  step.modelHx,step.modelDx,step.modelBx,
                  step.modelBconifx,step.modelBblx){
  mu1 <- errData1$muFSVda
  sigma1 <- errData1$sigmaFSVda
  mu2 <- errData2$muFSVda
  sigma2 <- errData2$sigmaFSVda
  
  #pST <- c(segIDx$pST1,segIDx$pST2,segIDx$pST3,segIDx$pST4,segIDx$pST5)
  pST <- c(0.1,0.25,0.3,0.25,0.1)
  st <- sample(rep(1:5,round(nSample*pST)),nSample,replace = T)
  set.seed(1234)
  
  muX <- mu1 + c(segIDx$BAtot, segIDx$D,segIDx$H,segIDx$BAconifPer,segIDx$BAblPer)
  if(any(muX<0)) muX[which(muX<0)] <- 1
  abc<-Fweibull.v(muX,diag(sigma1))
  shd<-data.frame(shape=abc$c,decay=abc$b^(-abc$c))
  rho<-cov2cor(sigma1)
  sampleX <- data.table(rmvweisd(nSample*10, shd$shape, shd$decay, rho))
  colnames(sampleX)<-c('BAtot','D','H','BAconifPer','BAblPer')
  
###filter data
  
  sampleX <- sampleX[H>1.5]
  sampleX <- sampleX[D>0.5]
  sampleX <- sampleX[BAtot>0.045]
  sampleX <- sampleX[H<35]
  sampleX <- sampleX[D<40]
  sampleX <- sampleX[BAtot<100]
  sampleX <- sampleX[BAconifPer<100]
  sampleX <- sampleX[BAblPer<100]
  sampleX <- sampleX[BAconifPer>0]
  sampleX <- sampleX[BAblPer>0]
  sampleX <- sampleX[1:nSample]
  sampleX$st <- st
  
  # sampleError <- data.table(mvrnorm(nSample,mu=mu1,Sigma=sigma1))
  
  # # segIDx <- dataSurV[segID==2]
  # sampleX <- data.table()
  # sampleX$H <- segIDx$H + sampleError$H
  # sampleX$D <- segIDx$D + sampleError$D
  # sampleX$BAtot <- segIDx$BAtot + sampleError$G
  # sampleX$BAconifPer <- segIDx$BAconifPer + sampleError$BAcon
  # sampleX$BAblPer <- segIDx$BAblPer + sampleError$BAbl
  # sampleX$BAconif <- segIDx$BAconifPer * sampleX$BAtot/100
  # sampleX$BAbl <- segIDx$BAblPer * sampleX$BAtot/100
  # 
  
  ###filter data
  # minH <- 1.5; minD <- 0.5; minB <- pi*(minD/2)^2/10000*2200
  # xx <- unique(c(which(sampleX$H< minH),which(sampleX$D< minD),which(sampleX$BAtot< minB)))

# to Check
  # sampleX[, c("BAconifPer", "BAblPer"):=
  #           as.list(fixBAper(unlist(.(BAconifPer,BAblPer)))), 
  #         by = seq_len(nrow(sampleX))]
  
  ###process pure forests!!!to check start
  # max.pro.est<-apply(segIDx[, c('BAconifPer','BAblPer')], 1, which.max)
  # segIDx$max.pro.est=max.pro.est
  # 
  # if(year1=="all" & tileX=="all"){
  #   logistic.model <- logisticPureF$all
  # }else{
  #   logistic.model <- logisticPureF[[paste0("y",year1)]][[paste0("t",tileX)]]
  #   # if(is.null(logistic.model)) logistic.model<- logisticPureF[[paste0("y",year1)]]$all
  #   if(is.null(logistic.model)) logistic.model<- logisticPureF$all
  # }
  # 
  # set.seed(1234)
  # sampleX$pureF <- runif(min(nSample,nrow(sampleX)),0,1)<predict(logistic.model,type="response",newdata = segIDx)
  # if(max.pro.est==1) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(100,0,0)]
  # if(max.pro.est==2) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,100,0)]
  # if(max.pro.est==3) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,0,100)]
  ###process pure forests!!!to check start
  
  sampleX[,BAconif:=BAconifPer*BAtot/100]
  sampleX[,BAbl:=BAblPer*BAtot/100]
  # sampleX[,st:=segIDx$st]
  sampleX[,V2:=segIDx$V2]
  sampleX[,ba2:=segIDx$ba2]
  sampleX[,h2:=segIDx$h2]
  sampleX[,dbh2:=segIDx$dbh2]
  # sampleX[,BApPer2:=segIDx$BApPer2]
  # sampleX[,BAspPer2:=segIDx$BAspPer2]
  # sampleX[,BAbPer2:=segIDx$BAbPer2]
  # sampleX[,segID:=segIDx$segID]
  
  # sampleX$lnVmod<-log(sampleX$Vmod)
  # sampleX$st<-factor(sampleX$st,levels = 1:5)     ##!!!!Xianglin
  # sampleX$st <- factor(sampleX$st)
  sampleX[,BAtot:=(BAconif+BAbl)]
  sampleX[,BAh:=BAtot*H]
  sampleX[,N:=pmin(20000,BAtot/(pi*(D/200)^2))]
  b = -1.605 ###coefficient of Reineke
  sampleX[,SDI:=N *(D/10)^b]
  sampleX$st <- st
  sampleX[,rootBAconif:=BAconif^0.5]
  sampleX[,BAconif2:=ba2*BAconifPer/100]
  sampleX[,BAbl2:=ba2*BAblPer/100]
  
  # full.model<-lm(lnVmod~H+D+lnBAp+lnBAsp+lnBAb+st,data=dataX)
  sampleX$st <- factor(sampleX$st)
  sampleX[,Hx := pmax(0.,predict(step.modelHx,newdata=sampleX))]
  sampleX[,Dx := pmax(0.,predict(step.modelDx,newdata=sampleX))]
  sampleX[,Bx := pmax(0.,predict(step.modelBx,newdata=sampleX))]
  # sampleX[,Vx := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,Bconifx := pmax(0.,predict(step.modelBconifx,newdata=sampleX))]
  sampleX[,Bblx := pmax(0.,predict(step.modelBblx,newdata=sampleX))]
  sampleX[,BconifPerx := Bconifx/Bx*100]
  sampleX[,BblPerx := Bblx/(Bx)*100]
  # sampleX[,rootBAp:=BAp^0.5]
  
  # if(length(xx)>0){
  #   ###estimates for negative values based on average growth  
  #   dH <- mean(sampleX[-xx]$Hx - sampleX[-xx]$H,na.rm=T)
  #   dD <- mean(sampleX[-xx]$Dx - sampleX[-xx]$D,na.rm=T)
  #   dB <- mean(sampleX[-xx]$Bx - sampleX[-xx]$BAtot,na.rm=T)
  #   sampleX[xx]$Hx <- sampleX[xx]$H + dH
  #   sampleX[xx]$Dx <- sampleX[xx]$D + dD
  #   sampleX[xx]$Bx <- sampleX[xx]$BAtot + dB
  #   sampleX[xx]$BconifPerx <- sampleX[xx]$BAconifPer
  #   sampleX[xx]$BblPerx <- sampleX[xx]$BAblPer
  # }
  
  nax <- unique(which(is.na(sampleX[,.(Hx,Dx,Bx,BconifPerx,BblPerx)]),arr.ind = T)[,1])
  if(length(nax)>0) sampleX <- sampleX[-nax]
  mux <- sampleX[,colMeans(cbind(Hx,Dx,Bx,BconifPerx,BblPerx))]
  sigmax <- sampleX[,cov(cbind(Hx,Dx,Bx,BconifPerx,BblPerx))]
  
  pMvnormx <- c(mux,as.vector(sigmax))
  
  ###sample from second measurement
  mux2 <- c(segIDx$h2 + mu2[3],segIDx$dbh2 + mu2[2],segIDx$BAtot2 + mu2[1],
            segIDx$BAconifPer2 + mu2[4],segIDx$BAblPer2 + mu2[5])
  sigmax2 <- sigma2[c(3,2,1,4,5),c(3,2,1,4,5)]
  
  pMvnormx2 <- c(mux2,as.vector(sigmax2))
  
  
  LL <- solve(sigmax + sigmax2,tol=1e-20)
  sigmaPost <- sigmax %*% LL %*% sigmax2
  muPost <- sigmax2 %*% LL %*% mux + 
    sigmax %*% LL %*% mux2
  
  pMvnormPost <- c(muPost,as.vector(sigmaPost))
  # ss= inv(inv(sigmax)+inv(sigmax2))
  # ff <- sigmaPost %*% (sigmax %^%(-1)) %*% mux + sigmaPost %*% (sigmax2 %^%(-1)) %*% mux2
  
  # return(list(muPrior=mux,muLik=mux2,muPost=as.vector(muPost),
  #             sigPrior=sigmax,sigLik=sigmax2,sigPost=sigmaPost))
  pars <- c(as.vector(pMvnormx),as.vector(pMvnormx2),as.vector(pMvnormPost))
  return(pars)
}


###function for structural variables forcast and uncertainty 
prForUnc <- function(segIDx,nSample,yearUnc,tileX){
  if(yearUnc=="all"){
    muUnc <- errData$all$muFSVda
    sigmaUnc <- errData$all$sigmaFSVda
    step.probitX <- step.probit$all
  }else{
    muUnc <- errData[[paste0("y",yearUnc)]][[paste0("t",tileX)]]$muFSVda
    sigmaUnc <- errData[[paste0("y",yearUnc)]][[paste0("t",tileX)]]$sigmaFSVda
    step.probitX <- step.probit[[paste0("y",yearUnc)]][[paste0("t",tileX)]]
  }
  
  pST <- predict(step.probitX,type='p',segIDx)   ### needs to be changed . We need to calculate with 2016 and 2019 data
  
  st <- sample(rep(1:5,round(nSample*pST)),nSample,replace = T)
  set.seed(1234)
  sampleError <- data.table(mvrnorm(nSample,mu=muUnc,Sigma=sigmaUnc))
  
  # segIDx <- dataSurV[segID==2]
  sampleX <- data.table()
  sampleX$H <- segIDx$H + sampleError$H
  sampleX$D <- segIDx$D + sampleError$D
  sampleX$BAtot <- segIDx$BAtot + sampleError$G
  sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
  sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
  sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
  sampleX$BAp <- segIDx$BApPer * sampleX$BAtot/100
  sampleX$BAsp <- segIDx$BAspPer * sampleX$BAtot/100
  sampleX$BAb <- segIDx$BAbPer * sampleX$BAtot/100
  ###filter data
  minH <- 1.5; minD <- 0.5; minB <- pi*(minD/2)^2/10000*2200
  xx <- unique(c(which(sampleX$H< minH),which(sampleX$D< minD),which(sampleX$BAtot< minB)))
  sampleX[,rootBAp:=BAp^0.5]
  # sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  # if(nrow(sampleX)<nSample){
  #   sample1 <- sampleX
  #   set.seed(123)
  #   sampleError <- data.table(mvrnorm(nSample*2,mu=muUnc,Sigma=sigmaUnc))
  #   # segIDx <- dataSurV[segID==2]
  #   sampleX <- data.table()
  #   sampleX$H <- segIDx$H + sampleError$H
  #   sampleX$D <- segIDx$D + sampleError$D
  #   sampleX$BAtot <- segIDx$BAtot + sampleError$G
  #   sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
  #   sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
  #   sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
  #   sampleX$BAp <- segIDx$BApPer * sampleX$BAtot/100
  #   sampleX$BAsp <- segIDx$BAspPer * sampleX$BAtot/100
  #   sampleX$BAb <- segIDx$BAbPer * sampleX$BAtot/100
  #   sampleX[H<1.3]$H <- 1.3
  #   sampleX[D<0]$D <- 0.1
  #   sampleX[BAtot<0]$BAtot <- 0.01
  #   sampleX[,rootBAp:=BAp^0.5]
  #   sampleX <- rbind(sample1,sampleX)
  #   sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  # }
  
  sampleX[, c("BApPer", "BAspPer", "BAbPer"):=
            as.list(fixBAper(unlist(.(BApPer,BAspPer,BAbPer)))), 
          by = seq_len(nrow(sampleX))]
  
  max.pro.est<-apply(segIDx[, c('BApPer','BAspPer','BAbPer')], 1, which.max)
  segIDx$max.pro.est=max.pro.est
  
  if(yearUnc=="all"){
    logistic.model <- logisticPureF$all
  }else{
    logistic.model <- logisticPureF[[paste0("y",yearUnc)]][[paste0("t",tileX)]]
  }
  
  set.seed(1234)
  sampleX$pureF <- runif(min(nSample,nrow(sampleX)),0,1)<predict(logistic.model,type="response",newdata = segIDx)
  if(max.pro.est==1) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(100,0,0)]
  if(max.pro.est==2) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,100,0)]
  if(max.pro.est==3) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,0,100)]
  
  sampleX[,BAp:=BApPer*BAtot/100]
  sampleX[,BAsp:=BAspPer*BAtot/100]
  sampleX[,BAb:=BAbPer*BAtot/100]
  
  # sampleX[,segID:=segIDx$segID]
  
  sampleX[,BAtot:=(BAp+BAsp+BAb)]
  sampleX[,BAh:=BAtot*H]
  sampleX[,N:=BAtot/(pi*(D/200)^2)]
  b = -1.605 ###coefficient of Reineke
  sampleX[,SDI:=N *(D/10)^b]
  sampleX$st <- st
  
  sampleX$st <- factor(sampleX$st)
  sampleX[,Hx := pmax(0.,predict(step.modelH,newdata=sampleX))]
  sampleX[,Dx := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,Bx := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,Bpx := pmax(0.,predict(step.modelBp,newdata=sampleX))]
  sampleX[,Bspx := pmax(0.,predict(step.modelBsp,newdata=sampleX))]
  sampleX[,Bdx := pmax(0.,predict(step.modelBd,newdata=sampleX))]
  sampleX[,BpPerx := Bpx/(Bpx+Bspx+Bdx)*100]
  sampleX[,BspPerx := Bspx/(Bpx+Bspx+Bdx)*100]
  sampleX[,BdPerx := Bdx/(Bpx+Bspx+Bdx)*100]
  
  ##estimates for negative values based on average growth  
  dH <- mean(sampleX[-xx]$Hx - sampleX[-xx]$H,na.rm=T)
  dD <- mean(sampleX[-xx]$Dx - sampleX[-xx]$D,na.rm=T)
  dB <- mean(sampleX[-xx]$Bx - sampleX[-xx]$BAtot,na.rm=T)
  sampleX[xx]$Hx <- sampleX[xx]$H + dH
  sampleX[xx]$Dx <- sampleX[xx]$D + dD
  sampleX[xx]$Bx <- sampleX[xx]$BAtot + dB
  sampleX[xx]$BpPerx <- sampleX[xx]$BApPer
  sampleX[xx]$BspPerx <- sampleX[xx]$BAspPer
  sampleX[xx]$BdPerx <- sampleX[xx]$BAbPer
  
  nax <- unique(which(is.na(sampleX),arr.ind = T)[,1])
  if(length(nax)>0) sampleX <- sampleX[-nax]
  mux <- sampleX[,colMeans(cbind(Hx,Dx,Bx,BpPerx,BspPerx,BdPerx))]
  sigmax <- sampleX[,cov(cbind(Hx,Dx,Bx,BpPerx,BspPerx,BdPerx))]
  
  pMvnormx <- c(mux,as.vector(sigmax))
  
  ###sample from second measurement
  pars <- as.vector(pMvnormx)
  return(pars)
}




####distribution function to correct the distribution
distr_correction <- function(Y,Xtrue,Xdiscr=FALSE){
  # y is the matrix of new sample - observations in lines, variables in columns
  # x is the matrix of the measured values = truth
  X <- Xtrue[sample(1:nrow(Xtrue),nrow(Xtrue),replace = TRUE),]
  m <- ncol(X)
  n <- nrow(X)
  ny <- nrow(Y)
  
  # Create cdf for each marginal distr. and use it to generate 
  # ny random variables from that distribution:
  Ynew <- matrix(0,nrow=ny,ncol=m)
  
  # Do post-processing variable by variable
  for(j in 1:m){
    xu <- unique(X[,j]) # choose the unique values in sample
    xu <- xu[order(xu)] # and sort them from smallest to biggest
    # create empirical cdf:
    cdf <- matrix(0,length(xu),1) # empty matrix
    for(i in 1:length(xu)){
      cdf[i] <- sum(X[,j] <= xu[i]) # count the sample values less than 
      # given value,
    }
    cdf <- cdf/(n+1) # cdf goes from 0 to 1, here are the inner points
    
    r <- matrix(ny,1)
    for(i in 1:nrow(Y)){
      r[i] <- sum(Y[,j]<=Y[i,j])/ny # Define the eCDF of each observation 
    }
    # use discr. or cont. inverse cdf to generate new sample
    if(Xdiscr[j] | length(Xdiscr)<m){ # discrete values 
      #interp <- approx(c(0,cdf), c(xu,max(xu)+1), r, method = "constant", yleft = xu[1], 
      #                 yright = xu[length(xu)], rule = 2:1)
      interp <- approx(c(0,cdf,1), c(min(xu),xu,max(xu)+1), r, method = "constant", yleft = xu[1], 
                       yright = xu[length(xu)], rule = 2:1)
    } else { # continuous values 
      #interp <- approx(c(0,cdf), c(xu,max(xu)+1), r, yleft = xu[1], 
      #                 yright = xu[length(xu)], rule = 2:1)
      interp <- approx(c(0,cdf,1), c(min(xu),xu,max(xu)*1.01), r, yleft = xu[1], 
                       yright = xu[length(xu)], rule = 2:1)
    }
    Ynew[,j] <- interp$y
  }
  return(Ynew)
}


library(mvtnorm)
#Function 1: random generation function of multivariate Weibull distribution
rmvweisd <- function(n, shape=1, decay=1, corr=diag(length(shape)))
{
  ## extract parameters, do sanity checks, deal with univariate case
  
  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)
  
  Ds = length(shape)
  if(Ds > D)
    warning("'shape' longer than width of 'corr', truncating to fit")
  if(Ds != D)
    shape = rep(shape, length.out=D)
  
  Dd = length(decay)
  if(Dd > D)
    warning("'decay' longer than width of 'corr', truncating to fit")
  if(Dd != D)
    decay = rep(decay, length.out=D)
  
  if(D == 1) rweisd(n, shape, decay)
  
  ## generate standard multivariate normal matrix, convert to CDF
  Z = rmvnorm(n, mean=rep(0,dim(corr)[1]), sigma = corr)
  cdf = pnorm(Z)
  
  ## convert to Weibull (WeiSD), return
  sapply(1:D, function(d) qweisd(cdf[,d], shape[d], decay[d]))
}

# Here we use rmvnorm function from mvtnorm package directly
# #Function 2: random generation function of multivariate normal distribtution
# rmvnormal <- function(n, mean=NULL, cov=NULL)
# {
# 	## munge parameters PRN and deal with the simplest univariate case
# 
# 	if(is.null(mean))
# 		if(is.null(cov))
# 			return(rnorm(n))
# 		else
# 			mean = rep(0, nrow(cov))
# 	else if (is.null(cov))
# 		cov = diag(length(mean))
# 
# 	## gather statistics, do sanity checks
# 
# 	D = length(mean)
# 	if (D != nrow(cov) || D != ncol(cov))
# 		stop("length of mean must equal nrow and ncol of cov")
# 
# 	E = eigen(cov, symmetric=TRUE)
# 	if (any(E$val < 0)){
# 	  stop("Numerically negative definite covariance matrix")
# 	  warning("Numerically negative definite covariance matrix")}
# 
# 	## generate values and return
# 	mean.term = mean
# 	covariance.term = E$vec %*% (t(E$vec) * sqrt(E$val))
# 	independent.term = matrix(rnorm(n*D), nrow=D)
# 
# 	drop(t(mean.term + (covariance.term %*% independent.term)))
# }

#Function 3: quantile function of Weibull distribution
qweisd <- function(p, shape, decay, lower.tail=TRUE, log.p=FALSE)
  qweibull(p, shape=shape, scale=decay^(-1/shape),
           lower.tail=lower.tail, log.p=log.p)


Fweibull<-function(dm,pdvar){
  
  a=0
  pdvar
  
  c = 0.1
  c1 = c
  g1 = gamma(1 + 1/c)
  g2 = gamma(1 + 2/c)
  # b = ((-a*g1)/g2) + ((((a/g2)**2)*(g2*g2-g1) + (dm*dm/g2))**0.5)
  b = (dm-a)/g1
  fOld = b*b * (g2 - g1**2) - pdvar
  
  c = c1 + 1
  g1 = gamma(1 + 1/c)
  g2 = gamma(1 + 2/c)
  # b = ((-a*g1)/g2) + ((((a/g2)**2)*(g2*g2-g1) + (dm*dm/g2))**0.5)
  b = (dm-a)/g1
  fNew = b*b * (g2 - g1**2) - pdvar
  fOld2 = fNew
  
  while (fOld*fNew > 0)
  {
    c1 = c
    c = c1 + 1
    g1 = gamma(1 + 1/c)
    g2 = gamma(1 + 2/c)
    # b = ((-a*g1)/g2) + ((((a/g2)**2)*(g2*g2-g1) + (dm*dm/g2))**0.5)
    b = (dm-a)/g1
    fNew = b*b * (g2 - g1**2) - pdvar
    # print(fNew)
  }
  
  inc = c - c1;
  while (abs(fNew) > 1e-8)
  {
    inc = -fNew * inc / (fNew - fOld)
    if(inc>1) inc=0.1
    if(inc < (-1)) inc=-0.1
    if((c + inc) >0) {c = c + inc} else {c=c/2}
    fOld = fNew
    g1 = gamma(1 + 1/c)
    g2 = gamma(1 + 2/c)
    # b = ((-a*g1)/g2) + ((((a/g2)**2)*(g2*g2-g1) + (dm*dm/g2))**0.5);
    b = (dm-a)/g1
    fNew =b*b * (g2 - g1**2) - pdvar
    # print(fOld)
    # print(fNew)
    # print(c)
  }
  abc<-cbind(a,b,c)
  abc
}

Fweibull.v<-function(dm,pdvar){
  if (length(dm)==1){abc<-Fweibull(dm,pdvar)}
  else{
    n=length(dm)
    abc<-data.frame(a=rep(NA,n),
                    b=rep(NA,n),
                    c=rep(NA,n))
    for (i in 1:n) {
      abc[i,]<-Fweibull(dm[i],pdvar[i])
    }
  }
  return(abc)
}


pSVDA_2steps <- function(segIDx,nSample,
                         errData1,errData2,errData3, ###error data of the 3 satellite measurements
                         errorPost, #bias error of posterior distribution to check ...
                         step.modelH,step.modelD,step.modelB,  #regression models could be  specific for each period
                         step.modelBconif,step.modelBbl,
                         step.modelH2,step.modelD2,
                         step.modelB2,  #regression models could be  specific for each period
                         step.modelBconif2,step.modelBbl2,
                         allData=F){
  
  kk <- pSVDA(segIDx,nSample,
              errData1,errData2,
              step.modelHx=step.modelH2,step.modelDx=step.modelD2,
              step.modelBx=step.modelB2,
              step.modelBconifx=step.modelBconif2,
              step.modelBblx=step.modelBbl2)
  xx <-segIDx
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
  errDataX <- list(); errDataX$muFSVda <- errorPost ###this is used 
  errDataX$sigmaFSVda <- matrix(round(kk[66:90],6),5,5)

  kk2 <- pSVDA(xx,nSample,
             errDataX,errData3,
             step.modelHx=step.modelH2,step.modelDx=step.modelD2,
             step.modelBx=step.modelB2,
             step.modelBconifx=step.modelBconif2,
             step.modelBblx=step.modelBbl2)
  if(allData){ #return the means and the varCov matrices in vector form
    return(c(segIDx$segID,kk,kk2))
  }else{
    #return just the means and the standard deviations  
    return(c(segIDx$segID,
             kk[1:5],#mean priort t2
             kk[31:35],#mean likel t2
             kk[61:65], #mean postt2
             kk2[1:5],#mean prior t3
             kk2[31:35],#mean likel t3
             kk2[61:65],#mean post t3
             sqrt(kk[seq(6,30,by=6)]),#sd priort t2
             sqrt(kk[seq(36,60,by=6)]),#sd likel t2
             sqrt(kk[seq(66,90,by=6)]), #sd postt2
             sqrt(kk2[seq(6,30,by=6)]),#sd priort t3
             sqrt(kk2[seq(36,60,by=6)]),#sd likel t3
             sqrt(kk2[seq(66,90,by=6)]) #sd postt3
    ))
  } 
}
