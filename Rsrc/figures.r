library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(raster)
library(ggpubr)
library(ggplot2)
library(devtools)
library(data.table)
library(ggridges)
library(Metrics)
devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/utilStuff/master/ErrorDecomposition/ErrorDecomposition.R")
# source("runSettings.r")
print("define tileX. example: tileX <- '35VLJ'")
tiles <- c("35VLJ", "34VEQ", "35WMN")
sites <- c("Tampere","Vaasa","Oulu")
pathX <- "C:/Users/minunno/Documents/research/assessCarbon/"
# CSCrun=TRUE
# if(CSCrun==TRUE){
#   pathX <- "/scratch/project_2000994/PREBASruns/assessCarbon/"
# }
namesTab <- c(paste0("B",c("m2019","s2019","DA2019")),
              paste0("H",c("m2019","s2019","DA2019")),
              paste0("D",c("m2019","s2019","DA2019")))
tabPrmse <- tabPbias <- matrix(NA,9,3,dimnames = list(namesTab,tiles))
tabRmse <- tabBias <- matrix(NA,9,3,dimnames = list(namesTab,tiles))

p_rmse <- function(obs,sim) sqrt( mean( (obs-sim)^2) ) / ( max(obs)-min(obs) )*100
p_bias <- function(obs,sim) mean( (obs-sim))/(max(obs)-min(obs) )*100
# msex <- function(sim,obs)  mean( (sim-obs)^2)
# rmsex <- function(sim,obs) sqrt( mean( (sim-obs)^2) )
# biasx <- function(sim,obs) mean( (sim-obs))


MSEall <- data.table()
stAll <- data.table()
pRMSE <- list()
nSample <- 100000
colX <- c("#0E95A5","#28B209","#DFB021","#ff8533")

load(paste0(pathX,"results/dataRes2019.rdata"))
for(i in 1:length(tiles)){
  tileX <- tiles[i]
  #####Figure 1 reseults 
  # load(paste0(pathX,"data2019res_",tileX,".rdata"))
  data2019res <- dataRes2019[Tile==tileX]
  plot(data2019res$Bs2019,data2019res$G)
  points(data2019res$Bda2019,data2019res$G,col=2,pch=20)
  points(data2019res$Bm2019,data2019res$G,col=3,pch=20)
  abline(0,1)
  
  nas <- intersect(which(!is.na(data2019res$Bm2019)),which(!is.na(data2019res$G)))
  nas <- intersect(nas,which(!is.na(data2019res$Bs2019)))
  nas <- intersect(nas,which(!is.na(data2019res$Bda2019)))
  tabPrmse[1,i] <- p_rmse(data2019res$G[nas],data2019res$Bm2019[nas])
  tabPrmse[2,i] <- p_rmse(data2019res$G[nas],data2019res$Bs2019[nas])
  tabPrmse[3,i] <- p_rmse(data2019res$G[nas],data2019res$Bda2019[nas])
  tabPbias[1,i] <- p_bias(data2019res$G[nas],data2019res$Bm2019[nas])
  tabPbias[2,i] <- p_bias(data2019res$G[nas],data2019res$Bs2019[nas])
  tabPbias[3,i] <- p_bias(data2019res$G[nas],data2019res$Bda2019[nas])
  tabRmse[1,i] <- rmse(data2019res$G[nas],data2019res$Bm2019[nas])
  tabRmse[2,i] <- rmse(data2019res$G[nas],data2019res$Bs2019[nas])
  tabRmse[3,i] <- rmse(data2019res$G[nas],data2019res$Bda2019[nas])
  tabBias[1,i] <- bias(data2019res$G[nas],data2019res$Bm2019[nas])
  tabBias[2,i] <- bias(data2019res$G[nas],data2019res$Bs2019[nas])
  tabBias[3,i] <- bias(data2019res$G[nas],data2019res$Bda2019[nas])

  nas <- intersect(which(!is.na(data2019res$Hm2019)),which(!is.na(data2019res$H)))
  nas <- intersect(nas,which(!is.na(data2019res$Hs2019)))
  nas <- intersect(nas,which(!is.na(data2019res$Hda2019)))
  tabPrmse[4,i] <- p_rmse(data2019res$H[nas],data2019res$Hm2019[nas])
  tabPrmse[5,i] <- p_rmse(data2019res$H[nas],data2019res$Hs2019[nas])
  tabPrmse[6,i] <- p_rmse(data2019res$H[nas],data2019res$Hda2019[nas])
  tabPbias[4,i] <- p_bias(data2019res$H[nas],data2019res$Hm2019[nas])
  tabPbias[5,i] <- p_bias(data2019res$H[nas],data2019res$Hs2019[nas])
  tabPbias[6,i] <- p_bias(data2019res$H[nas],data2019res$Hda2019[nas])
  tabRmse[4,i] <- rmse(data2019res$H[nas],data2019res$Hm2019[nas])
  tabRmse[5,i] <- rmse(data2019res$H[nas],data2019res$Hs2019[nas])
  tabRmse[6,i] <- rmse(data2019res$H[nas],data2019res$Hda2019[nas])
  tabBias[4,i] <- bias(data2019res$H[nas],data2019res$Hm2019[nas])
  tabBias[5,i] <- bias(data2019res$H[nas],data2019res$Hs2019[nas])
  tabBias[6,i] <- bias(data2019res$H[nas],data2019res$Hda2019[nas])
  
  nas <- intersect(which(!is.na(data2019res$Dm2019)),which(!is.na(data2019res$D)))
  nas <- intersect(nas,which(!is.na(data2019res$Ds2019)))
  nas <- intersect(nas,which(!is.na(data2019res$Dda2019)))
  tabPrmse[7,i] <- p_rmse(data2019res$D[nas],data2019res$Dm2019[nas])
  tabPrmse[8,i] <- p_rmse(data2019res$D[nas],data2019res$Ds2019[nas])
  tabPrmse[9,i] <- p_rmse(data2019res$D[nas],data2019res$Dda2019[nas])
  tabPbias[7,i] <- p_bias(data2019res$D[nas],data2019res$Dm2019[nas])
  tabPbias[8,i] <- p_bias(data2019res$D[nas],data2019res$Ds2019[nas])
  tabPbias[9,i] <- p_bias(data2019res$D[nas],data2019res$Dda2019[nas])
  tabRmse[7,i] <- rmse(data2019res$D[nas],data2019res$Dm2019[nas])
  tabRmse[8,i] <- rmse(data2019res$D[nas],data2019res$Ds2019[nas])
  tabRmse[9,i] <- rmse(data2019res$D[nas],data2019res$Dda2019[nas])
  tabBias[7,i] <- bias(data2019res$D[nas],data2019res$Dm2019[nas])
  tabBias[8,i] <- bias(data2019res$D[nas],data2019res$Ds2019[nas])
  tabBias[9,i] <- bias(data2019res$D[nas],data2019res$Dda2019[nas])
}

write.csv(t(rbind(tabPrmse,tabPbias)),file=paste0(pathX,"tab1_a.csv"))
write.csv(t(rbind(tabRmse,tabBias)),file=paste0(pathX,"tab1_b.csv"))

MSEall <- data.table()
for(i in 1:length(tiles)){
  tileX <- tiles[i]
  data2019res <- dataRes2019[Tile==tileX]
  
  #####Figure 1 reseults 
  # load(paste0(pathX,"data2019res_",tileX,".rdata"))
  
  # MSEs <- data.table(value=as.numeric(unlist(MSEdec("B_est", data2019res$G.est,data2019res$G.mea,method=1))[2:5]),
  #                    run="est",variable="B",components=c("sb","sdsd","lc","mse"),
  #                    rangeObs = max(data2019res$G.mea,na.rm=T)-min(data2019res$G.mea,na.rm=T))
  nas <- intersect(which(!is.na(data2019res$Bm2019)),which(!is.na(data2019res$G)))
  nas <- intersect(nas,which(!is.na(data2019res$Bs2019)))
  nas <- intersect(nas,which(!is.na(data2019res$Bda2019)))
  MSEs <- data.table(value=as.numeric(unlist(MSEdec("B", data2019res$Bm2019[nas],
                                 data2019res$G[nas],method=1))[2:5]),
                                run="m2019",variable="B",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$G,na.rm=T)-min(data2019res$G,na.rm=T))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("B",
                                data2019res$Bs2019[nas],data2019res$G[nas],method=1))[2:5]),
                                run="s2019",variable="B",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$G,na.rm=T)-min(data2019res$G,na.rm=T)))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("B", data2019res$Bda2019[nas],
                                 data2019res$G[nas],method=1))[2:5]),
                                run="DA2019",variable="B",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$G,na.rm=T)-min(data2019res$G,na.rm=T)))

  nas <- intersect(which(!is.na(data2019res$Dm2019)),which(!is.na(data2019res$D)))
  nas <- intersect(nas,which(!is.na(data2019res$Ds2019)))
  nas <- intersect(nas,which(!is.na(data2019res$Dda2019)))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("D", data2019res$Dm2019[nas],data2019res$D[nas],method=1))[2:5]),
                                run="m2019",variable="D",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$D,na.rm=T)-min(data2019res$D,na.rm=T)))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("D", data2019res$Ds2019[nas],data2019res$D[nas],method=1))[2:5]),
                                run="s2019",variable="D",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$D,na.rm=T)-min(data2019res$D,na.rm=T)))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("D", data2019res$Dda2019[nas],data2019res$D[nas],method=1))[2:5]),
                                run="DA2019",variable="D",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$D,na.rm=T)-min(data2019res$D,na.rm=T)))
  
  nas <- intersect(which(!is.na(data2019res$Hm2019)),which(!is.na(data2019res$H)))
  nas <- intersect(nas,which(!is.na(data2019res$Hs2019)))
  nas <- intersect(nas,which(!is.na(data2019res$Hda2019)))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("H", data2019res$Hm2019[nas],data2019res$H[nas],method=1))[2:5]),
                                run="m2019",variable="H",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$H,na.rm=T)-min(data2019res$H,na.rm=T)))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("H", data2019res$Hs2019[nas],data2019res$H[nas],method=1))[2:5]),
                                run="s2019",variable="H",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$H,na.rm=T)-min(data2019res$H,na.rm=T)))
  MSEs <- rbind(MSEs,data.table(value=as.numeric(unlist(MSEdec("H", data2019res$Hda2019[nas],data2019res$H[nas],method=1))[2:5]),
                                run="DA2019",variable="H",components=c("sb","sdsd","lc","mse"),
                                rangeObs = max(data2019res$H,na.rm=T)-min(data2019res$H,na.rm=T)))
  
  MSEs[components=="mse", RMSE:=sqrt(value)]
  
  # for(runX in unique(MSEs$run)){
  #   for(vX in unique(MSEs$variable)){
  #     MSEs[run==runX & variable==vX & !components%in%"mse",RMSE:=value/sum(value)]
  #     MSEs[run==runX & variable==vX & !components%in%"mse",]$RMSE <- MSEs[run==runX & variable==vX & !components%in%"mse",]$RMSE * MSEs[run==runX & variable==vX & components%in%"mse",]$RMSE
  #   }
  # }
  
  # MSEs$pRMSE <- MSEs$RMSE/MSEs$rangeObs*100
  
  
  MSEs$run <- factor(MSEs$run, levels = c("m2019","s2019","est","DA2019"))
  MSEs$components <- factor(MSEs$components, levels = c("sb","sdsd","lc","mse"))
  varX = "B"
  MSEs$Tile <- tileX
  MSEs$site <- sites[i]
    
  MSEall <- rbind(MSEall,MSEs)

  
  
  #####start proc data Figure 2 results 
  load(paste0(paste0(pathX,"/results/",tileX,"/stProbMod1.rdata")))

  xxAll <-  data.table(year=2019,value=rep(1:5,nSample*colMeans(stProb[,2:6])),
                       run="DA2019")
  xxAll <- rbind(xxAll,data.table(year=2019,value=rep(1:5,nSample*colMeans(stProbMod[,2:6])),
                                  run="m2019"))
  xxAll <- rbind(xxAll,data.table(year=2016,value=rep(1:5,nSample*colMeans(probit1[,1:5])),
                                  run="s2016"))
  xxAll <- rbind(xxAll,data.table(year=2019,value=rep(1:5,nSample*colMeans(probit2[,1:5])),
                                  run="s2019"))

  probit1 <- data.table(probit1)
  aa <- 1:1e5
  xyAll <- data.table(melt(probit1[aa,],variable.name = "siteClass",measure.vars = 1:5)
                      ,run="s2016")
  setnames(stProbMod,c("segID",1:5))
  ii <- data.table(melt(stProbMod[aa,2:6],measure.vars = 1:5,variable.name = "siteClass"),
             run="m2019")
  xyAll <- rbind(xyAll,ii)
  probit2 <- data.table(probit2)
  ii <- data.table(melt(probit2[aa,],variable.name = "siteClass",measure.vars = 1:5)
             ,run="s2019")
  xyAll <- rbind(xyAll,ii)
  stProb <- data.table(stProb)
  setnames(stProb,c("segID",1:5))
  ii <- data.table(melt(stProb[aa,2:6],measure.vars = 1:5,variable.name = "siteClass"),
                   run="DA2019")
  xyAll <- rbind(xyAll,ii)
  xyAll$run <- factor(xyAll$run,levels = c("s2016","m2019","s2019","DA2019"))
  xyAll$siteClass <- factor(xyAll$siteClass)

  ciao <- xyAll[,median(value),by=.(siteClass,run)]

  oo <- ciao[,spline(siteClass,V1),by=run]

  # ggplot(oo) +
  #   geom_line(aes(x = x, y = y, colour = run)) +
  #   coord_cartesian(ylim = c(0, 1.)) + theme(legend.title=element_blank())
  oo$tile <- tileX
  ciao$tile <- tileX
  xyAll$tile <- tileX
  xyAll$site <- sites[i] 
  
  stAll <- rbind(stAll,xyAll)
}
MSEall[variable=="B"]$variable <- "G"


# pRMSE <- ggplot(MSEall[!components %in% "mse"& !run %in% "est"], 
#                 aes(x=run, y=pRMSE, fill=components)) +
#   geom_bar(stat="identity")#, position = position_dodge(0.9)) +
# # scale_fill_manual(values = c("#00AFBB", "#E7B800"))
# RMSE <- ggplot(MSEall[!components %in% "mse"& !run %in% "est"], 
#                aes(x=run, y=RMSE, fill=components)) +
#   geom_bar(stat="identity")#, position = position_dodge(0.9)) +
# p + facet_grid(rows = vars(variable))
MSE <- ggplot(MSEall[!components %in% "mse" &
                       !run %in% "est"], 
               aes(x=run, y=value, fill=components)) +
  geom_bar(stat="identity")#, position = position_dodge(0.9)) +
# figRes1a <- pRMSE + facet_grid(rows = vars(variable), cols = vars(Tile)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# figRes1b <- RMSE + facet_grid(rows = vars(variable), cols = vars(Tile)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
MSEall$site <- factor(MSEall$site,levels = c("Tampere","Vaasa","Oulu"))

figRes1c <- MSE + facet_grid(rows = vars(variable), cols = vars(site),scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
figRes1c

# resFig1 <- figRes1c#ggarrange(figRes1a,figRes1b)
# ggsave(resFig1,filename = paste0(pathX,"/figures/resFig1.png"),device = "png")
# ggsave(figRes1a,filename = paste0(pathX,"/figures/resFig1pRMSE.png"),device = "png")
# ggsave(figRes1b,filename = paste0(pathX,"figures/resFig1rmse"),device = "png")
ggsave(figRes1c,filename = paste0(pathX,"/results/figures/resFig1mse.png"),device = "png",height=5,width=5)

# ggsave(figRes1b,filename = "resFig1rmse",device = "png")

# stPlot <- ggplot(stAll) + 
#   geom_line(aes(x = as.factor(siteClass), y = V1, colour = as.factor(run),group=as.factor(run))) +
#   # coord_cartesian(ylim = c(-0.10, 1.)) +
#   theme(legend.title=element_blank())
stAll$site <- factor(stAll$site,levels = c("Tampere","Vaasa","Oulu"))
stAll$run <- factor(stAll$run,levels = c("s2016","m2019","s2019","DA2019"))

stPlot <- ggplot(stAll, aes(x=siteClass, y=value)) +
    geom_boxplot(outlier.shape ="-") + ylab(NULL) + #scale_fill_manual(values=alpha(colX,0.3)) +
    #scale_color_manual(values=colX)+ 
  theme(legend.title=element_blank())

####Figure 2
figRes2 <- stPlot + facet_grid( cols = vars(site),rows = vars(run)) 
figRes2

ggsave(figRes2,filename = paste0(pathX,"/results/figures/resFig2.png"),
       device = "png",width = 7,h=7)

xpX <- dataPlot <- list()

####Figure 2.1
nX <- 123
xy <- data.table(siteClass=1:5,prob=as.numeric(probit1[nX,1:5]),run="s2016")
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(stProbMod[nX,2:6]),run="m2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(probit2[nX,1:5]),run="s2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(stProb[nX,2:6]),run="DA2019"))
xy$run <- factor(xy$run,levels = c("s2016","m2019","s2019","DA2019"))

pSC1 <- ggplot(data=xy, aes(x=siteClass, y=prob, color=run,group=run)) +
  geom_point() + geom_line() + theme(legend.title=element_blank())+
  scale_color_manual(values=colX) + ylab(NULL) + ggtitle("A")
pSC1

nX <- 12333
xy <- data.table(siteClass=1:5,prob=as.numeric(probit1[nX,1:5]),run="s2016")
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(stProbMod[nX,2:6]),run="m2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(probit2[nX,1:5]),run="s2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(stProb[nX,2:6]),run="DA2019"))
xy$run <- factor(xy$run,levels = c("s2016","m2019","s2019","DA2019"))

pSC2 <- ggplot(data=xy, aes(x=siteClass, y=prob, color=run,group=run)) +
  geom_point() + geom_line() + theme(legend.title=element_blank()) +
  scale_color_manual(values=colX) + ylab(NULL) + ggtitle("B")
pSC2
figRes2.1 <- ggarrange(pSC1,pSC2,common.legend = T)

ggsave(figRes2.1,filename = paste0(pathX,"/results/figures/resFig2.1.png"),
       device = "png",width = 7,height=7)






dataPlot <- data.table()

####Figure 3
for(i in 1:length(tiles)){
  tileX <- tiles[i]
  load(paste0(paste0(pathX,"/results/",tileX,"/dataForPlots.rdata")))
  
  dataAll$tile <- tileX
  dataPlot <- rbind(dataPlot,dataAll)
}
  
  vars <- unique(dataAll$varNam)
  units <- c("H (m)","D (cm)","BA (m2)", "cover (%)", "cover (%)", "cover (%)",
             "variance (m)","variance (cm)","variance (m2)",
             "variance (%)", "variance (%)", "variance (%)")
  
  ciao <- c("H", "B", "D","pine (%)","spruce (%)","birch (%)",
            "varH", "varD", "varB",
            "variance pine (%)","variance spruce (%)","variance birch (%)")
  
  my_labeller = as_labeller(
    c(
      H = 'H~(m)', 
      B = 'G~(m^2*ha^-1)', 
      `D` = 'D~(cm)'
    ), 
    default = label_parsed
  )
  my_labeller_var = as_labeller(
    c(
      varD = 'D~variance', 
      varH = 'H~variance', 
      varB = 'G~variance'
    ), 
    default = label_parsed
  )
  dataPlot$titles <- ciao[match(dataPlot$varNam,unique(dataPlot$varNam))]
  dataPlot$site <- sites[match(dataPlot$tile,tiles)]
  dataPlot$site <- factor(dataPlot$site,levels=c("Tampere","Vaasa","Oulu"))
  
  pFSV <- ggplot(dataPlot[varNam%in%vars[1:3]],
                                aes(x = value, y = run, fill = stat(x))) + 
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "", option = "C") +
     xlab(NULL)+ ylab(NULL)
  pvFSV <- ggplot(dataPlot[varNam%in%vars[7:9]],
                 aes(x = value, y = run, fill = stat(x))) + 
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "", option = "C") +
    xlab(NULL)+ ylab(NULL)
  
  pCover <- ggplot(dataPlot[varNam%in%vars[4:6]],
                 aes(x = value, y = run, fill = stat(x))) + 
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "", option = "C") +
    xlab(NULL)+ ylab(NULL)
  pvCover <- ggplot(dataPlot[varNam%in%vars[10:12]],
                   aes(x = value, y = run, fill = stat(x))) + 
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "", option = "C") +
    xlab(NULL)+ ylab(NULL)
  
    
    figRes3a <- pFSV + facet_grid(rows = vars(site), cols = vars(titles),labeller = my_labeller) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    figRes3b <- pCover + facet_grid(rows = vars(site), cols = vars(titles)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    figRes3c <- pvFSV + facet_grid(rows = vars(site), cols = vars(titles),labeller = my_labeller_var) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    figRes3d <- pvCover + facet_grid(rows = vars(site), cols = vars(titles)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    figRes3a
    figRes3b
    figRes3c
    figRes3d
    
    ggsave(figRes3a,filename = paste0(pathX,"results/figures/resFig3a.png"),device = "png",width = 7,h=5)
    ggsave(figRes3b,filename = paste0(pathX,"results/figures/resFig3b.png"),device = "png",width = 7,h=5)
    ggsave(figRes3c,filename = paste0(pathX,"results/figures/resFig3c.png"),device = "png",width = 7,h=5)
    ggsave(figRes3d,filename = paste0(pathX,"results/figures/resFig3d.png"),device = "png",width = 7,h=5)
    

##Maps and hist
    ####create Maps
    createMaps <- function(rastDA,rastS,rastM,varX,unitsX){
      # rast_df <- as.data.frame(rastDA, xy = TRUE)
      rast_df <- data.table(rasterToPoints(rastDA))
      fillX <- names(rast_df)[3]
      mapX <- ggplot() +
        geom_raster(data = rast_df, 
                    aes_string(x = "x", y = "y", 
                        fill = fillX)) + 
        geom_raster(data = rast_df, 
                    aes_string(x = "x", y = "y", 
                        alpha = fillX)) +  
        scale_fill_viridis_c() +  
        scale_alpha(range = c(0.15, 0.65), guide = "none") +  
        ggtitle(varX) + labs(fill = paste(varX, unitsX))
        coord_quickmap()
      
      f1 <- getValues(rastDA - rastS)
      f2 <- getValues(rastDA - rastM)
      dat1 <- data.table(counts= f1,difX="f1")
      dat2 <- data.table(counts= f2,difX="f2")
      dat <- rbind(dat1,dat2)
      dat$difX <- as.factor(dat$difX)
      
      histX <- ggplot(dat[sample(1:nrow(dat),1e5)], aes(x=counts, fill=difX)) + 
        geom_histogram() + ylab(NULL) + xlab(paste0(varX," deviation ",unitsX))
      plotX <- list(map = mapX,hist=histX)
      return(plotX)
    }
    
for(i in 1:length(tiles)){
  tileX <- tiles[i]
  pathLap <- paste0(paste0(pathX,"results/",tileX,"/rasters/"))
  pathCSC <- paste0("/scratch/project_2000994/PREBASruns/assessCarbon/rasters/Finland/AC_training_FI_",tileX,"/outRast/init2016/")
  pathX <- pathCSC
  
  Dda <- raster(paste0(pathX,"Dda2019.tif"))
  Dm <- raster(paste0(pathX,"Dm2019.tif"))
  Ds <- raster(paste0(pathX,"Ds2019.tif"))
  Bda <- raster(paste0(pathX,"Bda2019.tif"))
  Bm <- raster(paste0(pathX,"Bm2019.tif"))
  Bs <- raster(paste0(pathX,"Bs2019.tif"))
  Hda <- raster(paste0(pathX,"Hda2019.tif"))
  Hm <- raster(paste0(pathX,"Hm2019.tif"))
  Hs <- raster(paste0(pathX,"Hs2019.tif"))
  
  mapD <- createMaps(Dda,Ds,Dm,"D","(cm)")
  mapH <- createMaps(Dda,Ds,Dm,"H","(m)")
  mapB <- createMaps(Dda,Ds,Dm,"B","(m2/ha)")

  save(mapD,file = paste0(pathX,"results/mapsD.rdata"))
  save(mapH,file = paste0(pathX,"results/mapsH.rdata"))
  save(mapB,file = paste0(pathX,"results/mapsB.rdata"))
  # ggsave(mapH,filename = paste0(pathX,"/mapH.png"),device = "png")
  # ggsave(mapB,filename = paste0(pathX,"/mapB.png"),device = "png")
  print(tileX)
}

mapX <- list()
for(i in 1:length(tiles)){
  tileX <- tiles[i]
  pathLap <- paste0(pathX,"results/",tileX,"/rasters/")
  pathCSC <- paste0("/scratch/project_2000994/PREBASruns/assessCarbon/rasters/Finland/AC_training_FI_",tileX,"/outRast/init2016/")
  pathX <- pathCSC
  load(paste0(pathX,"/mapsD.rdata"))
  mapX[[tileX]] <- mapD
}
plotX <- ggarrange(mapX[[tiles[1]]]$map,
          mapX[[tiles[2]]]$map,
          mapX[[tiles[3]]]$map,
          mapX[[tiles[1]]]$hist,
          mapX[[tiles[2]]]$hist,
          mapX[[tiles[3]]]$hist,ncol=3,nrow=3
          ) 
 ggsave(Dx,filename = "/scratch/project_2000994/PREBASruns/assessCarbon/Dmaps.png",device = "png")

 