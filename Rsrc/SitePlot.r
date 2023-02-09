library(ggridges)
library(data.table)
library(ggplot2)
library(ggpubr)
load("C:/Users/minunno/Documents/research/assessCarbon/dataForPlots.rdata")


# colX <- c("#66ffd9","#80e5ff","#e6e600","#ff8533")
colX <- c("#0E95A5","#28B209","#e6e600","#ff8533")

vx <- dataAll[varNam=="B" & run=="m2019"]$value - dataAll[varNam=="B" & run=="s2019"]$value
segIDX <- dataAll[varNam=="B" & run=="m2019"]$segID[which(vx < -10 )[1]]
nX <- which(unique(dataAll$segID)==segIDX)
# nX <- 149#123
# segIDX <- unique(dataAll$segID)[nX]

subsetX <- dataAll[segID==segIDX]
subsetX
subsetX$run <- factor(subsetX$run,levels = c("s2016","m2019","s2019","DA2019"))
runX <- "s2016"
varX <- "B"
unitX <- "(m2ha-1)"
sdX <- "vB"
yearX=2016
nSample <- 100000
sampleX <- data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                  sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                      run=runX,year=yearX)

mean(sampleX$value)
yearX = 2019
runX <- "m2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "s2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "DA2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))

# ,""
sampleX$run <- factor(sampleX$run,levels = c("s2016","m2019","s2019","DA2019"))
pB1 <- ggplot(sampleX, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab(paste(varX,unitX)) + xlab("year") + scale_fill_manual(values=colX)

pB



runX <- "s2016"
varX <- "H"
unitX <- "(m)"
sdX <- "vH"
yearX=2016
nSample <- 100000
sampleX <- data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                  sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                      run=runX,year=yearX)

mean(sampleX$value)
yearX = 2019
runX <- "m2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "s2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "DA2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))

sampleX$run <- factor(sampleX$run,levels = c("s2016","m2019","s2019","DA2019"))
pH1 <- ggplot(sampleX, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab(paste(varX,unitX)) + xlab("year")+ scale_fill_manual(values=colX)




runX <- "s2016"
varX <- "D"
unitX <- "(cm)"
sdX <- "vD"
yearX=2016
nSample <- 100000
sampleX <- data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                  sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                      run=runX,year=yearX)

mean(sampleX$value)
yearX = 2019
runX <- "m2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "s2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "DA2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))

sampleX$run <- factor(sampleX$run,levels = c("s2016","m2019","s2019","DA2019"))
pD1 <- ggplot(sampleX, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab(paste(varX,unitX)) + xlab("year")+ scale_fill_manual(values=colX)





load("C:/Users/minunno/Documents/research/assessCarbon/stProbMod3.rdata")


xx <- data.table(year=2019,value=rep(1:5,nSample*stProb[nX,2:6]),
                 run="DA2019")
xx <- rbind(xx,data.table(year=2019,value=rep(1:5,nSample*as.numeric(stProbMod[nX,2:6])),
                          run="m2019"))
xx <- rbind(xx,data.table(year=2016,value=rep(1:5,nSample*probit1[nX,1:5]),
                          run="s2016"))
xx <- rbind(xx,data.table(year=2019,value=rep(1:5,nSample*probit2[nX,1:5]),
                          run="s2019"))

xx$run <- factor(xx$run,levels = c("s2016","m2019","s2019","DA2019"))


pSC <- ggplot(xx, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab("Site class") + xlab("year")+ scale_fill_manual(values=colX)

xy <- data.table(siteClass=1:5,prob=probit1[nX,1:5],run="s2016")
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(stProbMod[nX,2:6]),run="m2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=probit2[nX,1:5],run="s2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=stProb[nX,2:6],run="DA2019"))
xy$run <- factor(xy$run,levels = c("s2016","m2019","s2019","DA2019"))

pSC <- ggplot(data=xy, aes(x=siteClass, y=prob, fill=run)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values=alpha(colX,.3))

# pSite1 <- ggarrange(pB,pH,pD,common.legend = T,nrow=1)
# pSite1

#############site2
nX <- 123
segIDX <- unique(dataAll$segID)[nX]

subsetX <- dataAll[segID==segIDX]
subsetX
subsetX$run <- factor(subsetX$run,levels = c("s2016","m2019","s2019","DA2019"))
runX <- "s2016"
varX <- "B"
unitX <- "(m2ha-1)"
sdX <- "vB"
yearX=2016
nSample <- 100000
sampleX <- data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                  sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                      run=runX,year=yearX)

mean(sampleX$value)
yearX = 2019
runX <- "m2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "s2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "DA2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))

colX <- c("#66ffd9","#80e5ff","#e6e600","#ff8533")
# ,""
sampleX$run <- factor(sampleX$run,levels = c("s2016","m2019","s2019","DA2019"))
pB2 <- ggplot(sampleX, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab(paste(varX,unitX)) + xlab("year") + scale_fill_manual(values=colX)

pB



runX <- "s2016"
varX <- "H"
unitX <- "(m)"
sdX <- "vH"
yearX=2016
nSample <- 100000
sampleX <- data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                  sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                      run=runX,year=yearX)

mean(sampleX$value)
yearX = 2019
runX <- "m2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "s2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "DA2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))

sampleX$run <- factor(sampleX$run,levels = c("s2016","m2019","s2019","DA2019"))
pH2 <- ggplot(sampleX, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab(paste(varX,unitX)) + xlab("year")+ scale_fill_manual(values=colX)




runX <- "s2016"
varX <- "D"
unitX <- "(cm)"
sdX <- "vD"
yearX=2016
nSample <- 100000
sampleX <- data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                  sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                      run=runX,year=yearX)

mean(sampleX$value)
yearX = 2019
runX <- "m2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "s2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))
runX <- "DA2019"
sampleX <- rbind(sampleX,data.table(value=rnorm(nSample,mean=subsetX[run==runX & varNam==varX]$value,
                                                sd=sqrt(subsetX[run==runX & varNam==sdX]$value)),
                                    run=runX,year=yearX))

sampleX$run <- factor(sampleX$run,levels = c("s2016","m2019","s2019","DA2019"))
pD2 <- ggplot(sampleX, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab(paste(varX,unitX)) + xlab("year")+ scale_fill_manual(values=colX)


xx <- data.table(year=2019,value=rep(1:5,nSample*stProb[nX,2:6]),
                 run="DA2019")
xx <- rbind(xx,data.table(year=2019,value=rep(1:5,nSample*as.numeric(stProbMod[nX,2:6])),
                          run="m2019"))
xx <- rbind(xx,data.table(year=2016,value=rep(1:5,nSample*probit1[nX,1:5]),
                          run="s2016"))
xx <- rbind(xx,data.table(year=2019,value=rep(1:5,nSample*probit2[nX,1:5]),
                          run="s2019"))

xx$run <- factor(xx$run,levels = c("s2016","m2019","s2019","DA2019"))


pSC <- ggplot(xx, aes(x=as.factor(year), y=value, width = ..density.., fill=run)) +
  geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.3, scale = 2)+
  ylab("Site class") + xlab("year")+ scale_fill_manual(values=colX)

xy <- data.table(siteClass=1:5,prob=probit1[nX,1:5],run="s2016")
xy <- rbind(xy,data.table(siteClass=1:5,prob=as.numeric(stProbMod[nX,2:6]),run="m2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=probit2[nX,1:5],run="s2019"))
xy <- rbind(xy,data.table(siteClass=1:5,prob=stProb[nX,2:6],run="DA2019"))
xy$run <- factor(xy$run,levels = c("s2016","m2019","s2019","DA2019"))

pSC <- ggplot(data=xy, aes(x=siteClass, y=prob, fill=run)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values=alpha(colX,.3))

# pSite2 <- ggarrange(pB,pH,pD,pSC,common.legend = T)
# pSite2 <- ggarrange(pB,pH,pD,common.legend = T,nrow=1)
# pSite2


 
####plot site type all
# 
# xxAll <-  data.table(year=2019,value=rep(1:5,nSample*colMeans(stProb[,2:6])),
#                      run="DA2019")
# xxAll <- rbind(xxAll,data.table(year=2019,value=rep(1:5,nSample*colMeans(stProbMod[,2:6])),
#                                 run="m2019"))
# xxAll <- rbind(xxAll,data.table(year=2016,value=rep(1:5,nSample*colMeans(probit1[,1:5])),
#                                 run="s2016"))
# xxAll <- rbind(xxAll,data.table(year=2019,value=rep(1:5,nSample*colMeans(probit2[,1:5])),
#                                 run="s2019"))

aa <- 1:1e5
xyAll <- data.table(melt(probit1[aa,]),run="s2016")
ii <- as.matrix(stProbMod[aa,2:6]); ii <- data.table(melt(ii),run="m2019");ii$Var2 <- as.numeric(ii$Var2)
xyAll <- rbind(xyAll,ii)
xyAll <- rbind(xyAll,data.table(melt(probit2[aa,]),run="s2019"))
ii <- as.matrix(stProb[aa,2:6]); ii <- data.table(melt(ii),run="DA2019");ii$Var2 <- as.numeric(ii$Var2)
xyAll <- rbind(xyAll,ii)
xyAll$run <- factor(xyAll$run,levels = c("s2016","m2019","s2019","DA2019"))
setnames(xyAll,"Var2","siteClass")
xyAll$siteClass <- factor(xyAll$siteClass)

pAllST <- ggplot(xyAll, aes(x=siteClass, y=value, fill=run,color=run)) +
  geom_boxplot() + scale_fill_manual(values=alpha(colX,.3)) +
  scale_color_manual(values=colX) + ggtitle("Tile distribution")



####Diffrenvce between s2019 amd m2019
Bx <- data.table(res=dataAll[varNam=="B" & run=="m2019"]$value - dataAll[varNam=="B" & run=="s2019"]$value)
Bhist <-qplot(Bx$res, geom="histogram", main="B")
Hx <- data.table(res=dataAll[varNam=="H" & run=="m2019"]$value - dataAll[varNam=="H" & run=="s2019"]$value)
Hhist <-qplot(Hx$res, geom="histogram", main="H")
Dx <- data.table(res=dataAll[varNam=="D" & run=="m2019"]$value - dataAll[varNam=="D" & run=="s2019"]$value)
Dhist <-qplot(Dx$res, geom="histogram", main="D")
pDev <- ggarrange(Bhist,Hhist,Dhist)


# savePlots
pSite <- ggarrange(pB1,pD1,pH1,pB2,pD2,pH2,nrow=2,ncol=3,common.legend=T)
ggsave("pSite.jpeg",pSite,device = "jpeg")
# ggsave("site2.jpeg",pSite2,device = "jpeg")
ggsave("pDev.jpeg",pDev,device = "jpeg")
ggsave("pAllST.jpeg",pAllST,device = "jpeg")

pSC <- pSC + ggtitle("Pixel calculation (example)")
pST <- plot_grid(pSC, pAllST, labels=c("A", "B"), ncol = 2, nrow = 1)
# pST <- ggarrange(pSC,pAllST)
ggsave("pSC.jpeg",pST,device = "jpeg")
  
