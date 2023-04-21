library(abind)
library(aaply)
library(plyr)
library(data.table)
library(Rprebasso)

load("C:/Users/minunno/Documents/research/ForestCarbonMonitoring/FCM_CNN/data/Forest_Structure_layers/outRast/uniqueData1.rdata")

startSim=2015
pCROBAS <- pCROB
siteTypeX = 3


ciao <- calcVW(uniqueDataSplit,startSim,pCROBAS)
ciao
