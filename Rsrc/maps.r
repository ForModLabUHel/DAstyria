library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(raster)
library(ggpubr)
library(ggplot2)
library(devtools)
library(data.table)
library(ggridges)
library(parallel)

devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/utilStuff/master/ErrorDecomposition/ErrorDecomposition.R")
# source("runSettings.r")
print("define tileX. example: tileX <- '35VLJ'")
tileXs <- c("35VLJ", "34VEQ", "35WMN")
pathX <- "~/research/assessCarbon/results/"
# CSCrun=TRUE
# if(CSCrun==TRUE){
#   pathX <- "/scratch/project_2000994/PREBASruns/assessCarbon/"
# }
coresN <- 3
MSEall <- data.table()
stAll <- data.table()
pRMSE <- list()
nSample <- 100000
colX <- c("#0E95A5","#28B209","#DFB021","#ff8533")

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

# 
# # mclapply(1:length(tileXs),function(i){
for(i in 1:length(tileXs)){
  tileX <- tileXs[i]
  pathLap <- paste0("C:/Users/minunno/Documents/research/assessCarbon/results/",tileX,"/rasters/")
  pathCSC <- paste0("/scratch/project_2000994/PREBASruns/assessCarbon/rasters/Finland/AC_training_FI_",tileX,"/outRast/init2016/")
  pathX <- pathCSC

  Dda <- raster(paste0(pathX,"Dda2019.tif"))
  Dm <- raster(paste0(pathX,"Dm2019.tif"))
  Ds <- raster(paste0(pathX,"Ds2019.tif"))
  mapD <- createMaps(Dda,Ds,Dm,"D","(cm)")
  save(mapD,file = paste0(pathX,"/mapsD.rdata"))
rm(Dda,Dm,Ds,mapD);gc()
  print("mapD")
  Bda <- raster(paste0(pathX,"Bda2019.tif"))
  Bm <- raster(paste0(pathX,"Bm2019.tif"))
  Bs <- raster(paste0(pathX,"Bs2019.tif"))
  mapB <- createMaps(Bda,Bs,Bm,"B","(m2/ha)")
  save(mapB,file = paste0(pathX,"/mapsB.rdata"))
  print("mapB")
rm(Bda,Bm,Bs,mapB);gc()
  Hda <- raster(paste0(pathX,"Hda2019.tif"))
  Hm <- raster(paste0(pathX,"Hm2019.tif"))
  Hs <- raster(paste0(pathX,"Hs2019.tif"))
  mapH <- createMaps(Hda,Hs,Hm,"H","(m)")
  save(mapH,file = paste0(pathX,"/mapsH.rdata"))
  print("mapH")
rm(Hda,Hm,Hs,mapH);gc()

  # ggsave(mapH,filename = paste0(pathX,"/mapH.png"),device = "png")
  # ggsave(mapB,filename = paste0(pathX,"/mapB.png"),device = "png")
  print(tileX)
}
# },mc.cores = coresN)
# 
fileNames <- c("mapsD","mapsH","mapsB")
mapXs <- c("mapD","mapH","mapB")
# mclapply(1:3,function(jx){
  for(jx in 1:3){
    
  
mapX <- list()
for(i in 1:length(tileXs)){
  tileX <- tileXs[i]
  pathLap <- paste0("C:/Users/checcomi/Documents/research/assessCarbon/results/",tileX,"/rasters/")
  pathCSC <- paste0("/scratch/project_2000994/PREBASruns/assessCarbon/rasters/Finland/AC_training_FI_",tileX,"/outRast/init2016/")
  pathX <- pathCSC
  load(paste0(pathX,fileNames[jx],".rdata"))
  mapX[[tileX]] <- get(mapXs[jx])
}
ggsave(get(mapXs[jx])$map,filename = 
         paste0(pathX,fileNames[jx],".pdf"),device = "pdf")
ggsave(get(mapXs[jx])$hist,filename = 
         paste0(pathX,fileNames[jx],"_hist.pdf"),device = "pdf")

plotX <- ggarrange(mapX[[tileXs[1]]]$map,
          mapX[[tileXs[2]]]$map,
          mapX[[tileXs[3]]]$map,
          mapX[[tileXs[1]]]$hist,
          mapX[[tileXs[2]]]$hist,
          mapX[[tileXs[3]]]$hist,
          ncol=3,nrow=2
          )
# save(plotX,file = paste0("/scratch/project_2000994/PREBASruns/assessCarbon/",mapXs[jx],".rdata"))
 ggsave(plotX,filename = paste0("/scratch/project_2000994/PREBASruns/assessCarbon/",mapXs[jx],".pdf"),device = "pdf")
# },mc.cores = coresN)
  }




r <- raster(system.file("external/test.grd", package="raster"))

s1 <- stack(r, r*2)
names(s2) <- c('meuse', 'meuse x 2')

library(ggplot2)

figX <- list()
maps <- list()
dev <- list()
varXs <- c("D","H","B")
unitsX <- c("(cm)","(m)","(m2/ha)")
maxpixels <- 1000000
histX <- list()
pMap <-
# ## With raster
# for(i in 1:3){
#   maps[[i]] <- raster(paste0(pathX,varXs[i],"da2019.tif"))
#   dev$m2019 <- raster(paste0(pathX,varXs[i],"m2019.tif")) - raster(paste0(pathX,varXs[i],"da2019.tif"))
#   dev$s2019 <- raster(paste0(pathX,varXs[i],"s2019.tif")) - raster(paste0(pathX,varXs[i],"da2019.tif"))
# 
#   f1 <- getValues(dev$m2019)
#   f2 <- getValues(dev$s2019)
#   dat1 <- data.table(counts= f1,difX="f1")
#   dat2 <- data.table(counts= f2,difX="f2")
#   dat <- rbind(dat1,dat2)
#   dat$difX <- as.factor(dat$difX)
# 
#   histX[[i]] <- ggplot(dat[sample(1:nrow(dat),1e5)], aes(x=counts, fill=difX)) + 
#    geom_histogram() + ylab(NULL) + xlab(paste0(varXs[i]," offSet ",unitsX[i]))
# 
#   # maps[[i]] <- gplot(ciao,maxpixels=1000000) + geom_tile(aes(fill = value)) +
#   #   facet_wrap(~ varXs[i]) +
#   #   scale_fill_gradient(low = 'white', high = 'blue') +
#   #   coord_equal()
#   print(i)
# }
# 
# 
# maps <- list()
# dev <- list()
# varXs <- c("D","H","B")
# unitsX <- c("cm","m","m2/ha")
for(jx in 1:3){ #loop variables
  for(i in 1:length(tileXs)){ #loop tiles
    tileX <- tileXs[i]
    pathLap <- paste0("C:/Users/checcomi/Documents/research/assessCarbon/results/",tileX,"/rasters/")
    pathCSC <- paste0("/scratch/project_2000994/PREBASruns/assessCarbon/rasters/Finland/AC_training_FI_",tileX,"/outRast/init2016/")
    pathX <- pathCSC
    maps[[i]] <- raster(paste0(pathX,varXs[jx],"da2019.tif"))
    dev$m2019 <- raster(paste0(pathX,varXs[jx],"m2019.tif")) - raster(paste0(pathX,varXs[jx],"da2019.tif"))
    dev$s2019 <- raster(paste0(pathX,varXs[jx],"s2019.tif")) - raster(paste0(pathX,varXs[jx],"da2019.tif"))
    
    f1 <- getValues(dev$m2019)
    f2 <- getValues(dev$s2019)
    dat1 <- data.table(counts= f1,difX="f1")
    dat2 <- data.table(counts= f2,difX="f2")
    dat <- rbind(dat1,dat2)
    dat$difX <- as.factor(dat$difX)
    
    histX[[i]] <- ggplot(dat[sample(1:nrow(dat),1e5)], aes(x=counts, fill=difX)) + 
      geom_histogram(aes(y=..density..)) + ylab(NULL) + xlab(paste0(varXs[jx]," offSet ",unitsX[jx])) +
      theme(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             legend.title = element_blank())+
      scale_fill_discrete(labels = c("m2019-DA2019", "s2019-DA2019"))
    print(tileX)
    pMap[[i]] <- gplot(maps[[i]],maxpixels=maxpixels) + geom_tile(aes(fill = value)) +
      # facet_wrap(~ names) +
      # annotate(geom = 'text', x=1, y=1, label=tileX) + theme_void() +
      scale_fill_gradient(low = 'white', high = 'dark green',name = paste(varXs[jx],unitsX[jx])) +
      scale_alpha(range = c(0.15, 0.65), guide = "none") +
      coord_equal()+  theme(axis.title=element_blank(),
                            axis.text=element_blank(),
                            axis.ticks=element_blank()) +
      ggtitle(tileX) +
      theme(legend.position="top") + 
      ylab(NULL) + xlab(NULL)
  } #end loop tiles
  
  figX[[jx]] <- ggarrange(ggarrange(pMap[[1]], pMap[[2]],pMap[[3]],
                    ncol = 3, common.legend=T,label.x=tileXs), # first row with box and dot plots
                    ggarrange(histX[[1]], histX[[2]],histX[[3]],
                              ncol = 3, common.legend=T), # Second row with box and dot plots
                    nrow = 2) 
  
} #end loop variables
setwd("/scratch/project_2000994/PREBASruns/assessCarbon/")
save(figX,file = "maps.rdata")

####on laptop
setwd("~/research/assessCarbon/")
load("maps.rdata")
ggsave(figX[[1]],filename="D.png",device="png")
ggsave(figX[[2]],filename="H.png",device="png")
ggsave(figX[[3]],filename="B.png",device="png")
# 
mapAll <- ggarrange(figX[[1]],figX[[2]],figX[[3]],
                nrow = 3)
# 
ggsave(mapAll, width = 20, height = 30,device="png",
filename="~/research/assessCarbon/allMaps1.png")
