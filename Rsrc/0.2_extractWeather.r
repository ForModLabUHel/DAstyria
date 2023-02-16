library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/utilStuff/master/Clipick/extractClipickData.r")
source_url("https://raw.githubusercontent.com/ForModLabUHel/DAstyria/master/Rsrc/settings.r")

setwd(generalPath)

shapeRastOr <- raster("2015/FCM_STY_2015_G_10M_8BITS-Styria.tif")
shapeRast <- aggregate(shapeRastOr,fact=500)
projection(shapeRast)

myRaster <- projectRaster(shapeRast,crs =  "+proj=longlat +datum=WGS84")
plot(myRaster)

# extract the Altri weather data from clipick from 1.1.2000-2.1.2000 at 10 km grid
wDs <- getWD("Styria", myRaster, 10, 2014, 1, 1, 2022, 12, 31)

# Find the unique data tables
# weather data for the site is here 
weather_data <- unique(wDs)

length(wDs)
length(weather_data)


# raster of climID:s is here
climID_raster <- climRaster(weather_data, wDs, "Styria", myRaster, 10000)

###Converts in the orginal coordinate system of the input raster or shape File
climID_raster <- projectRaster(climID_raster,
                                       crs = crs(shapeRastOr),
                                       res = res(shapeRastOr),
                                method = "ngb") 

climID_raster <- crop(climID_raster,extent(shapeRastOr))
climID_raster <- raster(vals=values(climID_raster),ext=extent(shapeRastOr),crs=crs(shapeRastOr),
              nrows=dim(shapeRastOr)[1],ncols=dim(shapeRastOr)[2])

climID_raster <- mask(climID_raster,shapeRastOr)

writeRaster(climID_raster,filename = "weatherInputs/climID.tif",overwrite=T)

# # weather data as one data.table, climID as "id"
# weather_data_table <- rbindlist(weather_data, idcol="id")

prebasWeather(weather_data)

