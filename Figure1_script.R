#PLEASE REFER TO README.txt TO KNOW REQUIRED INPUTS

####MUST BE DEFINED BY USER #####
session_Fig1<-"your_session_directory_here"
session_Fig1<-"E:/BARCODE_MARONI/SCRIPTS_AND_DATA/FIGURE1"
setwd(session_Fig1)

library(raster)
#extent, raster, plot, shapefile
library(rgdal)
#readOGR
library(shapefiles)
#read.shapefile
library(mapplots)
#draw.shape

#plot elevation basemap
ext_basemap<-extent(-56.2,-53.0, 2.1, 5.85)
pal <- colorRampPalette(c("yellowgreen","chartreuse4","goldenrod","brown"))
basemap<-raster("shapefiles/sa_dem_30s.bil")   #basemap
plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3")

#plot basins limits after cropping for our extant
basins<-shapefile("shapefiles/sa_bas_30s_beta.shp")  #~1min
sub <- crop(basins,ext_basemap)
shapefile(sub, 'shapefiles/maroni_basins.shp',overwrite=TRUE) #maroni_basins.shp is a cropped version of sa_bas_30s
basins<- readOGR("shapefiles","maroni_basins")
maroni_basins<- subset(basins, AREA_SQKM>5) #We remove all the micro-basins on the coast, they are artefacts from the sat measures
plot(maroni_basins,add=T)

#REMINDER: If you want to access the Maroni in the sahapefile basin:  BASIN_ID[233]ou BASIN_ID==16680
#To access the small basin piece at the mouth:BASIN_ID[154]ou BASIN_ID==16555

#plot the rivers, again after cropping 
rivs<-shapefile("shapefiles/sarivs.shp") # takes a while, heavy file (~5min)
sub2 <- crop(rivs,ext_basemap)
shapefile(sub2, 'shapefiles/maroni_rivs.shp')
maroni_rivs<-read.shapefile('shapefiles/maroni_rivs')
draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)

#plot additional Maroni rivers and some lakes of Suriname
sur_waters<-shapefile("shapefiles/SUR_water_areas_dcw") 
sub3 <- crop(sur_waters,ext_basemap)
shapefile(sub3, 'shapefiles/maroni_sur_waters.shp')
maroni_sur_waters<-read.shapefile('shapefiles/maroni_sur_waters')
draw.shape(maroni_sur_waters, type = "lines", col = "dodgerblue3", lwd=0.8)

#plot a few more rivers from French Guiana
gf_waters<-shapefile("shapefiles/SUR_water_areas_dcw") 
sub4 <- crop(gf_waters,ext_basemap)
shapefile(sub4, 'shapefiles/maroni_gf_waters.shp')
maroni_gf_waters<-read.shapefile('shapefiles/maroni_gf_waters')
draw.shape(maroni_gf_waters, type = "lines", col = "dodgerblue3", lwd=0.8)

coordcam<-read.table("coord.txt", header=T, row.names=NULL,sep="\t")
points(coordcam$Lon,coordcam$Lat, pch=21, cex = coordcam$Size,col="black",bg="red")

save.image("Figure1_all_data.RData")

#The final manuscript version has been modified on Inkscape
rm(basins,coordcam,gf_waters,maroni_gf_waters,maroni_sur_waters,rivs)
rm(sub,sub2,sub3,sub4,sur_waters)

#objects needed for landscape pipeline: basemap, ext_basemap, pal, maroni_basins, maroni_rivs

save.image("Figure1_data.RData")

#### package summary and unload ####

library(NCmisc)
list.functions.in.file("Figure1_script.R")

detach("package:raster", unload=TRUE)
detach("package:rgdal", unload=TRUE)
detach("package:shapefiles", unload=TRUE)
detach("package:mapplots", unload=TRUE)
detach("package:sp", unload=TRUE)
detach("package:NCmisc", unload=TRUE)
(.packages())
