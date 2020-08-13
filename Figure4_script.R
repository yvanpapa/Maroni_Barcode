####MUST BE DEFINED BY USER #####
#Inputs needed: 
#Figure 1 environment
#Session 2 environment
###Maroni shapefile env is not enough anymore (raster update). Now you also need:
#shapefiles directory with sa_dem_30s.bil

library(rasterVis)
library(rgdal)
#readOGR
library(grid)
#grid.text

session2<-"your_session2_directory_here"
load(paste0(session2,"/env_landscapes_all_data.RData")) #~10 sec

session_Fig1<-"your_Fig1_directory_here"

session_Fig4<-"your_session_directory_here"
setwd(session_Fig4)

####
basemap<-raster("shapefiles/sa_dem_30s.bil")
crop(basemap,ext_basemap)->basemap2  
readOGR(paste0(session_Fig1,"/shapefiles"),"maroni_rivs")->maroni_rivs2
levelplot(basemap2,margin=F,colorkey=F,region=F)->rast1 #rast1= croppped south america
layer(sp.polygons(maroni_basins))+layer(sp.polygons(maroni_rivs2, col = "darkgrey", lwd=1))->lay1   #basins limits and rivers
extend(rnorm_xyz_total.masked,ext_basemap)->rnorm_xyz_total.masked2
extend(count_total.masked,ext_basemap)->count_total.masked2

dev.off()
dev.new()
#plot.new()
levelplot(rnorm_xyz_total.masked2,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label="All species"),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(1,1,2,1),newpage=F)

levelplot(count_total.masked2,col.regions=cols_tr,at=seq(0, 100, length.out=11), colorkey = list(space='right'),margin=F,main=list(label="Number of species computed by pixel"),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(2,1,2,1),newpage=F)
grid.text(x=.08,y=.91,label="A",gp=gpar(fontsize=22, col="BLACK",fontface="bold"))
grid.text(x=.55,y=.91,label="B",gp=gpar(fontsize=22, col="BLACK",fontface="bold"))
dev.print(pdf, paste0('Figure4.pdf')) # takes ~30sec to open with Adobe
dev.print(svg, paste0('Figure4.svg')) # takes ~1 min to open with inkscape
dev.print(png, paste0('Figure4.png'),width = 1024, height = 768)
dev.off()



####ALTERNATIVE COLOR GRADIENT VERSION
dev.off()
dev.new()
#plot.new()
levelplot(rnorm_xyz_total.masked2,col.regions=cols_gr,at=seq(0, 1, length.out=1001), colorkey = list(space='right'),margin=F,main=list(label="All species",cex=1.5),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(1,1,2,1),newpage=F)

levelplot(count_total.masked2,col.regions=cols_gr,at=seq(0, 100, length.out=1001), colorkey = list(space='right'),margin=F,main=list(label="Number of species computed by pixel"),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(2,1,2,1),newpage=F)
grid.text(x=.08,y=.91,label="A",gp=gpar(fontsize=22, col="BLACK",fontface="bold"))
grid.text(x=.55,y=.91,label="B",gp=gpar(fontsize=22, col="BLACK",fontface="bold"))
dev.print(pdf, paste0('Figure4_gr.pdf')) # takes ~30sec to open with Adobe
dev.print(svg, paste0('Figure4_gr.svg')) # takes ~1 min to open with inkscape
dev.print(png, paste0('Figure4_gr.png'),width = 1024, height = 768)
dev.off()

save.image("Figure4_data.RData")

######## LIGHT VERSION ##########################

dev.off()
dev.new()
#plot.new()
levelplot(rnorm_xyz_total.masked2,col.regions=cols_gr,at=seq(0, 1, length.out=1001), colorkey = list(space='right'),
          margin=F,main=list(label="Mean \"multispecies\" genetic landscape",cex=1.4),scales=list(draw=T,cex=1.2))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot)
dev.print(pdf, paste0('Figure4A_gr_light.pdf')) # takes ~30sec to open with Adobe
dev.print(svg, paste0('Figure4A_gr_light.svg')) # takes ~1 min to open with inkscape
dev.print(png, paste0('Figure4A_gr_light.png'),width = 1024, height = 768)
dev.off()

dev.off()
dev.new()
levelplot(count_total.masked2,col.regions=cols_tr,at=seq(0, 100, length.out=11), colorkey = list(space='right',labels=list(cex=1.8)),
          margin=F,main=list(label="Species landscapes coverage"),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
#rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
print(rast2)
#grid.text(x=.08,y=.91,label="A",gp=gpar(fontsize=22, col="BLACK",fontface="bold"))
#grid.text(x=.55,y=.91,label="B",gp=gpar(fontsize=22, col="BLACK",fontface="bold"))
dev.print(pdf, paste0('Figure4B_light.pdf')) # takes ~30sec to open with Adobe
dev.print(svg, paste0('Figure4B_light.svg')) # takes ~1 min to open with inkscape
dev.print(png, paste0('Figure4B_light.png'),width = 1024, height = 768)
dev.off()

save.image("Figure4_data.RData")


#The final manuscript version has been modified on Inkscape


