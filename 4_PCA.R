######### 4: PCA #########

#Inputs needed:
#df18.rds (output from 3_landscapes)
#PNGS_129 directory 
#the folder "coverage". Files in that folder were sorted by hand (formerly: OLD_MAPS_PICS)

####MUST BE DEFINED BY USER #####
session3<-"your_session3_directory_here"
session3<-"E:/BARCODE_MARONI/SCRIPTS_AND_DATA/LANDSCAPES_PIPELINE/3_PCA_PREP"
readRDS(paste0(session3,"/df18.rds"))->df18

session4<-"your_session_directory_here"
session4<-"E:/BARCODE_MARONI/SCRIPTS_AND_DATA/LANDSCAPES_PIPELINE/4_PCA"
setwd(session4)

#### 4.1: DO THE PCAs #####

library(ade4)
dir.create("outputs/PCA")

#### All locs ####
pca18<-dudi.pca(df18[1:129,], center=T, scale=F, scannf = F, nf=5) #~5 seconds
#do not do scatter(pca18,1,2), too many variables to plot!
par(mfrow=c(1,1))
#scatter(pca18,1,2)

barplot(pca18$eig[1:6])
png(paste0("outputs/PCA/","eig.png"))
barplot(pca18$eig[1:6])
dev.off()
pdf(paste0("outputs/PCA/","eig.pdf"))
barplot(pca18$eig[1:6])
dev.off()

s.label(pca18$li,1,2)
png(paste0("outputs/PCA/","labels.png"))
s.label(pca18$li,1,2)
dev.off()
pdf(paste0("outputs/PCA/","labels.pdf"))
s.label(pca18$li,1,2)
dev.off()

summary(pca18)
sink("outputs/PCA/summary.txt")
print(summary(pca18))
sink()

library(png)
lpix2 <- list()
for (nomfic in list.files("PNGS_129/", pattern = ".png")) {nomobj <- strsplit(nomfic, "[.]")[[1]][1]
toto <- readPNG(paste("PNGS_129/", nomfic, sep = ""))
toto <- list(toto)
names(toto) <- nomobj
lpix2 = c(lpix2, toto)}

thumbnails <- function(x, y, images, width, height){
  for (ii in seq_along(x)){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

#AXES 1 & 2
png(paste0("outputs/PCA/","PCA_12.png"))
plot(pca18$li[,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.035*diff(range(pca18$li[,1])),height = 0.05*diff(range(pca18$li[,2])))
dev.off()

pdf(paste0("outputs/PCA/","PCA_12.pdf"))
plot(pca18$li[,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.035*diff(range(pca18$li[,1])),height = 0.05*diff(range(pca18$li[,2])))
dev.off()

#AXES 2 & 3
png(paste0("outputs/PCA/","PCA_23.png"))
plot(pca18$li[,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[,2])-10,max(pca18$li[,2])+10)), ylim=(c(min(pca18$li[,3])-10, max(pca18$li[,3])+10)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.025*diff(range(pca18$li[,1])),height = 0.04*diff(range(pca18$li[,2])))
dev.off()

pdf(paste0("outputs/PCA/","PCA_23.pdf"))
plot(pca18$li[,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[,2])-10,max(pca18$li[,2])+10)), ylim=(c(min(pca18$li[,3])-10, max(pca18$li[,3])+10)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.025*diff(range(pca18$li[,1])),height = 0.04*diff(range(pca18$li[,2])))
dev.off()

#AXES 4 & 5
png(paste0("outputs/PCA/","PCA_45.png"))
plot(pca18$li[,4:5], t="n",axes=TRUE, xlim=(c(min(pca18$li[,4])-10,max(pca18$li[,4])+10)), ylim=(c(min(pca18$li[,5])-10, max(pca18$li[,5])+10)), xlab="PC4", ylab= "Axis 5")
abline(h = 0, v = 0)
thumbnails(pca18$li[,4], pca18$li[,5], lpix2, width = 0.02*diff(range(pca18$li[,1])),height = 0.04*diff(range(pca18$li[,2])))
dev.off()

pdf(paste0("outputs/PCA/","PCA_45.pdf"))
plot(pca18$li[,4:5], t="n",axes=TRUE, xlim=(c(min(pca18$li[,4])-10,max(pca18$li[,4])+10)), ylim=(c(min(pca18$li[,5])-10, max(pca18$li[,5])+10)), xlab="PC4", ylab= "Axis 5")
abline(h = 0, v = 0)
thumbnails(pca18$li[,4], pca18$li[,5], lpix2, width = 0.02*diff(range(pca18$li[,1])),height = 0.04*diff(range(pca18$li[,2])))
dev.off()

#Variables?
#plot(pca18$co[,1:2], t="p",axes=TRUE, xlim=(c(min(pca18$co[,1])-0.01,max(pca18$co[,1])+0.01)), ylim=(c(min(pca18$co[,2])-0.01, max(pca18$co[,2])+0.01)), xlab="PC1", ylab="PC2")
#abline(h = 0, v = 0)
#text(pca18$co[,1:2],labels=colnames(df18))

#### Latitude coordinates as variables ####
dir.create("outputs/PCA/coordinates_variables")

## SCATTER FILL FUNCTION SETUP ##
scatter_fill <- function (x, y, z,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),zlim=c(min(z),max(z)),
                          nlevels = 20, plot.title, plot.axes, 
                          key.title, key.axes, asp = NA, xaxs = "i", 
                          yaxs = "i", las = 1, 
                          axes = TRUE, frame.plot = axes, ...) 
{
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  
  # choose colors to interpolate
  levels <- seq(zlim[1],zlim[2],length.out = nlevels)
  col <- colorRampPalette(coord_colors)(nlevels)  
  colz <- col[cut(z,nlevels)]  
  #   
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
  
  rect(0, levels[-length(levels)], 1, levels[-1L],col=col,border=col) 
  if (missing(key.axes)) {if (axes){axis(4)}}
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  
  # points
  plot(x,y,type = "n",xaxt='n',yaxt='n',xlab="",ylab="",xlim=xlim,ylim=ylim,bty="n")
  points(x,y,col = colz,xaxt='n',yaxt='n',xlab="",ylab="",bty="n",...)
  
  ## options to make mapping more customizable
  
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

c("red","orange","yellow")->coord_colors

png(paste0("outputs/PCA/coordinates_variables/","projections_latitude12.png"))
scatter_fill(pca18$co[,1],pca18$co[,2],as.numeric(df18[131,]),main="Projections of latitudes",pch=".",cex=3,xlab="PC1",ylab="PC2",cex.main=1.8,cex.lab=1.5)
dev.off()
pdf(paste0("outputs/PCA/coordinates_variables/","projections_latitude12.pdf"))
scatter_fill(pca18$co[,1],pca18$co[,2],as.numeric(df18[131,]),main="Projections of latitudes",pch=".",cex=3,xlab="PC1",ylab="PC2",cex.main=1.8,cex.lab=1.5)
dev.off()

png(paste0("outputs/PCA/coordinates_variables/","legend_latitude.png"))
scatter_fill(as.numeric(df18[130,]),as.numeric(df18[131,]),as.numeric(df18[131,]),main="Colours latitudes basin",pch=".",cex=3,xlab="lon",ylab="lat",cex.main=1.8,cex.lab=1.5)
dev.off()
pdf(paste0("outputs/PCA/coordinates_variables/","legend_latitude.pdf"))
scatter_fill(as.numeric(df18[130,]),as.numeric(df18[131,]),as.numeric(df18[131,]),main="Colours latitudes basin",pch=".",cex=3,xlab="lon",ylab="lat",cex.main=1.8,cex.lab=1.5)
dev.off()

png(paste0("outputs/PCA/coordinates_variables/","projections_latitude23.png"))
scatter_fill(pca18$co[,2],pca18$co[,3],as.numeric(df18[131,]),main="Projections of latitudes",pch=".",cex=3,xlab="PC2",ylab="PC3",cex.main=1.8,cex.lab=1.5)
dev.off()
pdf(paste0("outputs/PCA/coordinates_variables/","projections_latitude23.pdf"))
scatter_fill(pca18$co[,2],pca18$co[,3],as.numeric(df18[131,]),main="Projections of latitudes",pch=".",cex=3,xlab="PC2",ylab="PC3",cex.main=1.8,cex.lab=1.5)
dev.off()


#### Longitude coordinates as variables ####
df18we<-df18[,order(df18[130,],-df18[131,])]
#View(df18we[,1:20])
#!!! pca18we gives exactly the same results as pca18, this has been thoroughly tested.
#The only difference is that it is ll ordered by longitude
pca18we<-dudi.pca(df18we[1:129,], center=T, scale=F, scannf = F, nf=5) #~5 seconds

c("green","blue","purple")->coord_colors
sizes<-c(1.5,1.5,1.5)

png(paste0("outputs/PCA/coordinates_variables/","projections_longitude12.png"))
scatter_fill((pca18we$co[,1]),pca18we$co[,2],as.numeric(df18[130,]),main="Projections of longitudes",pch=".",cex=3,xlab="PC1",ylab="PC2",cex.main=1.8,cex.lab=1.5)
dev.off()
pdf(paste0("outputs/PCA/coordinates_variables/","projections_longitude12.pdf"))
scatter_fill((pca18we$co[,1]),pca18we$co[,2],as.numeric(df18[130,]),main="Projections of longitudes",pch=".",cex=3,xlab="PC1",ylab="PC2",cex.main=1.8,cex.lab=1.5)
dev.off()

png(paste0("outputs/PCA/coordinates_variables/","legend_longitude.png"))
scatter_fill(as.numeric(df18we[130,]),as.numeric(df18we[131,]),as.numeric(df18we[130,]),main="Colours longitudes basin",pch=".",cex=3,xlab="lon",ylab="lat",cex.main=1.8,cex.lab=1.5)
dev.off()
pdf(paste0("outputs/PCA/coordinates_variables/","legend_longitude.pdf"))
scatter_fill(as.numeric(df18we[130,]),as.numeric(df18we[131,]),as.numeric(df18we[130,]),main="Colours longitudes basin",pch=".",cex=3,xlab="lon",ylab="lat",cex.main=1.8,cex.lab=1.5)
dev.off()

png(paste0("outputs/PCA/coordinates_variables/","projections_longitude23.png"))
scatter_fill(pca18we$co[,2],pca18we$co[,3],as.numeric(df18we[130,]),nlevels=1000,main="Projections of longitudes",pch=".",cex=3,xlab="PC2",ylab="PC3",cex.main=1.8,cex.lab=1.5)
dev.off()
pdf(paste0("outputs/PCA/coordinates_variables/","projections_longitude23.pdf"))
scatter_fill(pca18we$co[,2],pca18we$co[,3],as.numeric(df18we[130,]),nlevels=1000,main="Projections of longitudes",pch=".",cex=3,xlab="PC2",ylab="PC3",cex.main=1.8,cex.lab=1.5)
dev.off()

#### West_North_East ####
#NEED FOLDER OLD_MAP_PICS
dir.create("outputs/PCA/coverage/")

list.files("./coverages/West_North_East/")->West_North_East
unlist(lapply(strsplit(West_North_East, "norm_"), '[[', 2))->West_North_East
unlist(lapply(strsplit(West_North_East, "-1"), '[[', 1))->West_North_East

c()->West_North_East_nb
for(i in 1:length(West_North_East)){
  West_North_East_nb<-c(West_North_East_nb,which(rownames(df18)==West_North_East[i]))
}
sort(West_North_East_nb)->West_North_East_nb #to plot All_species last

s.label(pca18$li[West_North_East_nb,2:3])
png(paste0("outputs/PCA/coverage/","whole_coverage_labels23.png"))
s.label(pca18$li[West_North_East_nb,2:3])
dev.off()
pdf(paste0("outputs/PCA/coverage/","whole_coverage_labels23.pdf"))
s.label(pca18$li[West_North_East_nb,2:3])
dev.off()


thumbnails <- function(x, y, images, width,
                       height){
  for (ii in West_North_East_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","whole_coverage_23.png"))
plot(pca18$li[West_North_East_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_North_East_nb,2])-0.2,max(pca18$li[West_North_East_nb,2])+0.2)), ylim=(c(min(pca18$li[West_North_East_nb,3])-1, max(pca18$li[West_North_East_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width=3.2, height=3.8)
dev.off()

pdf(paste0("outputs/PCA/coverage/","whole_coverage_23.pdf"))
plot(pca18$li[West_North_East_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_North_East_nb,2])-0.2,max(pca18$li[West_North_East_nb,2])+0.2)), ylim=(c(min(pca18$li[West_North_East_nb,3])-1, max(pca18$li[West_North_East_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2,width=3.2, height=3.8)
dev.off()

#Axes1-2
s.label(pca18$li[West_North_East_nb,1:2])
png(paste0("outputs/PCA/coverage/","whole_coverage_labels12.png"))
s.label(pca18$li[West_North_East_nb,1:2])
dev.off()
pdf(paste0("outputs/PCA/coverage/","whole_coverage_labels12.pdf"))
s.label(pca18$li[West_North_East_nb,1:2])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in West_North_East_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","whole_coverage_12.png"))
plot(pca18$li[West_North_East_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_North_East_nb,1])-10,max(pca18$li[West_North_East_nb,1])+10)), ylim=(c(min(pca18$li[West_North_East_nb,2])-10, max(pca18$li[West_North_East_nb,2])+10)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2,width=3.2, height=3.8)
dev.off()

pdf(paste0("outputs/PCA/coverage/","whole_coverage_12.pdf"))
plot(pca18$li[West_North_East_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_North_East_nb,1])-10,max(pca18$li[West_North_East_nb,1])+10)), ylim=(c(min(pca18$li[West_North_East_nb,2])-10, max(pca18$li[West_North_East_nb,2])+10)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2,width=3.2, height=3.8)
dev.off()

#### West_East: #####

list.files("./coverages/West_East/")->West_East

unlist(lapply(strsplit(West_East, "norm_"), '[[', 2))->West_East
unlist(lapply(strsplit(West_East, "-1"), '[[', 1))->West_East

c()->West_East_nb
for(i in 1:length(West_East)){
  West_East_nb<-c(West_East_nb,which(rownames(df18)==West_East[i]))
}

s.label(pca18$li[West_East_nb,2:3])
png(paste0("outputs/PCA/coverage/","upper_maroni_labels23.png"))
s.label(pca18$li[West_East_nb,2:3])
dev.off()
pdf(paste0("outputs/PCA/coverage/","upper_maroni_labels23.pdf"))
s.label(pca18$li[West_East_nb,2:3])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in West_East_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","upper_maroni23.png"))
plot(pca18$li[West_East_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_East_nb,2])-0.2,max(pca18$li[West_East_nb,2])+0.2)), ylim=(c(min(pca18$li[West_East_nb,3])-1, max(pca18$li[West_East_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 2, height = 2)
dev.off()

pdf(paste0("outputs/PCA/coverage/","upper_maroni23.pdf"))
plot(pca18$li[West_East_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_East_nb,2])-0.2,max(pca18$li[West_East_nb,2])+0.2)), ylim=(c(min(pca18$li[West_East_nb,3])-1, max(pca18$li[West_East_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 2, height = 2)
dev.off()

#Axes1-2
s.label(pca18$li[West_East_nb,1:2])
png(paste0("outputs/PCA/coverage/","upper_maroni_labels12.png"))
s.label(pca18$li[West_East_nb,1:2])
dev.off()
pdf(paste0("outputs/PCA/coverage/","upper_maroni_labels12.pdf"))
s.label(pca18$li[West_East_nb,1:2])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in West_East_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","upper_maroni12.png"))
plot(pca18$li[West_East_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_East_nb,1])-0.02,max(pca18$li[West_East_nb,1])+0.02)), ylim=(c(min(pca18$li[West_East_nb,2])-0.2, max(pca18$li[West_East_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=2, height=2)
dev.off()

pdf(paste0("outputs/PCA/coverage//","upper_maroni12.pdf"))
plot(pca18$li[West_East_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[West_East_nb,1])-0.02,max(pca18$li[West_East_nb,1])+0.02)), ylim=(c(min(pca18$li[West_East_nb,2])-0.2, max(pca18$li[West_East_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=2, height=2)
dev.off()

##### East ####

list.files("./coverages/East/")->East

unlist(lapply(strsplit(East, "norm_"), '[[', 2))->East
unlist(lapply(strsplit(East, "-1"), '[[', 1))->East

c()->East_nb
for(i in 1:length(East)){
  East_nb<-c(East_nb,which(rownames(df18)==East[i]))
}

s.label(pca18$li[East_nb,2:3])
png(paste0("outputs/PCA/coverage/","east_maroni_label23.png"))
s.label(pca18$li[East_nb,2:3])
dev.off()
pdf(paste0("outputs/PCA/coverage/","east_maroni_label23.pdf"))
s.label(pca18$li[East_nb,2:3])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in East_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","east_maroni23.png"))
plot(pca18$li[East_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[East_nb,2])-0.2,max(pca18$li[East_nb,2])+0.2)), ylim=(c(min(pca18$li[East_nb,3])-1, max(pca18$li[East_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.2, height = 1)
dev.off()

pdf(paste0("outputs/PCA/coverage/","east_maroni23.pdf"))
plot(pca18$li[East_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[East_nb,2])-0.2,max(pca18$li[East_nb,2])+0.2)), ylim=(c(min(pca18$li[East_nb,3])-1, max(pca18$li[East_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.2, height = 1)
dev.off()

#Axes1-2
s.label(pca18$li[East_nb,1:2])
png(paste0("outputs/PCA/coverage/","east_maroni_labels12.png"))
s.label(pca18$li[East_nb,1:2])
dev.off()
pdf(paste0("outputs/PCA/coverage/","east_maroni_labels12.pdf"))
s.label(pca18$li[East_nb,1:2])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in East_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","east_maroni12.png"))
plot(pca18$li[East_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[East_nb,1])-0.02,max(pca18$li[East_nb,1])+0.02)), ylim=(c(min(pca18$li[East_nb,2])-0.2, max(pca18$li[East_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=0.4, height=0.3)
dev.off()

pdf(paste0("outputs/PCA/coverage//","east_maroni12.pdf"))
plot(pca18$li[East_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[East_nb,1])-0.02,max(pca18$li[East_nb,1])+0.02)), ylim=(c(min(pca18$li[East_nb,2])-0.2, max(pca18$li[East_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=0.4, height=0.3)
dev.off()

##### North-East #####

list.files("./coverages/North_East/")->north_east

unlist(lapply(strsplit(north_east, "norm_"), '[[', 2))->north_east
unlist(lapply(strsplit(north_east, "-1"), '[[', 1))->north_east

c()->north_east_nb
for(i in 1:length(north_east)){
  north_east_nb<-c(north_east_nb,which(rownames(df18)==north_east[i]))
}

s.label(pca18$li[north_east_nb,2:3])
png(paste0("outputs/PCA/coverage/","north_east_maroni_label23.png"))
s.label(pca18$li[north_east_nb,2:3])
dev.off()
pdf(paste0("outputs/PCA/coverage/","north_east_maroni_label23.pdf"))
s.label(pca18$li[north_east_nb,2:3])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in north_east_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","north_east_maroni23.png"))
plot(pca18$li[north_east_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[north_east_nb,2])-0.2,max(pca18$li[north_east_nb,2])+0.2)), ylim=(c(min(pca18$li[north_east_nb,3])-1, max(pca18$li[north_east_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 1, height = 3)
dev.off()

pdf(paste0("outputs/PCA/coverage/","north_east_maroni23.pdf"))
plot(pca18$li[north_east_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[north_east_nb,2])-0.2,max(pca18$li[north_east_nb,2])+0.2)), ylim=(c(min(pca18$li[north_east_nb,3])-1, max(pca18$li[north_east_nb,3])+1)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 1, height = 3)
dev.off()

#Axes1-2
s.label(pca18$li[north_east_nb,1:2])
png(paste0("outputs/PCA/coverage/","north_east_maroni_labels12.png"))
s.label(pca18$li[north_east_nb,1:2])
dev.off()
pdf(paste0("outputs/PCA/coverage/","north_east_maroni_labels12.pdf"))
s.label(pca18$li[north_east_nb,1:2])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in north_east_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","north_east_maroni12.png"))
plot(pca18$li[north_east_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[north_east_nb,1])-0.02,max(pca18$li[north_east_nb,1])+0.02)), ylim=(c(min(pca18$li[north_east_nb,2])-0.2, max(pca18$li[north_east_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=0.8, height=1.5)
dev.off()

pdf(paste0("outputs/PCA/coverage//","north_east_maroni12.pdf"))
plot(pca18$li[north_east_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[north_east_nb,1])-0.02,max(pca18$li[north_east_nb,1])+0.02)), ylim=(c(min(pca18$li[north_east_nb,2])-0.2, max(pca18$li[north_east_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=0.8, height=1.5)
dev.off()

#### Others ####

list.files("./coverages/Others/")->others

unlist(lapply(strsplit(others, "norm_"), '[[', 2))->others
unlist(lapply(strsplit(others, "-1"), '[[', 1))->others

c()->others_nb
for(i in 1:length(others)){
  others_nb<-c(others_nb,which(rownames(df18)==others[i]))
}

s.label(pca18$li[others_nb,2:3])
png(paste0("outputs/PCA/coverage/","others_maroni_label23.png"))
s.label(pca18$li[others_nb,2:3])
dev.off()
pdf(paste0("outputs/PCA/coverage/","others_maroni_label23.pdf"))
s.label(pca18$li[others_nb,2:3])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in others_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","others_maroni23.png"))
plot(pca18$li[others_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[others_nb,2]),max(pca18$li[others_nb,2])+0.2)), ylim=(c(min(pca18$li[others_nb,3]), max(pca18$li[others_nb,3]))), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.025, height = 0.08)
dev.off()

pdf(paste0("outputs/PCA/coverage/","others_maroni23.pdf"))
plot(pca18$li[others_nb,2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[others_nb,2]),max(pca18$li[others_nb,2])+0.2)), ylim=(c(min(pca18$li[others_nb,3]), max(pca18$li[others_nb,3]))), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.025, height = 0.08)
dev.off()

#Axes1-2
s.label(pca18$li[others_nb,1:2])
png(paste0("outputs/PCA/coverage/","others_maroni_labels12.png"))
s.label(pca18$li[others_nb,1:2])
dev.off()
pdf(paste0("outputs/PCA/coverage/","others_maroni_labels12.pdf"))
s.label(pca18$li[others_nb,1:2])
dev.off()

thumbnails <- function(x, y, images, width,
                       height){
  for (ii in others_nb){
    rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                ybottom= y[ii] - 1*height,
                xright=x[ii] + 1*width,
                ytop= y[ii] + 1*height, interpolate=FALSE)
  }
}

png(paste0("outputs/PCA/coverage/","others_maroni12.png"))
plot(pca18$li[others_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[others_nb,1])-0.02,max(pca18$li[others_nb,1])+0.02)), ylim=(c(min(pca18$li[others_nb,2])-0.2, max(pca18$li[others_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=0.015, height=0.04)
dev.off()

pdf(paste0("outputs/PCA/coverage//","others_maroni12.pdf"))
plot(pca18$li[others_nb,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[others_nb,1])-0.02,max(pca18$li[others_nb,1])+0.02)), ylim=(c(min(pca18$li[others_nb,2])-0.2, max(pca18$li[others_nb,2])+0.2)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width=0.015, height=0.04)
dev.off()


save.image("env_4_PCA_all_data.RData")


