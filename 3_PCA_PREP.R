######### 3: PCA #########

#Inputs needed:
#2_Landscapes environment
#IF YOU DONT WANT TO SPEND 2H RUNNING IT, YOU WILL NEED DF6.RDS
#Also need df6_2 if dont want to wait 15min

####MUST BE DEFINED BY USER #####
session2<-"your_2_Landscapes_directory_here"
load(paste0(session2,"/env_landscapes.RData")) 

session3<-"your_session_directory_here"
setwd(session3)

####ALL PACKAGES##### not including their dependencies
library(data.table)
#data.table
library(MCRestimate)
#replace.NA
library(png)
#png (save rasters)
library(dplyr)
#bind_cols
library(raster)
#"mask", "rasterFromXYZ", "rasterToPoints"

dir.create("outputs")

#### 3.1: COMPILE DF FOR PCA (CAN BE VERY SLOW) ####
#lz is the list which contains the 179 df with x, y and z
norm_listxyz->lz
names(lz)<-v_species[1:length(lz)]
which(lapply(lz,class)=="NULL")->remove
lz2<-lz[-remove]  #we remove the empty dfs->128 now
#! From now on, number IDs don't correspond anymore
lapply(lz2, function(i) paste(i$x,i$y,sep="/"))->xy #create list of columns with each xy coordinate together in one case only
lz3 <- mapply(cbind, lz2, "xy"=xy, SIMPLIFY=F) #a bit slow ~10 sec
 
as.data.frame(n_means)->means #means=multispecies df
means$xy<-paste(means$x,means$y,sep="/")
t(means)->tmeans
colnames(tmeans)<-tmeans[4,]
data.table((t(as.data.frame(tmeans[3,]))))->tmeans
tmeans<-cbind(data.table(Species="All_species",tmeans)) #tmeans=line with all the multisp values

tmeans->df6  
#View(df6[,1:20])


###!!!! IF YOU DONT WANT TO SPEND 2H RUNNING THIS, LOAD DF6.RDS!!!###
#print("CREATE COMPLETE DF WITH ALL SPECIES AND MULTISP ##VERY SLOW, LIKE 1-2 HOURS")
#for (i in 1:length(lz3)){
 # lz3[[i]]->dfi
#  t(dfi)->dfi
#  colnames(dfi)<-dfi[4,]
#  data.table(t(as.data.frame(dfi[3,])))->dfi
#  dfi<-cbind(data.table(Species=names(lz3[i])),dfi)
#  merge(data.table(df6),data.table(dfi),all=T,by=intersect(colnames(df6),colnames(dfi)))->df6
#  print(i)
#}

###LOAD DF6.RDATA
readRDS(file = "df6.rds")->df6
### KEEP GOING FROM HERE, YOU'RE SAFE NOW

#View(df6[,1:20])
print("TOTAL DF FOR PCA WITH NAs COMPLETED")

#### 3.2: FORMAT THE DF FOR PCA ####

as.data.frame(df6)->df6 #keeping it as data.table causes issues
rownames(df6)<-df6$Species
df6[,-1]->df6
#View(df6[,1:20])
#df6=complete DF with maximum extant, all species and multispecies (=129 rows).

##now we replace the NAs

###LOAD DF6_2.RDS IF YOU DONT WANT TO WAIT 15MIN
#library(MCRestimate)
#replace.NA(df6, as.numeric(df6[5,]), byRow = F)->df6_2
#print("replace all NAs by multisp value, takes ~10-15min")
##View(df6[,1:20])
##View(df6_2[,1:20])

###LOAD DF6_2.RDATA
readRDS(file = "df6_2.rds")->df6_2
### KEEP GOING FROM HERE, YOU'RE SAFE NOW

df6_3 <- as.data.frame(sapply(df6_2, as.numeric)) #MUST convert everything in numeric, must be df6_2 or we lose the rownames
rownames(df6_3)<-rownames(df6_2)
#View(df6_3[,1:20])

#I put multispecies in last row to plot it last
rbind(df6_3[1:4,],df6_3[6:129,],df6_3[5,])->df6_4
#View(df6_4[,1:20])

print("8TH STEP: REPLACE NAS COMPLETED")

#NO NEED TO RE-RUN EVERY TIME/
#Save rasters to display on PCA #par(bg=NA) is for transparence (=no white background)
names(list_rasters)<-v_species[1:length(list_rasters)] #not sure if necessary but I'm paranoid
#Remove rasters with null values
which(lapply(list_rasters,class)=="NULL")->remove
lr2<-list_rasters[-remove]
par(mfrow=c(3,5))
###The loop bellow is just to take a quick look
for(i in 1:15){
  plot(lr2[[i]])
}
par(mfrow=c(3,5))
for(i in 16:30){
  plot(lr2[[i]])
}

print("saving landscape rasters for pca")
dir.create("outputs/PNGS_129/")
#library(png)
for(i in 1:length(lr2)){
  png(paste0("outputs/PNGS_129/",names(lr2[i]),".png"))
  par(bg=NA)
  plot(lr2[[i]],col=cols,zlim = c(0,1),legend=F,axes=F,bty="n",box=FALSE) # box is used to avoid the black "margins" around!
  print(i)
  dev.off()}

png(paste0("outputs/PNGS_129/","Z_All_species.png"))
par(bg=NA)
mar=c(0,0,0,0)
plot(rnorm_xyz_total.masked,col=cols,zlim = c(0,1),legend=F,axes=F,bty="n",box=FALSE)
dev.off()
#/END OF NO NEED TO RE-RUN EVERYTIME

print("9th STEP: LANDSCAPE RASTERS SAVED")

#### 3.4: CROP DF VALUES TO MARONI BASIN BOUNDARIES ####

print("10th STEP: CROP DF VALUES TO MARONI BASIN")

#df6_4: complete df with all species and multisp, NAs replaced by multisp values, multisp line last in df for better plotting

df6_4->df8
df8->df18
#those numbers correspond to the different trials I did, they got dragged in the code

#1: add values of x and y for each pixel

as.numeric(gsub('/(.*)',"", colnames(df18)))->x
as.numeric(gsub(".*/", "", colnames(df18)))->y
rbind(df18,x=x)->df18_2
rbind(df18_2,y=y)->df18_2

#if the bellow is TRUE you are on the right path:
#max(y)=5.670833=NORTH
#min(y)=2.670833=SOUTH
#min(x)=-55.4375=WEST
#max(x)=-53.1625=EAST

#To scan the basin "by row, from left to right":
#We want to sort first by y in descending order then by x in ascending order

df18_3<-df18_2[,order(-df18_2[131,],df18_2[130,])]

#TESTS VISUALISATION:
as.matrix(df18_3)->m_3
#View(m_3[,1:20])
#scatter_fill(m_3[130,],m_3[131,],m_3[4,],nlevels=10,main="TEST",pch=".",cex=8)

#We need to take each sp, change them in raster, then Maroni-crop and reput them together.
#I think this is because of the resolution

mylistxyz <- vector("list", 129)
for (i in 1:129){
  data.frame(x=m_3[130,],y=m_3[131,],Z=m_3[i,])->mylistxyz[[i]]}

l<- mylistxyz
names(l)<-rownames(df18_3[1:129,])

l2 <- vector("list", 129)
for (i in 1:129){
  rasterFromXYZ(l[[i]])->rr
  rr.masked <- mask(x=rr, mask=myshp.r)
  rr.masked->l2[[i]]
}

l3<- vector("list", 129)
for (i in 1:129){
  rasterToPoints(l2[[i]])->l3[[i]]}

l4<- vector("list", 129)
for (i in 1:129){
  l4[[i]]<-l3[[i]][,3]}
names(l)->names(l4)

#library(dplyr)
bind_cols(l4)->df18_4
#viou(df18_4)  ! very heavy

cbind(df18_4,l3[[1]][,1],l3[[1]][,2])->df18_5
colnames(df18_5)<-c(names(l),"x","y")
#View(df18_5)

df18_6<-as.data.frame(t(df18_5))
#View(df18_6[,1:20]) 

as.matrix(df18_6)->m_6

print("DF HAS BEEN CROPPED TO HAVE NO PIXEL FROM OUTSIDE OF THE BASIN BOUNDARIES")


#### 3.5: SCATTER FILL FUNCTION SETUP ####
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
  col <- colorRampPalette(c("blue","deepskyblue","darkturquoise","green","yellow","orange","red"))(nlevels)  
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
#### 3.6 USE SCATTER FILL FUNCTION AS A TEST FOR OUR NEW DF WITH BOUNDARIES#####
dev.off()
scatter_fill(m_6[130,],m_6[131,],m_6[129,],nlevels=10,main="TEST",pch=".",cex=8)

df18_6->df18

print("DATA FRAME READY FOR PCA")
save.image("env_3_PCA_prep_all_data.RData")

#### package summary and unload ####

library(NCmisc)
list.functions.in.file("3_PCA_PREP.R")

###TO DO: DETACH UNUSED PACKAGES

#needed for next steps: 
saveRDS(df18, file = "df18.rds")
rm(list = ls())
