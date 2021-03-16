###2: PHYLIN ANALYSIS AND MAPS###

#Script 1 modified for Figure S3: Simulate a distribution restricted to the outer parts of the sampling area for a few selected species 

#Inputs needed:
#S2 Dist_tables environment
#Figure 1 environment (Figure1_data.RData)
#Maroni shapefile environment (was drawn by hand, formerly a) 
###Maroni shapefile env is not enough anymore (raster update). Now you also need:
#shapefiles directory with sa_dem_30s.bil

####MUST BE DEFINED BY USER #####
session1<-"my/session/path"
load(paste0(session1,"/env_1_Dist_tables.RData")) #~10 sec
#3 objects

session_Fig1<-"my/session/path"
load(paste0(session_Fig1,"/Figure1_data.RData")) 
#4 objects: basemap, ext_basemap, pal, maroni_basins

session2<-"my/session/path"
setwd(session2)

load("maroni_shapefile.RData")
#1 object: maroni_shapefile

####ALL PACKAGES##### not including their dependencies
library(ade4)
#"mantel.rtest"
library(dplyr)
#"group_by","summarise"
library(mapplots)
#"draw.shape"
library(phylin)
#"midpoints","extract.val","idw"
library(raster)
#"mask","rasterFromXYZ","rasterize","setValues"
library(seqinr)
#"col2alpha"
library(xlsx)
#write.xlsx


#### 2.1: SET UP ENVIRONMENT AND REMOVE SPECIES WITH <2 LOCS ####
dir.create("outputs")
error<-c()  #vector to list non-plotted species
filter1<-c()
list()->norm_listxyz #for nomalized multispecies
list()->list_rasters

#data (grid) #table containing grid centroids defining the interpolation area (7955 cells)
basemap<-raster("shapefiles/sa_dem_30s.bil")
as.data.frame(basemap,xy=T)->grid_tot
grid_tot2<-grid_tot[,-3]

coord_tot <- data.frame()#coordinates for multispecies (same for raw and normalized)
cols<-c("darkblue","blue","deepskyblue","darkturquoise","green","yellow","orange","orangered","red","firebrick") #colors for figures
sapply(cols,col2alpha,alpha=0.7)->cols_tr #make colors transparent

print("FIRST STEP START:remove species <2 locs")
seq1<-1:length(list_mat_dist)
for (i in seq1) { #append species to error file if less than 3 specimens
  
  if  (length(list_att_tables[[i]]$Latitude)<=2) { #       (diff(abs(list_att_tables[[i]]$Latitude))<=0.0001)
    append(error,paste0(i," ",v_species[i]," less than 3 specimens"))->error
    errorfile<-file("outputs/file_errors.txt")
    writeLines(error, errorfile)
    close(errorfile)
    append(filter1,i)->filter1
    
  }} 
append(error,paste0("n=",length(error)))->error
errorfile<-file("outputs/file_errors.txt")
writeLines(error, errorfile)
close(errorfile)
print("FIRST STEP COMPLETE:remove species <2 locs")


#### 2.2: PLOT ZERODIST MAPS ONLY ####
dir.create("outputs/zerodist")
seq1 [! seq1 %in% filter1]->seq2

zerodists<-c()
file_zero_dist<-c()
error2<-c()

for(i in seq2){
  if (all(list_mat_dist[[i]]==0)){
    append(zerodists,i)->zerodists
  }}

for(i in zerodists)
{
  #distance matrix i input=list_mat_dist
  matdist<-as.matrix(list_mat_dist[[i]])
  
  #table of coordinates i + field names
  coord<-as.data.frame(list_att_tables[[i]])
  rownames(coord)<-coord$Sample_ID
  coord$lin<-coord$numcode
  coord<-coord[,-c(1:2,5)]
  coord<-coord[,c(2,1,3)]
  colnames(coord)<-c("x","y","lin")
  
  #make new boundaries for grid based on the most extreme coordinates of specimens
  min(coord$x)->minx
  min(coord$y)->miny
  max(coord$x)->maxx
  max(coord$y)->maxy
  
  grid2<- subset(grid_tot2, x >= minx & y >= miny)
  grid2<- subset(grid2, x <= maxx & y <= maxy)
  #test grid
  #plot(grid2, cex=0.01, asp=1, main="Grid of pixels for interpolation",xlab="Longitude", ylab="Latitude")
  
  #If grid2 has 0 pixels because specimens too close -> error
  if  (length(grid2$x)==0) { #       (diff(abs(list_att_tables[[i]]$Latitude))<=0.0001)
    append(error2,paste0(i," ",v_species[i]," grid has zero pixels"))->error2
    error2file<-file("outputs/file_errors2.txt")
    writeLines(error2, error2file)
    close(error2file)
  }
  
  else {
    
    #IDW
    
    mp <- midpoints(coord[,1:2])
    #r.dist <- dist(coord[,1:2])
    
    mp$z <- extract.val(matdist, mp[,1:2])
    int <- idw(mp[,5], mp[,3:4], grid2) 
    label_raw<-'K2P distances'
    label_norm<-'Normalized K2P distances'
    
    #heatmap of K2P as raster
    xyz<-cbind(grid2,int)
    
    if  (length(unique(xyz$x))==1 | length(unique(xyz$y))==1 ) { #       #error2 if grid has only one pixel
      append(error2,paste0(i," ",v_species[i]," grid has only one pixel"))->error2
      error2file<-file("outputs/file_errors2.txt")
      writeLines(error2, error2file)
      close(error2file)
    } else {
      ## bellow if we want to plot the raw values, but I just want the normalized ones
      #rasterFromXYZ(xyz)->rxyz
      #créer un "masque" pour croper au bassin
      #crop <- setValues(rxyz, NA)
      #myshp.r <- rasterize(a, crop)
      #nouveau raster K2P cropé
      #rxyz.masked <- mask(x=rxyz, mask=myshp.r)
      
      #Plot pure K2p distS
      #pdf(paste0("outputs/zerodist/raw_",v_species[i],".pdf"))
      #plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main=v_species[i])
      #plot(maroni_basins,add=T) #Bassins versants noir
      #library(mapplots)
      #draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
      #plot(rxyz.masked,col=cols_tr,add=T,legend=F)
      #points(coord[,1:2], pch=16)
      #plot(rxyz.masked,legend.only=TRUE,col=cols,legend.args=list(text=label_raw, side=2, font=2, line=0.3, cex=1))
      #dev.off()
      
      #Do plot with ormalized distances
      #note: since all distances are null it is not actually normalized, but we call them like
      #that for consistency
      norm_int<-int
      
      #heatmap of normalized interpolation of K2P as raster. 
      norm_xyz<-cbind(grid2,norm_int)
      
      rasterFromXYZ(norm_xyz)->norm_rxyz
      # Crop=dummy raster with a spatial extension equal to the cropped raster, but full of NA values
      crop <- setValues(norm_rxyz, NA)
      #Rasterize the catchment boundaries, with NA outside the catchment boundaries
      myshp.r <- rasterize(maroni_shapefile, crop)
      norm_rxyz.masked <- mask(x=norm_rxyz, mask=myshp.r)
      
      #Plot normalized K2P dist on pdfs
      pdf(paste0("outputs/zerodist/norm_",v_species[i],".pdf"))
      plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main=v_species[i])
      plot(maroni_basins,add=T) #Bassins versants noir
      draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
      plot(norm_rxyz.masked,col=cols_tr,add=T,legend=F,zlim = c(0,1))
      points(coord[,1:2], pch=16)
      plot(norm_rxyz.masked,legend.only=TRUE,col=cols,zlim = c(0,1),legend.args=list(text=label_norm, side=2, font=2, line=0.3, cex=1))
      dev.off()
      
      #Plot normalized K2P dist on pngs
      png(paste0("outputs/zerodist/norm_",v_species[i],".png"))
      plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main=v_species[i])
      plot(maroni_basins,add=T) #Bassins versants noir
      draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
      plot(norm_rxyz.masked,col=cols_tr,add=T,legend=F,zlim = c(0,1))
      points(coord[,1:2], pch=16)
      plot(norm_rxyz.masked,legend.only=TRUE,col=cols,zlim = c(0,1),legend.args=list(text=label_norm, side=2, font=2, line=0.3, cex=1))
      dev.off()
      
      
      #Add xyz and coords to list and df for multispecies:
      coord_tot<-rbind(coord_tot,coord)
      #listxyz[[i]]<-xyz #for raw multispecies
      norm_listxyz[[i]]<-norm_xyz #for norm_multispecies
      #assign(paste0("raster_",v_species[i]),norm_rxyz.masked) #for raster multivariate analysis
      list_rasters[[i]]<-norm_rxyz.masked ##for raster plotting of landscapes
      print(paste0("zerodits ",i,"/",length(v_species)," ",v_species[i]))
      
      #write report
      append(file_zero_dist,paste0(i," ",v_species[i]))->file_zero_dist
      file_matdist_zero<-file("outputs/zerodist/file_matdist_zero.txt")
      writeLines(file_zero_dist, file_matdist_zero)
      close(file_matdist_zero)
      
    }}}
#Report:
append(file_zero_dist,paste0("n=",length(file_zero_dist)))->file_zero_dist
file_matdist_zero<-file("outputs/zerodist/file_matdist_zero.txt")
writeLines(file_zero_dist, file_matdist_zero)
close(file_matdist_zero)
print("SECOND STEP COMPLETE: zerodist landscapes")


#### 2.3: CALCULATE MIDPOINTS, REAL DISTANCES, MANTEL TEST AND CORRECTIONS ON ALL REMAININGS####
seq2 [! seq2 %in% zerodists]->seq3

matdists<-list()
coords<-list()
grids<-list()
mps<-list()
r.dists<-list()
index<-c()
pvalues<-c()
filter2<-c()
mantel_res<-c()

for(i in seq3) 
{
  #distance matrix i input=list_mat_dist
  matdists[[i]]<-as.matrix(list_mat_dist[[i]])
  
  #table of coordinates i + field names
  coord<-as.data.frame(list_att_tables[[i]])
  rownames(coord)<-coord$Sample_ID
  coord$lin<-coord$numcode
  coord<-coord[,-c(1:2,5)]
  coord<-coord[,c(2,1,3)]
  colnames(coord)<-c("x","y","lin")
  coords[[i]]<-coord
  
  #make new boundaries for grid based on the most extreme coordinates of specimens
  min(coord$x)->minx
  min(coord$y)->miny
  max(coord$x)->maxx
  max(coord$y)->maxy
  
  grid2<- subset(grid_tot2, x >= minx & y >= miny)
  grid2<- subset(grid2, x <= maxx & y <= maxy)
  grids[[i]]<-grid2
  #test grid
  #plot(grid2, cex=0.01, asp=1, main="Grid of pixels for interpolation",xlab="Longitude", ylab="Latitude")
  
  #If grid2 has 0 pixels because specimens too close -> error2
  if  (length(grid2$x)==0) { #       (diff(abs(list_att_tables[[i]]$Latitude))<=0.0001)
    append(error2,paste0(i," ",v_species[i]," grid has zero pixels"))->error2
    error2file<-file("outputs/file_errors2.txt")
    writeLines(error2, error2file)
    close(error2file)
    append(filter2,i)->filter2
  }
  
  else {
    
    mps[[i]] <- midpoints(coord[,1:2])
    r.dists[[i]] <- dist(coord[,1:2])
    mantel.rtest(r.dists[[i]],list_mat_dist[[i]],nrepet=9999)->mantel
    append(index,i)->index
    append(pvalues,mantel$pvalue)->pvalues
    
    append(mantel_res,paste0(i," ",v_species[i]," ",mantel))->mantel_res
    mantel_resfile<-file("outputs/file_mantel_res.txt")
    writeLines(mantel_res, mantel_resfile)
    close(mantel_resfile)
    
    #for multispecies
    coord_tot<-rbind(coord_tot,coord)
    
    print(paste0("Mantel ",i,"/",length(v_species)," ",v_species[i]))
    
  }}

#Mantel report
append(mantel_res,paste0("n=",length(mantel_res)/7))->mantel_res
mantel_resfile<-file("outputs/file_mantel_res.txt")
writeLines(mantel_res, mantel_resfile)
close(mantel_resfile) 

#Correction for multitest, if significant, extract the residuals from a linear regression, 
#attribute the residuals to the midpoints
p.adjust(pvalues, method = "fdr")->fdr

dfpval<-cbind(v_species[index],index,pvalues,fdr)

colnames(dfpval)<-c("Species","index",paste0("pvalue ",length(which(pvalues<0.05))),paste0("fdr ",length(which(fdr<0.05))))
#library(xlsx)
write.xlsx(dfpval, "outputs/file_p-values.xlsx")
print("THIRD STEP COMPLETE: Mantel tests")

#### 2.4: PLOT EVERYTHING ON FDR ####
dir.create("outputs/fdr")
file_fdr_sign<-c()
file_fdr_not_sign<-c()
for(i in index)
{
  if  (fdr[which(dfpval[,2]==i)]<0.05) { #       
    d.real <- as.matrix(r.dists[[i]]) #real distances to matrix
    fit <- lm(as.vector(matdists[[i]]) ~ as.vector(d.real)) #linear regression bw geo and genetic
    resid <- matrix(fit$residuals, nrow(coords[[i]]), nrow(coords[[i]]))
    dimnames(resid) <- dimnames(matdists[[i]])
    
    #interpolation of residuals
    mps[[i]]$z <- extract.val(resid, mps[[i]][,1:2])
    int <- idw(mps[[i]][,5], mps[[i]][,3:4], grids[[i]])
    
    #fig legend and write report
    label_raw<-'Residuals of K2P vs. real distances'
    label_norm<-'Normalized residuals of K2P vs. real distances'
    
    append(file_fdr_sign,paste0(i," ",v_species[i]))->file_fdr_sign
    file_fdr_significant<-file("outputs/fdr/file_fdr_significant.txt")
    writeLines(file_fdr_sign, file_fdr_significant)
    close(file_fdr_significant)
    
  } else {
    #library(phylin)
    mps[[i]]$z <- extract.val(matdists[[i]], mps[[i]][,1:2])
    int <- idw(mps[[i]][,5], mps[[i]][,3:4], grids[[i]])
    label_raw<-'K2P distances'
    label_norm<-'Normalized K2P distances'
    
    #write report
    append(file_fdr_not_sign,paste0(i," ",v_species[i]))->file_fdr_not_sign
    file_fdr_not_significant<-file("outputs/fdr/file_fdr_not_significant.txt")
    writeLines(file_fdr_not_sign, file_fdr_not_significant)
    close(file_fdr_not_significant)
  }
  
  #heatmap of K2P or residuals as raster
  xyz<-cbind(grids[[i]],int)
  
  if  (length(unique(xyz$x))==1 | length(unique(xyz$y))==1 ) { #       #error2 si la grid n'a qu'un pixel
    append(error2,paste0(i," ",v_species[i]," grid has only one pixel"))->error2
    error2file<-file("outputs/file_errors2.txt")
    writeLines(error2, error2file)
    close(error2file)
  } else {
    
    #Only if we awant to plot RAW K2p dist or RAW residuals
    #rasterFromXYZ(xyz)->rxyz
    #crop <- setValues(rxyz, NA)
    #myshp.r <- rasterize(a, crop)
    #rxyz.masked <- mask(x=rxyz, mask=myshp.r)
    
    #Plot pure K2p dist or residuals
    #pdf(paste0("outputs/fdr/raw_",v_species[i],".pdf"))
    #plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main=v_species[i])
    #plot(maroni_basins,add=T) #Bassins versants noir
    #draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
    #plot(rxyz.masked,col=cols_tr,add=T,legend=F)
    #points(coords[[i]][,1:2], pch=16)
    #plot(rxyz.masked,legend.only=TRUE,col=cols,legend.args=list(text=label_raw, side=2, font=2, line=0.3, cex=1))
    #dev.off()
    
    #SAME WITH NORMALIZED INTERPOLATION OF RESIDUALS
    
    ##Normaliser directement sur les résidus (et non l'interpol) engendre des erreurs, je sais plus pk
    ##norm<-function(resid){(resid-min(resid))/(max(resid)-min(resid))}
    ##norm(resid)->norm_resid
    ##mp->norm_mp
    ##norm_mp$z <- extract.val(norm_resid, norm_mp[,1:2])
    ##norm_int <- idw(norm_mp[,5], norm_mp[,3:4], grid2)
    
    norm_int<-(int-min(int))/(max(int)-min(int))
    
    #heatmap of normalized interpolation of K2P or residuals as raster. 
    norm_xyz<-cbind(grids[[i]],norm_int)
    
    if  (length(unique(norm_xyz$x))==1 | length(unique(norm_xyz$y))==1 ) { #       #error2 si la grid n'a qu'un pixel
      append(error2,paste0(i," ",v_species[i]," norm grid has only one pixel"))->error2
      error2file<-file("outputs/file_errors2.txt")
      writeLines(error2, error2file)
      close(error2file)
    } else {
      
      rasterFromXYZ(norm_xyz)->norm_rxyz
      crop <- setValues(norm_rxyz, NA)
      myshp.r <- rasterize(maroni_shapefile, crop)
      norm_rxyz.masked <- mask(x=norm_rxyz, mask=myshp.r)
      
      #Plot normalized K2P dist or residuals
      pdf(paste0("outputs/fdr/norm_",v_species[i],".pdf"))
      plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main=v_species[i])
      plot(maroni_basins,add=T) #Bassins versants noir
      draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
      plot(norm_rxyz.masked,col=cols_tr,add=T,legend=F,zlim = c(0,1))
      points(coords[[i]][,1:2], pch=16)
      plot(norm_rxyz.masked,legend.only=TRUE,col=cols,zlim = c(0,1),legend.args=list(text=label_norm, side=2, font=2, line=0.3, cex=1))
      dev.off()
      
      #Plot normalized K2P dist or residuals on png
      png(paste0("outputs/fdr/norm_",v_species[i],".png"))
      plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main=v_species[i])
      plot(maroni_basins,add=T) #Bassins versants noir
      draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
      plot(norm_rxyz.masked,col=cols_tr,add=T,legend=F,zlim = c(0,1))
      points(coords[[i]][,1:2], pch=16)
      plot(norm_rxyz.masked,legend.only=TRUE,col=cols,zlim = c(0,1),legend.args=list(text=label_norm, side=2, font=2, line=0.3, cex=1))
      dev.off()
      
      
      #Add xyz and coords to list and df for multispecies:
      #listxyz[[i]]<-xyz#for raw multispecies
      norm_listxyz[[i]]<-norm_xyz #for norm multispecies
      #assign(paste0("raster_",v_species[i]),norm_rxyz.masked)
      list_rasters[[i]]<-norm_rxyz.masked #for raster plotting of landscapes
      print(paste0("fdr ",i,"/",length(v_species)," ",v_species[i]))
      
    }}}

#write reports:
append(file_fdr_not_sign,paste0("n=",length(file_fdr_not_sign)," including non-plotted errors"))->file_fdr_not_sign
file_fdr_not_significant<-file("outputs/fdr/file_fdr_not_significant.txt")
writeLines(file_fdr_not_sign, file_fdr_not_significant)
close(file_fdr_not_significant)

append(file_fdr_sign,paste0("n=",length(file_fdr_sign)," including non-plotted errors"))->file_fdr_sign
file_fdr_significant<-file("outputs/fdr/file_fdr_significant.txt")
writeLines(file_fdr_sign, file_fdr_significant)
close(file_fdr_significant)

print("FOURTH STEP COMPLETE: fdr+all landscapes")
