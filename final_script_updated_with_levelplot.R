####1: CREATE DIST TABLES FROM DNA AND LOCS ----
#Need: 
#1)Maroni_loc.txt
#2)all_seq.fas
#3)Fig1_env.RData (last one only for plotting maps)
#4) a folder "shapefiles" with sa_dem_30s files to avoid Error in .local(.Object, ...) (new raster version)
#5) tHE FOLDER OLD_MAP_PICS if you want PCA per coverage
#Useful tip if all_seq.fas comes fom BOLD: use GBOL.*?[|] to remove the GBOL characters
#Other useful tip: Check carefully that species names match between Maroni loc and all_seq! (ex: G. brevispinnis...)
setwd("E:/BARCODE/TRAVAIL V4/1905 MAJ FINALE/LANDSCAPE FINAL SCRIPT")
library(ape)

dir.create("outputs")
dir.create("outputs/tables")
dir.create("outputs/tables/matdist")
dir.create("outputs/tables/att_tables")

#!!!!This part until PHYLIN ANALYSIS AND FIGS compute distance matrices and other useful tables for
#each species. If you don't want to compute it everytime you produce
#maps just save the results as env_all_tables.Rdata and load it before PHYLIN ANALYSIS AND FIGS

#dfloc: df with 5 c: Species, Lat, Lon, field num, exact site, ordered by field.num
dfloc<-read.table("Maroni_locv5.txt",header = T)
dfloc<- dfloc[order(dfloc$field.num),] 

#dfdna: input=dna sequences; output=df with 2 c: 1)Species|field.numb 2)seq as character strings
#Do not forget to change spaces into underscores before
#AND remove BOLD identifiers with GBOL.*?[|]
#library(ape)
read.dna("all_seq_v5.fas",format="fasta",as.character=T)->myDNA #heavy
dfdna<-as.data.frame(myDNA) #heavy
colnames(dfdna)<-seq(1:length(dfdna[1,]))
dfdna$seq <- apply( dfdna[ , colnames(dfdna) ] , 1 , paste , collapse = "" )
dfdna <- as.data.frame(dfdna[,-(1:length(dfdna)-1),drop = FALSE]) #remove number of columns -1
dfdna<-cbind(rownames(dfdna),dfdna)
dfdna <- data.frame(lapply(dfdna, as.character), stringsAsFactors=FALSE)
colnames(dfdna)<-c("specimen","seq")

#list_split: create list with splitted species name and field.num
strsplit(as.character(dfdna$specimen),split="\\|")->list_split

#df_split: df with 2c: species and field.num
as.data.frame(t(as.data.frame(list_split)))->df_split
df_split <- data.frame(lapply(df_split, as.character), stringsAsFactors=FALSE)

#dflocdna: merge df_split and dfloc to get df with 5 c: sp, seq, field.num, Lat, Lon, exact site
dflocdna<-data.frame(cbind(df_split[,1],dfdna$seq,df_split[,2]),stringsAsFactors = FALSE)
dflocdna<- dflocdna[order(dflocdna$X3),]
dflocdna<-cbind(dflocdna,dfloc$Lat,dfloc$Lon,dfloc$Exact_Site)
colnames(dflocdna)<-c("sp","seq","field.num","Lat","Lon","Exact_Site")

#list_dflocdna_sp: list of n=number of specimens elements. Each element is a data frame with seq, coord, field num for each species.
#Because of lapply, outputs as many identical df for each species as there are specimens in each species
lapply(dflocdna$sp,function(i) assign(paste("df", i, sep = "."), dflocdna[dflocdna$sp==i,]))->list_dflocdna_sp

#list_unique_dflocdna_sp:list of n=number of species elements. We just removed the identical dfs
unique(list_dflocdna_sp)->list_unique_dflocdna_sp


#Changer les latitudes egales en +0.00001 pour chaque sp
for (j in 1:length(list_unique_dflocdna_sp)) {
  k=0.00001
  for (i in 1:length(list_unique_dflocdna_sp[[j]]$Lat)) {
    if (list_unique_dflocdna_sp[[j]]$Lat[i] %in% list_unique_dflocdna_sp[[j]]$Lat[-i]) {
      list_unique_dflocdna_sp[[j]]$Lat[i]<-list_unique_dflocdna_sp[[j]]$Lat[i]+k
      k=k+0.00001
    }}}
rm(i,j,k)

#Ok, mtnt assigner un numero de population à chaque Latitude dans chaque df
#lapply(list_unique_dflocdna_sp,function(i) transform(i,pop=as.numeric(factor(Lat))))->list_unique_dflocdna_sp

#list_seqmats: list of n= nb of sp matrices. Each matrix contains dna sequences with nb of rows= nb of specimens/sp
lapply(list_unique_dflocdna_sp,function(i) t(sapply(strsplit(i[,2],""), tolower)))->list_seqmats

#list_sp_field.nums (formerly "a"): list of n= nb of sp vectors. Each vector contains the field.num of each specimen
lapply(list_unique_dflocdna_sp,function(i) print(i[,3]))->list_sp_field.nums

#list_seqmats2: attribute field.num to each sequence
list()->list_seqmats2
for (i in 1:length(list_sp_field.nums)) {
  list_sp_field.nums[[i]]->names
  list_seqmats[[i]]->m
  names->rownames(m)
  list_seqmats2[[i]]<-m
}
rm(i,m,names)

#apply dist.dna in each group of sequences in list_seqmats 2 -> list_mat_dist
lapply(list_seqmats2,function(i) dist.dna(as.DNAbin(i)))->list_mat_dist

unique(dflocdna$sp)->v_species

#save genetic distance as matrices:
for (i in 1:length(list_mat_dist)) {
  write.table(as.matrix(list_mat_dist[[i]]), paste0("outputs/tables/matdist/matdist.",v_species[i],".txt"), sep="\t",row.names = T,quote=F)
}
rm(i)

#save attribute tables: (dont worry too much about 'numcode' (formerly Popcode. Just a quick way to index unique pairs of coords for later))
list()->list_att_tables
for(i in 1:length(list_unique_dflocdna_sp)) {
  list_unique_dflocdna_sp[[i]]->att
  att[,-c(1,2)]->att
  rownames(att)->att[,5]
  c("Sample_ID","Latitude","Longitude","Exact_Site","numcode")->colnames(att)
  att<-att[,c(5,1:4)]
  list_att_tables[[i]]<-att
}
rm(i,att)
names(list_att_tables)<-v_species

for (i in 1:length(list_att_tables)) {
  write.table(list_att_tables[[i]], paste0("outputs/tables/att_tables/att.",v_species[i],".txt"), sep="\t",row.names = F,quote=F)
}
rm(i)

#no more needed
#rm(df_split,dfdna,dfloc,dflocdna,list_dflocdna_sp,list_seqmats,list_seqmats2,list_sp_field.nums,list_split,list_unique_dflocdna_sp,myDNA)

save.image("env_all_tables.RData")

#needed for next step:
#rm(list_att_tables,list_mat_dist,v_species)


###2: PHYLIN ANALYSIS AND MAPS ----
#If starting from here, load Fig1_env.RData and env_all_tables.RData!!

#load("env_all_tables.RData")  
load("Fig1_env.RData") #~10 sec


#library(geometry)
#library(sp)
#library(rgdal)
#library(ape)
#library(maps)
#library(mapdata)
#library(mapplots)
#library(scales)
#library(ggplot2)
#library(graphics)
#library(colorRamps)
#require(ggmap)
#require(shapefiles)
#require(rgdal)
#library(ggplot2)
#library(ggspatial)
#library(seqinr)
#library(maptools)
#library(smatr) #for major axis regression
#library(dplyr)
#library(xlsx)

#### 2.1: SET UP ENVIRONMENT AND REMOVE SPECIES WITH <2 LOCS ####
library(sp) #required for raster
library(raster) #used for most shapefiles and raster manipulation / projection
library(seqinr) #used for col2alpha
library(phylin) #used for IDW, midpoints, extract.val
library(mapplots) #used for draw.shape rivers
library(ade4) #for mantel.rtest and PCA
library(xlsx) #write.xlsx
library(dplyr) #%>%
library(MCRestimate) #to replace NAs in PCA df
library(png)

basemap<-raster("shapefiles/sa_dem_30s.bil")
as.data.frame(basemap,xy=T)->grid_tot
rm(e,basins,blank_raster,maroni_raster,poly_maroni_basin,riv14,rivs,sub,sub2)
#total=9 objects, faudra changer le nom de a

ext_basemap<-extent(-56.2,-53.0, 2.1, 5.85)#extent of basemap, formerly e
error<-c()  #Vecteur pour lister les especes non-plottés
filter1<-c()
#list()->listxyz #for raw multispecies
list()->norm_listxyz #for nomalized multispecies
list()->list_rasters

coord_tot <- data.frame()#coordinates for multispecies (same for raw and normalized)
cols<-c("darkblue","blue","deepskyblue","darkturquoise","green","yellow","orange","orangered","red","firebrick") #colors for figures
#library(seqinr) #col2alpha
sapply(cols,col2alpha,alpha=0.7)->cols_tr #make colors transparent
#Idea for alternative scale colors: 
#colfunc<-colorRampPalette(c("royalblue","springgreen","yellow","red"))
#plot(blabla,col=(colfunc(10))) #it seems by using 10 colors instead of 7 we might visualize the three "spots of divergence" better!

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
  
  #data (grid) #table containing grid centroids defining the interpolation area (7955 cells)
  grid_tot2<-grid_tot[,-3]
  
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
    
    #library(phylin)
    mp <- midpoints(coord[,1:2])
    #r.dist <- dist(coord[,1:2])
    
    mp$z <- extract.val(matdist, mp[,1:2])
    int <- idw(mp[,5], mp[,3:4], grid2) 
    label_raw<-'K2P distances'
    label_norm<-'Normalized K2P distances'
    
    #heatmap of K2P as raster
    xyz<-cbind(grid2,int)
    
    if  (length(unique(xyz$x))==1 | length(unique(xyz$y))==1 ) { #       #error2 si la grid n'a qu'un pixel
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
      
      #Faire plots avec legendes normalisées
      norm_int<-int
      
      #heatmap of normalized interpolation of K2P as raster. 
      norm_xyz<-cbind(grid2,norm_int)
      
      rasterFromXYZ(norm_xyz)->norm_rxyz
      crop <- setValues(norm_rxyz, NA)
      myshp.r <- rasterize(a, crop)
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
  
  #data (grid) #table containing grid centroids defining the interpolation area (7955 cells)
  #need Fig 1 environment with raster basemap for that
  #for some reason the line bellow gives an error now!!! Error in .local(.Object, ...) the only fix I found is using grid_tot from Fig1_Env
  #as.data.frame(basemap,xy=T)->grid_tot
  grid_tot2<-grid_tot[,-3]
  
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
    #library(ade4)
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
      myshp.r <- rasterize(a, crop)
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

#### 2.5: PLOT NORM MULTISPECIES FDR ####
###Calculate average Z for each point and map it

norm <- do.call("rbind", norm_listxyz)
#library(dplyr)
norm %>%                     #applies the following function to dfn
  group_by(x, y) %>%
  summarise(meanZ = mean(Z))->n_means #formerly just n

colnames(n_means)<-c("x","y","Z")

rasterFromXYZ(n_means)->rnorm_xyz_total
crop <- setValues(rnorm_xyz_total, NA)
myshp.r <- rasterize(a, crop)
rnorm_xyz_total.masked <- mask(x=rnorm_xyz_total, mask=myshp.r)

#with levelplot
#prepare shared rasters and layers for plots (can put at start of the script)
crop(basemap,ext_basemap)->basemap2  
myTheme <- BTCTheme()
myTheme$panel.background$col = 'dodgerblue3' #for the sea
levelplot(basemap2,col.regions=(pal(1000)),margin=F,colorkey=FALSE)->rast1 #rast1= croppped south america
layer(sp.polygons(maroni_basins))+layer(sp.lines(maroni_rivs2, col = "dodgerblue3", lwd=1))->lay1 #basins limits and rivers
extend(rnorm_xyz_total.masked,ext_basemap)->rnorm_xyz_total.masked2

dev.off()
dev.new()
#plot.new()
levelplot(rnorm_xyz_total.masked2,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label="All species",cex=1.5),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(1,1,2,1),newpage=F)
#dev.print(pdf, paste0('outputs/fdr/All_species_norm.pdf'))
#dev.print(png, paste0('outputs/fdr/All_species_norm.png'),width = 1024, height = 768)
#dev.off()

#GRADIENT
colorRampPalette(cols_tr,alpha=T)->colfunc
colfunc(1000)->cols_gr

pdf(paste0("outputs/fdr/All_species_norm_grad.pdf"))
plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main="All species")
plot(maroni_basins,add=T) #Bassins versants noir
draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
plot(rnorm_xyz_total.masked,col=cols_gr,add=T,legend=F,zlim = c(0,1))
points(coord_tot[,1:2], pch=16)
plot(rnorm_xyz_total.masked,legend.only=TRUE,col=cols_gr,zlim = c(0,1),legend.args=list(text='Mean of Z values (Normalized K2P or normalized residuals of K2P vs. real distances)', side=2, font=2, line=0.3, cex=1))
dev.off()

png(paste0("outputs/fdr/All_species_norm_grad.png"))
plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main="All species")
plot(maroni_basins,add=T) #Bassins versants noir
draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
plot(rnorm_xyz_total.masked,col=cols_gr,add=T,legend=F,zlim = c(0,1))
points(coord_tot[,1:2], pch=16)
plot(rnorm_xyz_total.masked,legend.only=TRUE,col=cols_gr,zlim = c(0,1),legend.args=list(text='Mean of Z values (Normalized K2P or normalized residuals of K2P vs. real distances)', side=2, font=2, line=0.3, cex=1))
dev.off()

print("FIFTH STEP COMPLETE: multispecies landscape")

#### 2.6: HEATMAP OF NUMBER OF SPECIMENS AND TOTAL ERROR REPORT ####

norm->count
count %>%                     
  group_by(x, y) %>%
  summarise(lengthZ = length(Z))->count

colnames(count)<-c("x","y","Z")

rasterFromXYZ(count)->count_total
crop <- setValues(count_total, NA)
myshp.r <- rasterize(a, crop)
count_total.masked <- mask(x=count_total, mask=myshp.r)

pdf(paste0("outputs/count_nb_species.pdf"))
plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main="Number of species computed by pixel")
plot(maroni_basins,add=T) #Bassins versants noir
draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
plot(count_total.masked,col=cols_tr,add=T,legend=F)
points(coord_tot[,1:2], pch=16)
plot(count_total.masked,legend.only=TRUE,col=cols,zlim = c(0,100),legend.args=list(text='Number', side=2, font=2, line=0.3, cex=1))
dev.off()

png(paste0("outputs/count_nb_species.png"))
plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main="Number of species computed by pixel")
plot(maroni_basins,add=T) #Bassins versants noir
draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
plot(count_total.masked,col=cols_tr,add=T,legend=F)
points(coord_tot[,1:2], pch=16)
plot(count_total.masked,legend.only=TRUE,col=cols,zlim = c(0,100),legend.args=list(text='Number', side=2, font=2, line=0.3, cex=1))
dev.off()


#with levelplot
#prepare shared rasters and layers for plots (can put at start of the script)
crop(basemap,ext_basemap)->basemap2  
myTheme <- BTCTheme()
myTheme$panel.background$col = 'dodgerblue3' #for the sea
levelplot(basemap2,col.regions=(pal(1000)),margin=F,colorkey=FALSE)->rast1 #rast1= croppped south america
layer(sp.polygons(maroni_basins))+layer(sp.lines(maroni_rivs2, col = "dodgerblue3", lwd=1))->lay1 #basins limits and rivers
extend(rnorm_xyz_total.masked,ext_basemap)->rnorm_xyz_total.masked2
extend(count_total.masked,ext_basemap)->count_total.masked2

dev.off()
dev.new()
#plot.new()
levelplot(rnorm_xyz_total.masked2,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label="All species",cex=1.5),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(1,1,2,1),newpage=F)

levelplot(count_total.masked2,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 100, length.out=11), colorkey = list(space='right'),margin=F,main=list(label="Number of species computed by pixel",cex=1.5),scales=list(draw=T,cex=1.5))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
r2l2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(2,1,2,1),newpage=F)
dev.print(pdf, paste0('outputs/all_sp_and_count.pdf'))
dev.print(png, paste0('outputs/all_sp_and_count.png'),width = 1024, height = 768)
dev.off()


#GRADIENT
pdf(paste0("outputs/count_nb_species_grad.pdf"))
plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main="Number of species computed by pixel")
plot(maroni_basins,add=T) #Bassins versants noir
draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
plot(count_total.masked,col=cols_gr,add=T,legend=F)
points(coord_tot[,1:2], pch=16)
plot(count_total.masked,legend.only=TRUE,col=cols_gr,zlim = c(0,100),legend.args=list(text='Number', side=2, font=2, line=0.3, cex=1))
dev.off()

png(paste0("outputs/count_nb_species_grad.png"))
plot(basemap,col=(pal(1000)), ext=ext_basemap,colNA="dodgerblue3",legend=FALSE,main="Number of species computed by pixel")
plot(maroni_basins,add=T) #Bassins versants noir
draw.shape(maroni_rivs, type = "lines", col = "dodgerblue3", lwd=1)
plot(count_total.masked,col=cols_gr,add=T,legend=F)
points(coord_tot[,1:2], pch=16)
plot(count_total.masked,legend.only=TRUE,col=cols_gr,zlim = c(0,100),legend.args=list(text='Number', side=2, font=2, line=0.3, cex=1))
dev.off()
###

#### 2.7: ERROR REPORTS, CLEAN-UP, SAVE IMAGE ####

append(error2,paste0("n=",length(error2)))->error2
error2file<-file("outputs/file_errors2.txt")
writeLines(error2, error2file)
close(error2file)

print("SIXTH STEP COMPLETE: heatmap number species")

#Clean objects from script tables
#rm(df_split,dfdna,dfloc,dflocdna,list_dflocdna_sp,list_seqmats,list_seqmats2,list_sp_field.nums,list_split,list_unique_dflocdna_sp,myDNA)
#rm(list_att_tables,list_mat_dist,v_species)

#Clean objects from landscape maps 
#rm(error,filter1,list_rasters,coord_tot,cols_tr,errorfile,seq1,seq2,zerodists,file_zero_dist,error2,matdist,grid_tot2)
#rm(minx,miny,maxx,maxy,grid2,error2file,mp,int,label_raw,label_norm,rxyz,crop,rxyz.masked,norm_int,norm_xyz,norm_rxyz,norm_rxyz.masked)
#rm(file_matdist_zero,seq3,matdists,coords,grids,mps,r.dists,index,pvalues,filter2,mantel_res,mantel_resfile,dfpval)
#rm(file_fdr_sign,file_fdr_not_sign,d.real,fit,resid,file_fdr_not_significant,norm,rnorm_xyz_total,count,count_total.masked)
#rm(basemap,coord,count_total,ext_basemap,grid_tot,list_att_tables,list_mat_dist,mantel,maroni_basins,maroni_rivs)
#rm(fdr,file_fdr_significant,i)


save.image("env_landscapes.RData")
#Needed for PCA
#rm(a,myshp.r,n_means,norm_listxyz,rnorm_xyz_total.masked,xyz,cols,v_species,pal)


######### 3: PCA --------------------------------------------------------------------------------------

#### 3.1: COMPILE DF FOR PCA (SUPER SLOW) ####

#load("env_landscapes.RData")
library(data.table)
library(MCRestimate)

#lz est la liste qui contient les 179 df avec x, y et z
norm_listxyz->lz
names(lz)<-v_species[1:length(lz)]
which(lapply(lz,class)=="NULL")->remove
lz2<-lz[-remove]
#! From now on, number IDs don't correspond anymore
lapply(lz2, function(i) paste(i$x,i$y,sep="/"))->xy #bit slow, create list of columns with each xy coordinate together in one case only
lz3 <- mapply(cbind, lz2, "xy"=xy, SIMPLIFY=F) #a bit slow ~10 sec

#TEST
#paste(lz2[[1]]$x,lz2[[1]]$y,sep="/")->xyt

#lz4<-lz3[-rm3] #REMOVE SPECIES WITH AREA <3. That-s why df and n don't have the same number of coordinates

#En utilisant la carte multispecies  comme référence, remplace les NA par la valeur correspondante du multispecies. 
#Toutes les cartes devraient du coup avoir la même couverture

#lz2=lz-toutes les df qui sont NULL -> 128 df  #! From now on, number IDs don't correspond anymore
#lz3=lz2+ "x/y" column 
#lz4=lz3 filtrée pour la taille. On ne garde que les 22 plus grandes aires
#lr2=list of 128 plotable rasters

#library(data.table)

#means=multispecies df #I think it's n_means, verifier
as.data.frame(n_means)->means #means=multispecies df
means$xy<-paste(means$x,means$y,sep="/")
t(means)->tmeans
colnames(tmeans)<-tmeans[4,]
data.table((t(as.data.frame(tmeans[3,]))))->tmeans
tmeans<-cbind(data.table(Species="All_species",tmeans)) #tmeans=line with all the multisp values

tmeans->df6  #it's called df6 because it was the 6th test
#View(df6[,1:20])

print("CREATE COMPLETE DF WITH ALL SPECIES AND MULTISP ##VERY SLOW, LIKE 1-2 HOURS")
for (i in 1:length(lz3)){
  lz3[[i]]->dfi
  t(dfi)->dfi
  colnames(dfi)<-dfi[4,]
  data.table(t(as.data.frame(dfi[3,])))->dfi
  dfi<-cbind(data.table(Species=names(lz3[i])),dfi)
  merge(data.table(df6),data.table(dfi),all=T,by=intersect(colnames(df6),colnames(dfi)))->df6
  print(i)
}

print("TOTAL DF FOR PCA WITH NAs COMPLETED")
save.image("env_PCAtable.RData")

#### 3.2: FORMAT THE DF FOR PCA ----

#load("env_PCAtable.RData")

#View(df6[,1:20])
as.data.frame(df6)->df6 #ça deconne de ouf avec data.table
rownames(df6)<-df6$Species
df6[,-1]->df6
#View(df6[,1:20])
#df=complete DF with maximum extant, all species and multispecies (=128 rows).

#ok, ya plus qu'a remplacer les nas
#library(MCRestimate)
replace.NA(df6, as.numeric(df6[5,]), byRow = F)->df6_2
print("replace all NAs by multisp value, takes ~10-15min")
#View(df6[,1:20])
#View(df6_2[,1:20])

df6_3 <- as.data.frame(sapply(df6_2, as.numeric)) #MUST convert everything in numeric, must be df6_2 or we lose the rownames
rownames(df6_3)<-rownames(df6_2)
#View(df6_3[,1:20])

#I put multispecies in last row to plot it last
rbind(df6_3[1:4,],df6_3[6:129,],df6_3[5,])->df6_4
View(df6_4[,1:20])

print("8TH STEP: REPLACE NAS COMPLETED")

#### 3.3: SAVE RASTERS AS PICS (optional) ####

#PAS BESOIN DE REFAIRE CHAQUE FOIS/
#Save rasters to display on PCA #par(bg=NA) is for transparence (=no white background)
names(list_rasters)<-v_species[1:length(list_rasters)] #pas sûr que c'est necessaire si j'ai bien corrigé le script
#Remove rasters with null values
which(lapply(list_rasters,class)=="NULL")->remove
lr2<-list_rasters[-remove]
par(mfrow=c(3,5))
###Just taking a look
#for(i in 1:15){
#  plot(lr2[[i]])
#}
#par(mfrow=c(3,5))
#for(i in 16:30){
#  plot(lr2[[i]])
#}

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
#/PAS BESOIN DE REFAIRE CHAQUE FOIS

print("9th STEP: LANDSCAPE RASTERS SAVED")

#### 3.4: CROP DF VALUES TO MARONI BASIN BOUNDARIES ####

print("10th STEP: CROP DF VALUES TO MARONI BASIN")

#df6_4: complete df with all species and multisp, NAs replaced by multisp values, multisp line last in df for better plotting

df6_4->df8
df8->df18
#viou(df18)  #those numbers correspond to the different trials I did

#1: add values of x and y for each pixel

as.numeric(gsub('/(.*)',"", colnames(df18)))->x
as.numeric(gsub(".*/", "", colnames(df18)))->y
rbind(df18,x=x)->df18_2
rbind(df18_2,y=y)->df18_2

#max(y)=5.670833=NORTH
#min(y)=2.670833=SOUTH
#min(x)=-55.4375=WEST
#max(x)=-53.1625=EAST

#To scan the basin "by row, from left to right":
#We want to sort first by y in descending order then by x in ascending order

df18_3<-df18_2[,order(-df18_2[131,],df18_2[130,])]
#viou(df18_3)

#TESTS VISULALISATION:
as.matrix(df18_3)->m_3
#viou(m_3)
#scatter_fill(m_3[130,],m_3[131,],m_3[4,],nlevels=10,main="TEST",pch=".",cex=8)

#Relou mais on va reprendre chaque espece et les transformer en raster pour pouvoir les cropper
#puis les remettre ensemble. Le probleme venait de la résolution je crois

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
#viou(df18_4)

cbind(df18_4,l3[[1]][,1],l3[[1]][,2])->df18_5
colnames(df18_5)<-c(names(l),"x","y")
#View(df18_5)

df18_6<-as.data.frame(t(df18_5))
#viou(df18_6)

as.matrix(df18_6)->m_6

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
#### 3.6 USE SCATTER FILL FUNCTION #####

scatter_fill(m_6[130,],m_6[131,],m_6[129,],nlevels=10,main="TEST",pch=".",cex=8)

df18_6->df18

print("DATA FRAME READY FOR PCA")
save.image("env_PCAtable_final.RData")

#### 3.7: DO THE PCAs #####
#load("env_PCAtable_final.RData")
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

s.label(pca18$li,1,2)
png(paste0("outputs/PCA/","labels.png"))
s.label(pca18$li,1,2)
dev.off()

summary(pca18)
sink("outputs/PCA/summary.txt")
print(summary(pca18))
sink()

library(png)
lpix2 <- list()
for (nomfic in list.files("outputs/PNGS_129/", pattern = ".png")) {nomobj <- strsplit(nomfic, "[.]")[[1]][1]
toto <- readPNG(paste("outputs/PNGS_129/", nomfic, sep = ""))
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
#text(pca18$co[,1:2],labels=colnames(df8))

#### Latitude coordinates as variables ####
png(paste0("outputs/PCA/","projections_latitude12.png"))
scatter_fill(pca18$co[,1],pca18$co[,2],as.numeric(df18[131,]),main="PCA PROJECTION OF LATITUDES (AXES 1-2)",pch=".",cex=3)
dev.off()
pdf(paste0("outputs/PCA/","projections_latitude12.pdf"))
scatter_fill(pca18$co[,1],pca18$co[,2],as.numeric(df18[131,]),main="PCA PROJECTION OF LATITUDES (AXES 1-2)",pch=".",cex=3)
dev.off()

png(paste0("outputs/PCA/","legend_latitude.png"))
scatter_fill(as.numeric(df18[130,]),as.numeric(df18[131,]),as.numeric(df18[131,]),main="Colors latitudes basin",pch=".",cex=3)
dev.off()
pdf(paste0("outputs/PCA/","legend_latitude.pdf"))
scatter_fill(as.numeric(df18[130,]),as.numeric(df18[131,]),as.numeric(df18[131,]),main="Colors latitudes basin",pch=".",cex=3)
dev.off()

png(paste0("outputs/PCA/","projections_latitude23.png"))
scatter_fill(pca18$co[,2],pca18$co[,3],as.numeric(df18[131,]),main="PCA PROJECTION OF LATITUDES (AXES 2-3)",pch=".",cex=3)
dev.off()
pdf(paste0("outputs/PCA/","projections_latitude23.pdf"))
scatter_fill(pca18$co[,2],pca18$co[,3],as.numeric(df18[131,]),main="PCA PROJECTION OF LATITUDES (AXES 2-3)",pch=".",cex=3)
dev.off()

#### West_North_East ####
#NEED FOLDER OLD_MAP_PICS
dir.create("outputs/PCA/coverage/")

list.files("./OLD_MAPS_PICS/West_North_East/")->West_North_East
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

list.files("./OLD_MAPS_PICS/West_East/")->West_East

unlist(lapply(strsplit(West_East, "norm_"), '[[', 2))->West_East
unlist(lapply(strsplit(West_East, "-1"), '[[', 1))->West_East

c()->West_East_nb
for(i in 1:length(West_East)){
  West_East_nb<-c(West_East_nb,which(rownames(df8)==West_East[i]))
}

s.label(pca18$li[West_East_nb,2:3])
png(paste0("outputs/PCA/coverage/","upper_maroni_labels23.png"))
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

list.files("./OLD_MAPS_PICS/East/")->East

unlist(lapply(strsplit(East, "norm_"), '[[', 2))->East
unlist(lapply(strsplit(East, "-1"), '[[', 1))->East

c()->East_nb
for(i in 1:length(East)){
  East_nb<-c(East_nb,which(rownames(df8)==East[i]))
}

s.label(pca18$li[East_nb,2:3])
png(paste0("outputs/PCA/coverage/","east_maroni_label23.png"))
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

list.files("./OLD_MAPS_PICS/North_East/")->north_east

unlist(lapply(strsplit(north_east, "norm_"), '[[', 2))->north_east
unlist(lapply(strsplit(north_east, "-1"), '[[', 1))->north_east

c()->north_east_nb
for(i in 1:length(north_east)){
  north_east_nb<-c(north_east_nb,which(rownames(df8)==north_east[i]))
}

s.label(pca18$li[north_east_nb,2:3])
png(paste0("outputs/PCA/coverage/","north_east_maroni_label23.png"))
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

list.files("./OLD_MAPS_PICS/Others/")->others

unlist(lapply(strsplit(others, "norm_"), '[[', 2))->others
unlist(lapply(strsplit(others, "-1"), '[[', 1))->others

c()->others_nb
for(i in 1:length(others)){
  others_nb<-c(others_nb,which(rownames(df8)==others[i]))
}

s.label(pca18$li[others_nb,2:3])
png(paste0("outputs/PCA/coverage/","others_maroni_label23.png"))
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

#### West_North_East + West_East ####
df18we<-df18[,order(df18[130,],-df18[131,])]
#viou(df18we)
#low values=west, high=east

pca18we<-dudi.pca(df18we[1:129,], center=T, scale=F, scannf = F, nf=5) #~5 seconds
#do not do scatter(pca18we,1,2), too many variables to plot!
par(mfrow=c(1,1))
#scatter(pca18we,1,2)
barplot(pca18we$eig[1:6])
s.label(pca18we$li,1,2)
summary(pca18we)

library(png)
lpix2 <- list()
for (nomfic in list.files("outputs/PNGS_129/", pattern = ".png")) {nomobj <- strsplit(nomfic, "[.]")[[1]][1]
toto <- readPNG(paste("outputs/PNGS_129/", nomfic, sep = ""))
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

plot(pca18we$li[,1:2], t="n",axes=TRUE, xlim=(c(min((pca18we$li[,1]))-10,max((pca18we$li[,1]))+10)), ylim=(c(min(pca18we$li[,2])-10, max(pca18we$li[,2])+10)), xlab="PC1", ylab="PC2")
abline(h = 0, v = 0)
thumbnails((pca18we$li[,1]), pca18we$li[,2], lpix2, width = 0.03*diff(range(pca18we$li[,1])),height = 0.07*diff(range(pca18we$li[,2])))

plot(pca18we$li[,2:3], t="n",axes=TRUE, xlim=(c(min(pca18we$li[,2])-10,max(pca18we$li[,2])+10)), ylim=(c(min(pca18we$li[,3])-10, max(pca18we$li[,3])+10)), xlab="PC2", ylab="PC3")
abline(h = 0, v = 0)
thumbnails(pca18we$li[,2], pca18we$li[,3], lpix2, width = 0.03*diff(range(pca18we$li[,1])),height = 0.06*diff(range(pca18we$li[,2])))

plot(pca18we$li[,4:5], t="n",axes=TRUE, xlim=(c(min(pca18we$li[,4])-10,max(pca18we$li[,4])+10)), ylim=(c(min(pca18we$li[,5])-10, max(pca18we$li[,5])+10)), xlab="PC4", ylab= "Axis 5")
abline(h = 0, v = 0)
thumbnails(pca18we$li[,4], pca18we$li[,5], lpix2, width = 0.03*diff(range(pca18we$li[,1])),height = 0.06*diff(range(pca18we$li[,2])))

#Variables?
#plot(pca18we$co[,1:2], t="p",axes=TRUE, xlim=(c(min(pca18we$co[,1])-0.01,max(pca18we$co[,1])+0.01)), ylim=(c(min(pca18we$co[,2])-0.01, max(pca18we$co[,2])+0.01)), xlab="PC1", ylab="PC2")
#abline(h = 0, v = 0)
#text(pca18we$co[,1:2],labels=colnames(df8))

#### Longitude coordinates as variables ####
png(paste0("outputs/PCA/","projections_longitude12.png"))
scatter_fill((pca18we$co[,1]),pca18we$co[,2],as.numeric(df18[130,]),main="PCA PROJECTIONS OF LONGITUDES (AXES1-2)",pch=".",cex=3)
dev.off()
pdf(paste0("outputs/PCA/","projections_longitude12.pdf"))
scatter_fill((pca18we$co[,1]),pca18we$co[,2],as.numeric(df18[130,]),main="PCA PROJECTIONS OF LONGITUDES (AXES1-2)",pch=".",cex=3)
dev.off()

png(paste0("outputs/PCA/","legend_longitude.png"))
scatter_fill(as.numeric(df18we[130,]),as.numeric(df18we[131,]),as.numeric(df18we[130,]),main="Colors longitudes basin",pch=".",cex=3)
dev.off()
pdf(paste0("outputs/PCA/","legend_longitude.pdf"))
scatter_fill(as.numeric(df18we[130,]),as.numeric(df18we[131,]),as.numeric(df18we[130,]),main="Colors longitudes basin",pch=".",cex=3)
dev.off()

png(paste0("outputs/PCA/","projections_longitude23.png"))
scatter_fill(pca18we$co[,2],pca18we$co[,3],as.numeric(df18we[130,]),nlevels=1000,main="PCA PROJECTIONS OF LONGITUDES (AXES2-3)",pch=".",cex=3)
dev.off()
pdf(paste0("outputs/PCA/","projections_longitude23.pdf"))
scatter_fill(pca18we$co[,2],pca18we$co[,3],as.numeric(df18we[130,]),nlevels=1000,main="PCA PROJECTIONS OF LONGITUDES (AXES2-3)",pch=".",cex=3)
dev.off()

save.image("env_PCA_final.RData")

##### 4: K-MEAN CLUSTER #####
#load("env_PCA_final.RData")

library(adegenet)
dir.create('outputs/kmean_clusters')

#These functions implement the clustering procedure used in Discriminant Analysis of Principal Components 
#(DAPC, Jombart et al. 2010). This procedure consists in running successive K-means with an increasing number 
#of clusters (k), after transforming data using a principal component analysis (PCA). 
#For each model, a statistical measure of goodness of fit (by default, BIC) is computed, 
#which allows to choose the optimal k


#plot(pca18$li[,1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
#abline(h = 0, v = 0)
#thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.07*diff(range(pca18$li[,2])))

#### 4.1 FIND.CLUSTERS ####
#### Test K=5 clusters ####
find.clusters(df18[1:129,],n.pca=3,max.n.clust=20,dudi=pca18,n.iter=1e7,n.clust=5)->clusters_5
clusters_5

dir.create('outputs/kmean_clusters/clusters_5')
sink("outputs/kmean_clusters/clusters_5/clusters_5.txt")
print(clusters_5)
sink()

for (i in 1:(length(clusters_5$size))) {
  which(clusters_5$grp==i)
  names(clusters_5$grp[which(clusters_5$grp==i)])
  
  plot(pca18$li[which(clusters_5$grp==i),1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_5$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.07*diff(range(pca18$li[,2])))
  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_5/grp',i,'_12.pdf'))
  dev.print(png, paste0('./outputs/kmean_clusters/clusters_5/grp',i,'_12.png'),width = 1024, height = 768)
  
  
}

for (i in 1:(length(clusters_5$size))) {
  which(clusters_5$grp==i)
  names(clusters_5$grp[which(clusters_5$grp==i)])
  
  plot(pca18$li[which(clusters_5$grp==i),2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[,2])-10,max(pca18$li[,2])+10)), ylim=(c(min(pca18$li[,3])-10, max(pca18$li[,3])+10)), xlab="PC1", ylab="PC2")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_5$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.06*diff(range(pca18$li[,2])))
  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_5/grp',i,'_23.pdf'))
  dev.print(png, paste0('./outputs/kmean_clusters/clusters_5/grp',i,'_23.png'),width = 1024, height = 768)
  
  
}



#### Test K=17 clusters ####
find.clusters(df18[1:129,],n.pca=3,max.n.clust=20,dudi=pca18,n.iter=1e7,n.clust=17)->clusters_17
clusters_17

dir.create('outputs/kmean_clusters/clusters_17')
sink("outputs/kmean_clusters/clusters_17/clusters_17.txt")
print(clusters_17)
sink()

for (i in 1:(length(clusters_17$size))) {
  which(clusters_17$grp==i)
  names(clusters_17$grp[which(clusters_17$grp==i)])
  
  plot(pca18$li[which(clusters_17$grp==i),1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_17$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.07*diff(range(pca18$li[,2])))
  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_17/grp',i,'_12.pdf'))
  dev.print(png, paste0('./outputs/kmean_clusters/clusters_17/grp',i,'_12.png'),width = 1024, height = 768)
  
  
}

for (i in 1:(length(clusters_17$size))) {
  which(clusters_17$grp==i)
  names(clusters_17$grp[which(clusters_17$grp==i)])
  
  plot(pca18$li[which(clusters_17$grp==i),2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[,2])-10,max(pca18$li[,2])+10)), ylim=(c(min(pca18$li[,3])-10, max(pca18$li[,3])+10)), xlab="PC1", ylab="PC2")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_17$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.06*diff(range(pca18$li[,2])))
  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_17/grp',i,'_23.pdf'))
  dev.print(png, paste0('./outputs/kmean_clusters/clusters_17/grp',i,'_23.png'),width = 1024, height = 768)
  
  
}

#### K=11 clusters ####
find.clusters(df18[1:129,],n.pca=3,max.n.clust=20,dudi=pca18,n.iter=1e7,n.clust=11)->clusters_11
clusters_11

dir.create('outputs/kmean_clusters/clusters_11b')
sink("outputs/kmean_clusters/clusters_11b/clusters_11.txt")
print(clusters_11)
sink()

for (i in 1:(length(clusters_11$size))) {
  which(clusters_11$grp==i)
  names(clusters_11$grp[which(clusters_11$grp==i)])
  
  plot(pca18$li[which(clusters_11$grp==i),1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_11$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.07*diff(range(pca18$li[,2])))
  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_12.pdf'))
  dev.print(png, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_12.png'),width = 1024, height = 768)
  
  
}

for (i in 1:(length(clusters_11$size))) {
  which(clusters_11$grp==i)
  names(clusters_11$grp[which(clusters_11$grp==i)])
  
  plot(pca18$li[which(clusters_11$grp==i),2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[,2])-10,max(pca18$li[,2])+10)), ylim=(c(min(pca18$li[,3])-10, max(pca18$li[,3])+10)), xlab="PC2", ylab="PC3")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_11$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.06*diff(range(pca18$li[,2])))
  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_23.pdf'))
  dev.print(png, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_23.png'),width = 1024, height = 768)
  
  
}

save.image('env_K_clusters.RData')


##### 4.2 PCA PROJECTIONS GROUPED BY K=11 CLUSTERS #####

par (mfrow=c(2,3))
for (i in 1:(length(clusters_11$size))) {
  which(clusters_11$grp==i)
  names(clusters_11$grp[which(clusters_11$grp==i)])
  
  plot(pca18$li[which(clusters_11$grp==i),1:2],main=paste0('K=',i), t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_11$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.07*diff(range(pca18$li[,2])))
  #dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_12.pdf'))
  #dev.print(png, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_12.png'),width = 1024, height = 768)
  
  
}


for (i in 1:(length(clusters_11$size))) {
  which(clusters_11$grp==i)
  names(clusters_11$grp[which(clusters_11$grp==i)])
  
  plot(pca18$li[which(clusters_11$grp==i),2:3],main=paste0('K=',i),t="n",axes=TRUE, xlim=(c(min(pca18$li[,2])-10,max(pca18$li[,2])+10)), ylim=(c(min(pca18$li[,3])-10, max(pca18$li[,3])+10)), xlab="PC2", ylab="PC3")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_11$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.06*diff(range(pca18$li[,2])))
  #dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_23.pdf'))
  #dev.print(png, paste0('./outputs/kmean_clusters/clusters_11b/grp',i,'_23.png'),width = 1024, height = 768)
  
  
}

save.image('env_K_clusters.RData')

#### 4.3 PLOT 50 MAPS / FIGURE FOR K=11 ####
#SAME SCRIPT FOR K=17 CAN BE FOUND IN E:\BARCODE\TRAVAIL V4\1905 MAJ FINALE\SCRIPTS\tests_figures\k_clusters_par_50.R

#setwd("E:/BARCODE/TRAVAIL V4/1905 MAJ FINALE/LANDSCAPE FINAL SCRIPT")
#load('env_K_clusters.RData')
#load('heavy scripts/maroni_rivs2.RData')
library(rasterVis)
library(pdftools)
#Inkscape ne peut pas ouvrir de pdfs plus gros que 2Mb sur mon laptop (???) et les png ne se
#sauvegardent pas bien dans le script, donc on test convertir pdf en png avec pdftools


crop(basemap,ext_basemap)->basemap2
myTheme <- BTCTheme()
myTheme$panel.background$col = 'dodgerblue3' #for the sea
dir.create("./outputs/kmean_clusters5")

levelplot(basemap2,col.regions=(pal(1000)),margin=F,colorkey=FALSE)->rast1
layer(sp.polygons(maroni_basins))->lay1
#layer(sp.polygons(maroni_basins))+layer(sp.lines(maroni_rivs2, col = "dodgerblue3", lwd=1))->lay1 #lourd avec les rivieres mais peut etre refaire plus tard

v_species[200]<-"All_species"
rnorm_xyz_total.masked->lr2[[129]]  ###NEED TO PUT ALL SPECIES RASTER IN LR2
names(lr2)<-c(names(lr2)[1:128],"All_species")
list_att_tables[[200]]<-coord_tot[,c(2,1,2,1)]


#### Figure test ####
dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10
grps=4

for (grp in c(1,4)){
  #par(pin=c(5, 4))
  for (i in 1:length(which(clusters_11$grp==grp))){
    which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
    which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
    print(coordnb)
    lr2[[spnb]]->landscape
    extend(landscape,ext_basemap)->landscape
    #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=1,y=grid::unit(4, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=gsub("_","\n",names(lr2[spnb])),cex=0.8,y=grid::unit(-0.1, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.8))->r2l2
    r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
    #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
    print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
    #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
    t=t+1
  }}
#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters5/figureTEST.pdf'))
dev.print(png, paste0('./outputs/kmean_clusters5/figureTEST.png'),width=width, height = height, unit = "cm",res=300)



#### Figure 1 ####
dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10
#Groups 1-5
for (grp in c(1:5)){
  #par(pin=c(5, 4))
  for (i in 1:length(which(clusters_11$grp==grp))){
    which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
    which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
    lr2[[spnb]]->landscape
    extend(landscape,ext_basemap)->landscape
    gsub("_"," ",names(lr2[spnb]))->name
    #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
    r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
    #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
    print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
    #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
    t=t+1
  }}

#Group 6
grp=6
#par(pin=c(5, 4))
for (i in 1){
  which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
  which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
  lr2[[spnb]]->landscape
  extend(landscape,ext_basemap)->landscape
  gsub("_"," ",names(lr2[spnb]))->name
  #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
  r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
  #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
  print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
  #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
  t=t+1
}

#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters5/figure1.pdf'))
dev.print(svg, paste0('./outputs/kmean_clusters5/figure1.svg'))
#dev.print(png, paste0('./outputs/kmean_clusters5/figure1.png'),width=width+1, height = height+1, unit = "cm",res=300) # marche pas for some reason
pdf_convert('./outputs/kmean_clusters5/figure1.pdf',format = "png",filenames='./outputs/kmean_clusters5/figure1.png',dpi=300)

#### Figure 2 ####

dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10

#Group 6 (2-40)
grp=6
#par(pin=c(5, 4))
for (i in 2:40){
  which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
  which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
  lr2[[spnb]]->landscape
  extend(landscape,ext_basemap)->landscape
  gsub("_"," ",names(lr2[spnb]))->name
  #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
  r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
  #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
  print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
  #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
  t=t+1
}


#Groups 7-8
for (grp in c(7:8)){
  #par(pin=c(5, 4))
  for (i in 1:length(which(clusters_11$grp==grp))){
    which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
    which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
    lr2[[spnb]]->landscape
    extend(landscape,ext_basemap)->landscape
    gsub("_"," ",names(lr2[spnb]))->name
    #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
    r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
    #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
    print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
    #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
    t=t+1
  }}

#Group 13 (1:22)
grp=9
#par(pin=c(5, 4))
for (i in 1:3){
  which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
  which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
  lr2[[spnb]]->landscape
  extend(landscape,ext_basemap)->landscape
  gsub("_"," ",names(lr2[spnb]))->name
  #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
  r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
  #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
  print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
  #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
  t=t+1
}

#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters5/figure2.pdf'))
dev.print(svg, paste0('./outputs/kmean_clusters5/figure2.svg'))
#dev.print(png, paste0('./outputs/kmean_clusters5/figure2.png'),width=width+1, height = height+1, unit = "cm",res=300) # marche pas for some reason
pdf_convert('./outputs/kmean_clusters5/figure2.pdf',format = "png",filenames='./outputs/kmean_clusters5/figure2.png',dpi=300)


#### Figure 3 ####

dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10

#Group 10 (4:13)
grp=9
#par(pin=c(5, 4))
for (i in 4:13){
  which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
  which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
  lr2[[spnb]]->landscape
  extend(landscape,ext_basemap)->landscape
  gsub("_"," ",names(lr2[spnb]))->name
  


  #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
  r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
  #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
  print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
  #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
  t=t+1
}


#Groups 10-11
for (grp in c(10:11)){
  #par(pin=c(5, 4))
  for (i in 1:length(which(clusters_11$grp==grp))){
    which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
    which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
    lr2[[spnb]]->landscape
    extend(landscape,ext_basemap)->landscape
    gsub("_"," ",names(lr2[spnb]))->name
    #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
    r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
    #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
    print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
    #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
    t=t+1
  }}


#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters5/figure3.pdf'))
dev.print(svg, paste0('./outputs/kmean_clusters5/figure3.svg'))
dev.print(png, paste0('./outputs/kmean_clusters5/figure3.png'),width=width+1, height = height+1, unit = "cm",res=300) # marche pas for some reason
pdf_convert('./outputs/kmean_clusters5/figure3.pdf',format = "png",filenames='./outputs/kmean_clusters5/figure3.png',dpi=300)


#### 4.4 FIND.CLUSTERS USING AXES 2 AND 3 ONLY ####

#### FIND.CLUSTERS.NOAXIS1 ####
source("find.clusters.noaxis1.r")
#### TEST ####
#find.clusters.noaxis1(df18[1:129,],n.pca=3,max.n.clust=20,dudi=pca18,n.iter=1e7,n.start=100)
###DROPs IN BIC VALUES SEEM TO INDICATE A OPTIMAL K BETWEEN 11 AND 16


#### K=11 clusters ####
find.clusters.noaxis1(df18[1:129,],n.pca=3,max.n.clust=20,dudi=pca18,n.iter=1e7,n.clust=11,n.start=100)->clusters_11
clusters_11

dir.create('outputs/kmean_clusters/clusters_11_a23')
sink("outputs/kmean_clusters/clusters_11_a23/clusters_11.txt")
print(clusters_11)
sink()


#Axes 1 & 2, not relevant in that case
#for (i in 1:(length(clusters_11$size))) {
 # which(clusters_11$grp==i)
  #names(clusters_11$grp[which(clusters_11$grp==i)])
  #
#  plot(pca18$li[which(clusters_11$grp==i),1:2], t="n",axes=TRUE, xlim=(c(min(pca18$li[,1])-10,max(pca18$li[,1])+10)), ylim=(c(min(pca18$li[,2])-10, max(pca18$li[,2])+10)), xlab="PC1", ylab="PC2")
#  
#  abline(h = 0, v = 0)
#  thumbnails <- function(x, y, images, width = 2.8,
#                         height = 4.8){
#    for (ii in which(clusters_11$grp==i)){
#      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
#                  ybottom= y[ii] - 1*height,
#                  xright=x[ii] + 1*width,
#                  ytop= y[ii] + 1*height, interpolate=FALSE)
#    }
#  }
#  
#  
#  thumbnails(pca18$li[,1], pca18$li[,2], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.07*diff(range(pca18$li[,2])))
#  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_11_a23/grp',i,'_12.pdf'))
#  dev.print(png, paste0('./outputs/kmean_clusters/clusters_11_a23/grp',i,'_12.png'),width = 1024, height = 768)
#  
#  
#}
par(mfrow=c(2,3)) ### DELETE THIS LINE TO PLOT SEPARATE CLEANS PLOTS, OR MODIFY DO MAKE THEM CLEANER
for (i in 1:(length(clusters_11$size))) {
  which(clusters_11$grp==i)
  names(clusters_11$grp[which(clusters_11$grp==i)])
  
  plot(pca18$li[which(clusters_11$grp==i),2:3], t="n",axes=TRUE, xlim=(c(min(pca18$li[,2])-10,max(pca18$li[,2])+10)), ylim=(c(min(pca18$li[,3])-10, max(pca18$li[,3])+10)), xlab="PC2", ylab="PC3")
  
  abline(h = 0, v = 0)
  thumbnails <- function(x, y, images, width = 2.8,
                         height = 4.8){
    for (ii in which(clusters_11$grp==i)){
      rasterImage(lpix2[[ii]], xleft=x[ii] - 1*width,
                  ybottom= y[ii] - 1*height,
                  xright=x[ii] + 1*width,
                  ytop= y[ii] + 1*height, interpolate=FALSE)
    }
  }
  
  
  thumbnails(pca18$li[,2], pca18$li[,3], lpix2, width = 0.03*diff(range(pca18$li[,1])),height = 0.06*diff(range(pca18$li[,2])))
  dev.print(pdf, paste0('./outputs/kmean_clusters/clusters_11_a23/grp',i,'_23.pdf'))
  dev.print(png, paste0('./outputs/kmean_clusters/clusters_11_a23/grp',i,'_23.png'),width = 1024, height = 768)
  
  
}

save.image('env_K_clusters_23.RData')

#### 4.5 PLOT 50 MAPS / FIGURE FOR K=11 ####
#SAME SCRIPT FOR K=17 CAN BE FOUND IN E:\BARCODE\TRAVAIL V4\1905 MAJ FINALE\SCRIPTS\tests_figures\k_clusters_par_50.R

#setwd("E:/BARCODE/TRAVAIL V4/1905 MAJ FINALE/LANDSCAPE FINAL SCRIPT")
#load('env_K_clusters_23.RData')
#load('heavy scripts/maroni_rivs2.RData')
library(rasterVis)
library(pdftools)
#Inkscape ne peut pas ouvrir de pdfs plus gros que 2Mb sur mon laptop (???) et les png ne se
#sauvegardent pas bien dans le script, donc on test convertir pdf en png avec pdftools


crop(basemap,ext_basemap)->basemap2
myTheme <- BTCTheme()
myTheme$panel.background$col = 'dodgerblue3' #for the sea
dir.create("./outputs/kmean_clusters_a23")

levelplot(basemap2,col.regions=(pal(1000)),margin=F,colorkey=FALSE)->rast1
layer(sp.polygons(maroni_basins))->lay1
#layer(sp.polygons(maroni_basins))+layer(sp.lines(maroni_rivs2, col = "dodgerblue3", lwd=1))->lay1 #lourd avec les rivieres mais peut etre refaire plus tard

v_species[200]<-"All_species"
rnorm_xyz_total.masked->lr2[[129]]  ###NEED TO PUT ALL SPECIES RASTER IN LR2
names(lr2)<-c(names(lr2)[1:128],"All_species")
list_att_tables[[200]]<-coord_tot[,c(2,1,2,1)]


#### Figure test ####
dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10
grps=4

for (grp in c(1,4)){
  #par(pin=c(5, 4))
  for (i in 1:length(which(clusters_11$grp==grp))){
    which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
    which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
    print(coordnb)
    lr2[[spnb]]->landscape
    extend(landscape,ext_basemap)->landscape
    #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=1,y=grid::unit(4, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=gsub("_","\n",names(lr2[spnb])),cex=0.8,y=grid::unit(-0.1, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.8))->r2l2
    r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
    #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
    print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
    #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
    t=t+1
  }}
#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters_a23/figureTEST.pdf'))
dev.print(png, paste0('./outputs/kmean_clusters_a23/figureTEST.png'),width=width, height = height, unit = "cm",res=300)



#### Figure 1 ####
dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10
#Groups 1-3
for (grp in c(1:3)){
  #par(pin=c(5, 4))
  for (i in 1:length(which(clusters_11$grp==grp))){
    which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
    which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
    lr2[[spnb]]->landscape
    extend(landscape,ext_basemap)->landscape
    gsub("_"," ",names(lr2[spnb]))->name
    #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
    r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
    #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
    print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
    #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
    t=t+1
  }}


#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters_a23/figure1.pdf'))
dev.print(svg, paste0('./outputs/kmean_clusters_a23/figure1.svg'))
#dev.print(png, paste0('./outputs/kmean_clusters_a23/figure1.png'),width=width+1, height = height+1, unit = "cm",res=300) # marche pas for some reason
pdf_convert('./outputs/kmean_clusters_a23/figure1.pdf',format = "png",filenames='./outputs/kmean_clusters_a23/figure1.png',dpi=300)

#### Figure 2 ####

dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10


#Groups 4-11
for (grp in c(4:11)){
  #par(pin=c(5, 4))
  for (i in 1:length(which(clusters_11$grp==grp))){
    which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
    which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
    lr2[[spnb]]->landscape
    extend(landscape,ext_basemap)->landscape
    gsub("_"," ",names(lr2[spnb]))->name
    #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
    rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
    r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
    #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
    print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
    #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
    t=t+1
  }}

#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters_a23/figure2.pdf'))
dev.print(svg, paste0('./outputs/kmean_clusters_a23/figure2.svg'))
#dev.print(png, paste0('./outputs/kmean_clusters_a23/figure2.png'),width=width+1, height = height+1, unit = "cm",res=300) # marche pas for some reason
pdf_convert('./outputs/kmean_clusters_a23/figure2.pdf',format = "png",filenames='./outputs/kmean_clusters_a23/figure2.png',dpi=300)


#### Figure 3 ####

dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10

#Group 11 (4:32)
grp=11
#par(pin=c(5, 4))
for (i in 4:32){
  which((names(lr2))==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->spnb
  which(v_species==(names(clusters_11$grp[which(clusters_11$grp==grp)])[i]))->coordnb
  lr2[[spnb]]->landscape
  extend(landscape,ext_basemap)->landscape
  gsub("_"," ",names(lr2[spnb]))->name
  
  
  
  #levelplot(landscape,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = list(space='right'),margin=F,main=list(label=gsub("_"," ",names(lr2[spnb]))))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  levelplot(landscape,par.settings = myTheme,sub=list(label=paste0("K=",grp),cex=0.8,y=grid::unit(4.5, "mm")),col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=sub(" ","\n",name),cex=0.8,y=grid::unit(-0.01, "mm")),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
  rast2+layer(sp.points(SpatialPoints(list_att_tables[[coordnb]][,c(4,3)]),pch=21,col='black',fill="white",cex=0.5))->r2l2
  r2l2+as.layer(rast1+lay1,under=T)->plot  #wtf? corriger?
  #rast2+as.layer(rast1+lay1,under=T)+lay2->plot
  print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE,panel.width = list(2.25, "cm"),panel.height = list(4.5, "cm"))
  #print(plot,split=c(xs[t], ys[t], nx, ny),newpage=FALSE)
  t=t+1
}


#title(paste0("K=",grp), outer=TRUE,line = -5,cex.main=3,adj=0)
dev.print(pdf, paste0('./outputs/kmean_clusters_a23/figure3.pdf'))
dev.print(svg, paste0('./outputs/kmean_clusters_a23/figure3.svg'))
dev.print(png, paste0('./outputs/kmean_clusters_a23/figure3.png'),width=width+1, height = height+1, unit = "cm",res=300) # marche pas for some reason
pdf_convert('./outputs/kmean_clusters_a23/figure3.pdf',format = "png",filenames='./outputs/kmean_clusters_a23/figure3.png',dpi=300)



#### TEST ####
#LDA en utilisant les groupes k-mean comme facteur 



#### TEST DIRECT CLUSTERING ####
#distances euclidiennes entre tes coordonnées d'individus selon les axes sélectionnés 
#(2-3 par exemple)-> Clustering (voir si les deux approches concordent)
library(ade4)
library(pvclust)

#coord_inds_lda<- pca18$li[,2:3]

#Euclidean distances between individuals on axes 2 and 3 of PCA
#dist(coord_inds_23)->dist_coord_ints_23
#hclust(dist_coord_ints_23,method="mcquitty")->hclust_a23
#plot(hclust_a23)

#as.data.frame(t(coord_inds_23))->tcoord_inds_23
#pvclust(tcoord_inds_23,method.hclust="mcquitty",method.dist="euclidean",nboot=10000,r=1)->pvclust_a23
#if I dont put r=1, I get an error
#see https://www.biostars.org/p/95477/
#plot(pvclust_a23)

setwd("E:/BARCODE/190924 IN PROGRESS")
#load("E:/BARCODE/190924 IN PROGRESS/.RData")

coord_inds<- pca18$li[, 2:3]
cors_coph<-c()
as.matrix(dist.quant (coord_inds, method=1,diag=T,upper=T))->dist_coord_ints_23
methods=c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median","centroid")
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
for (m in methods) {
  pvclust(dist_coord_ints_23,method.hclust=m,method.dist="euclidean",nboot=1000)->pvclust_a23
  dev.off()
  width = 60
  height = 30
  dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
  plot(pvclust_a23)
  pvrect(pvclust_a23,alpha=0.95)
  dev.print(pdf, paste0('./pvclust/a23_', m,'.pdf'))
  dev.print(svg, paste0('./pvclust/a23_', m,'.svg'))
  dev.print(png, paste0('./pvclust/a23_', m,'.png'),width=width+1, height = height+1, unit = "cm",res=300) 
  
  #based on https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre
  cophenetic(pvclust_a23$hclust)->coph
  cor(as.dist(dist_coord_ints_23),coph)->cor_coph
  
  append(cors_coph,paste0("a23_",m,"=",cor_coph))->cors_coph
  cors_cophfile<-file("./pvclust/cors_coph.txt")
  writeLines(cors_coph, cors_cophfile)
  close(cors_cophfile)
  }

#Same with axes 1 to 3
coord_inds<- pca18$li[, 1:3]
#cors_coph<-c()
as.matrix(dist.quant (coord_inds, method=1,diag=T,upper=T))->dist_coord_ints_13
methods=c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median","centroid")
#dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
for (m in methods) {
  pvclust(dist_coord_ints_13,method.hclust=m,method.dist="euclidean",nboot=1000)->pvclust_a13
  dev.off()
  width = 60
  height = 30
  dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
  plot(pvclust_a13)
  pvrect(pvclust_a13,alpha=0.95)
  dev.print(pdf, paste0('./pvclust/a13_', m,'.pdf'))
  dev.print(svg, paste0('./pvclust/a13_', m,'.svg'))
  dev.print(png, paste0('./pvclust/a13_', m,'.png'),width=width+1, height = height+1, unit = "cm",res=300) 
  
  #based on https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre
  cophenetic(pvclust_a13$hclust)->coph
  cor(as.dist(dist_coord_ints_13),coph)->cor_coph
  
  append(cors_coph,paste0("a13_",m,"=",cor_coph))->cors_coph
  cors_cophfile<-file("./pvclust/cors_coph.txt")
  writeLines(cors_coph, cors_cophfile)
  close(cors_cophfile)
}


pvclust(dist_coord_ints_23,method.hclust=m,method.dist="euclidean",nboot=1000)->pvclust_a23
dev.off()
width = 60
height = 30
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
plot(pvclust_a23)
pvrect(pvclust_a23,alpha=0.95)
dev.print(pdf, paste0('./pvclust/a23_mcquitty.pdf'))
dev.print(svg, paste0('./pvclust/a23_mcquitty.svg'))
dev.print(png, paste0('./pvclust/a23_mcquitty.png'),width=width+1, height = height+1, unit = "cm",res=300) 

#based on https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre
cophenetic(pvclust_a23$hclust)->coph
cor(as.dist(dist_coord_ints_23),coph)->cor_coph

append(cors_coph,paste0("a23_mcquitty=",cor_coph))->cors_coph
cors_cophfile<-file("./pvclust/cors_coph.txt")
writeLines(cors_coph, cors_cophfile)
close(cors_cophfile)

#Same with axes 1 to 3
coord_inds<- pca18$li[, 1:3]
as.matrix(dist.quant (coord_inds, method=1,diag=T,upper=T))->dist_coord_ints_13
pvclust(dist_coord_ints_13,method.hclust="mcquitty",method.dist="euclidean",nboot=1000)->pvclust_a13

dev.off()
width = 60
height = 30
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
plot(pvclust_a13)
pvrect(pvclust_a13,alpha=0.95)

#test correlation distance

setwd("E:/BARCODE/190924 IN PROGRESS")
load("E:/BARCODE/190924 IN PROGRESS/.RData")

coord_inds<- pca18$li[, 2:3]
cors_coph<-c()
as.matrix(dist.quant (coord_inds, method=1,diag=T,upper=T))->dist_coord_ints_23
methods=c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median","centroid")
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
for (m in methods) {
  pvclust(dist_coord_ints_23,method.hclust=m,method.dist="correlation",nboot=1000)->pvclust_a23
  dev.off()
  width = 60
  height = 30
  dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
  plot(pvclust_a23)
  pvrect(pvclust_a23,alpha=0.95)
  dev.print(pdf, paste0('./pvclust_corrdist/a23_', m,'.pdf'))
  dev.print(svg, paste0('./pvclust_corrdist/a23_', m,'.svg'))
  dev.print(png, paste0('./pvclust_corrdist/a23_', m,'.png'),width=width+1, height = height+1, unit = "cm",res=300) 
  
  #based on https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre
  cophenetic(pvclust_a23$hclust)->coph
  cor(as.dist(dist_coord_ints_23),coph)->cor_coph
  
  append(cors_coph,paste0("a23_",m,"=",cor_coph))->cors_coph
  cors_cophfile<-file("./pvclust_corrdist/cors_coph.txt")
  writeLines(cors_coph, cors_cophfile)
  close(cors_cophfile)
}

#Same with axes 1 to 3
coord_inds<- pca18$li[, 1:3]
#cors_coph<-c()
as.matrix(dist.quant (coord_inds, method=1,diag=T,upper=T))->dist_coord_ints_13
methods=c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median","centroid")
#dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
for (m in methods) {
  pvclust(dist_coord_ints_13,method.hclust=m,method.dist="correlation",nboot=1000)->pvclust_a13
  dev.off()
  width = 60
  height = 30
  dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
  plot(pvclust_a13)
  pvrect(pvclust_a13,alpha=0.95)
  dev.print(pdf, paste0('./pvclust_corrdist/a13_', m,'.pdf'))
  dev.print(svg, paste0('./pvclust_corrdist/a13_', m,'.svg'))
  dev.print(png, paste0('./pvclust_corrdist/a13_', m,'.png'),width=width+1, height = height+1, unit = "cm",res=300) 
  
  #based on https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre
  cophenetic(pvclust_a13$hclust)->coph
  cor(as.dist(dist_coord_ints_13),coph)->cor_coph
  
  append(cors_coph,paste0("a13_",m,"=",cor_coph))->cors_coph
  cors_cophfile<-file("./pvclust_corrdist/cors_coph.txt")
  writeLines(cors_coph, cors_cophfile)
  close(cors_cophfile)
}


#### 4.6 multispecies aveage by cluster ####

library(sp) #required for raster
library(raster) #used for most shapefiles and raster manipulation / projection
library(seqinr) #used for col2alpha
library(phylin) #used for IDW, midpoints, extract.val
library(mapplots) #used for draw.shape rivers
library(ade4) #for mantel.rtest and PCA
library(xlsx) #write.xlsx
library(dplyr) #%>%
library(MCRestimate) #to replace NAs in PCA df
library(png)
library(rasterVis)
library(pdftools)
library(lattice)
library(pvclust)

dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
#plot.new()
t=1 #t=nb species
#xs=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
#ys=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
#nx=5
#ny=5
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10

#lz est la liste qui contient les 179 df avec x, y et z

#Need to add All_species in lz
lz->lzm

as.data.frame(n_means)->lzm[[180]]
names(lzm)<-c(names(lz),"All_species")

rm(clusters)
clusters<-list()
#in notepad++, enter this: [\r\n]+ and replace by lzm$
clusters[[1]] <- rbind(lzm$Guyanancistrus_brevispinis,lzm$Pimelodella_leptosoma)
clusters[[2]] <- rbind(lzm$Auchenipterus_nuchalis, lzm$Serrasalmus_rhombeus, lzm$Hypostomus_gymnorhynchus, lzm$Pimelodus_ornatus)
clusters[[3]] <- rbind(lzm$Gasteropelecus_sternicla, lzm$Helogenes_marmoratus, lzm$Eigenmannia_virescens, lzm$Leporinus_fasciatus)
clusters[[4]] <- rbind(lzm$Leporinus_maculatus,lzm$Brycon_pesu,lzm$Acestrorhynchus_microlepis,lzm$Leporinus_granti)
clusters[[5]] <- rbind(lzm$Bryconops_affinis,lzm$Cyphocharax_helleri)
clusters[[6]] <- rbind(lzm$Characidium_zebra,lzm$Roeboexodon_guyanensis,lzm$Bryconops_caudomaculatus,lzm$Harttia_guianensis)
clusters[[7]] <- rbind(lzm$Prochilodus_rubrotaeniatus,lzm$Chalceus_macrolepidotus,lzm$Curculionichthys_sp._Maroni)
clusters[[8]] <- rbind(lzm$Cteniloricaria_platystoma,lzm$Parodon_guyanensis)
clusters[[9]] <- rbind(lzm$Ageneiosus_inermis,lzm$Triportheus_brachipomus)
clusters[[10]] <- rbind(lzm$Caenotropus_maculosus,lzm$Jupiaba_keithi,lzm$Hemiodus_huraulti,lzm$Cynopotamus_essequibensis,lzm$Semaprochilodus_varii)
clusters[[11]] <- rbind(lzm$Anostomus_brevior,lzm$Pristobrycon_striolatus)
clusters[[12]] <- rbind(lzm$Pimelabditus_moli,lzm$Bryconamericus_guyanensis,lzm$Copella_carsevennensis,lzm$Moenkhausia_oligolepis)
clusters[[13]] <- rbind(lzm$Myloplus_rhomboidalis,lzm$Erythrinus_erythrinus,lzm$Melanocharacidium_sp._2,lzm$Hypopomus_artedi,lzm$Moenkhausia_aff._colletti,lzm$Hypopygus_lepturus,lzm$Gymnotus_coropinae,lzm$Jupiaba_abramoides)
clusters[[14]] <- rbind(lzm$Pseudancistrus_barbatus,lzm$Crenicichla_multispinosa,lzm$Nannostomus_bifasciatus,lzm$Pimelodella_cf._cristata,lzm$Rhamphichthys_rostratus)
clusters[[15]] <- rbind(lzm$Moenkhausia_intermedia,lzm$Tetragonopterus_georgiae,lzm$Cleithracara_maronii,lzm$Aequidens_tetramerus,lzm$Myloplus_ternetzi,lzm$Hoplias_malabaricus,lzm$Jupiaba_meunieri,lzm$Platydoras_costatus)
clusters[[16]] <- rbind(lzm$Doras_micropoeus,lzm$Cynodon_meionactis,lzm$Hoplias_aimara)
clusters[[17]] <- rbind(lzm$Cichla_ocellaris,lzm$Jupiaba_maroniensis,lzm$Loricaria_aff._nickeriensis,lzm$Metaloricaria_paucidens)
clusters[[18]] <- rbind(lzm$Hemiodus_unimaculatus,lzm$Knodus_heteresthes)
clusters[[19]] <- rbind(lzm$Pimelodella_geryi,lzm$Loricaria_cataphracta,lzm$Chasmocranus_brevior,lzm$Harttiella_crassicauda,lzm$Corydoras_aff._breei,lzm$Harttiella_lucifer,lzm$All_species)
clusters[[20]] <- rbind(lzm$Hemisorubim_platyrhynchos,lzm$Cyphocharax_biocellatus,lzm$Heptapterus_tapanahoniensis,lzm$Peckoltia_otali,lzm$Corydoras_aff._guianensis,lzm$Aphyocharacidium_melandetum,lzm$Serrapinnus_gracilis)
clusters[[21]] <- rbind(lzm$Corydoras_geoffroy,lzm$Potamotrygon_marinae,lzm$Thayeria_ifati)
clusters[[22]] <- rbind(lzm$Imparfinis_sp.,lzm$Melanocharacidium_dispilomma)
clusters[[23]] <- rbind(lzm$Ituglanis_amazonicus,lzm$Myloplus_planquettei)
clusters[[24]] <- rbind(lzm$Moenkhausia_moisae,lzm$Myloplus_rubripinnis)
clusters[[25]] <- rbind(lzm$Phenacogaster_wayana,lzm$Bryconops_melanurus,lzm$Hypostomus_plecostomus)
clusters[[26]] <- rbind(lzm$Brycon_falcatus,lzm$Schizodon_fasciatus,lzm$Serrasalmus_eigenmanni)
clusters[[27]] <- rbind(lzm$Geophagus_surinamensis,lzm$Plagioscion_squamosissimus,lzm$Pseudacanthicus_serratus,lzm$Doras_carinatus,lzm$Pachypops_aff._fourcroi,lzm$Chasmocranus_longior,lzm$Tometes_lebaili)
clusters[[28]] <- rbind(lzm$Cephalosilurus_nigricaudus,lzm$Hemiancistrus_medians)


for (cluster in clusters) {

cluster %>%                     #applies the following function to dfn
  group_by(x, y) %>%
  summarise(meanZ = mean(Z))->cluster_means #formerly just n

colnames(cluster_means)<-c("x","y","Z")


rasterFromXYZ(cluster_means)->cluster_xyz_total
crop <- setValues(cluster_xyz_total, NA)
myshp.r <- rasterize(a, crop)
cluster_xyz_total.masked <- mask(x=cluster_xyz_total, mask=myshp.r)

#with levelplot
#prepare shared rasters and layers for plots (can put at start of the script)
#crop(basemap,ext_basemap)->basemap2  
myTheme <- BTCTheme()
myTheme$panel.background$col = 'dodgerblue3' #for the sea
levelplot(basemap2,col.regions=(pal(1000)),margin=F,colorkey=FALSE)->rast1 #rast1= croppped south america
layer(sp.polygons(maroni_basins))->lay1 #basins limits and rivers
extend(cluster_xyz_total.masked,ext_basemap)->cluster_xyz_total.masked2

#plot.new()
levelplot(cluster_xyz_total.masked2,par.settings = myTheme,col.regions=cols_tr,at=seq(0, 1, length.out=11), colorkey = F,margin=F,main=list(label=paste0("Cluster ",t),cex=1),scales=list(draw=F))->rast2  #colorkey(title=label_norm,title.gpar=list(fontsize=20))
#rast2+layer(sp.points(SpatialPoints(coord_tot[,1:2]),pch=21,col='black',fill="white",cex=1))->r2l2
rast2+as.layer(rast1+lay1,under=T)->plot
print(plot,split=c(xs[t], ys[t], nx, ny),newpage=F)
t=t+1
}

dev.print(pdf, paste0('./pvclust/pvclust_maps.pdf'))
dev.print(png, paste0('./pvclust/pvclust_maps.png'),width = 1024, height = 768)
dev.off()

