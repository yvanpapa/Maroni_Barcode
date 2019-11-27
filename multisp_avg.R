#### 6 MUTISPECIES AVERAGE BY PVCLUST CLUSTER ####

#Inputs needed:
#lz, n_means (in 3_PCA_PREP)

session3<-"your_session3_directory_here"
session3<-"E:/BARCODE_MARONI/SCRIPTS_AND_DATA/LANDSCAPES_PIPELINE/3_PCA_PREP"
load(paste0(session3,"/env_3_PCA_prep_all_data.RData")) #~10 sec

Fig4<-"your_Fig4_directory_here"
Fig4<-"E:/BARCODE_MARONI/SCRIPTS_AND_DATA/FIGURE4"
load(paste0(Fig4,"/Figure4_data.RData")) #~10 sec

session6<-"your_session_directory_here"
session6<-"E:/BARCODE_MARONI/SCRIPTS_AND_DATA/LANDSCAPES_PIPELINE/6_multisp_avg"
setwd(session6)

library(dplyr)
#group_by
library(ade4)
library(rasterVis)

dev.off()
width = 25
height = 16
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
t=1 #t=nb species
xs=c(1:10,1:10,1:10,1:10,1:10)
ys=rep(1:5,each=10)
ny=5
nx=10

#lz is the ist which contains 179 df with x, y and z
#lz est la liste qui contient les 179 df avec x, y et z

#Need to add All_species in lz
lz->lzm
#We add all_species to this list
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
  myshp.r <- rasterize(maroni_shapefile, crop)
  cluster_xyz_total.masked <- mask(x=cluster_xyz_total, mask=myshp.r)
  
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

dev.print(pdf, paste0('./pvclust_maps.pdf'))
dev.print(png, paste0('./pvclust_maps.png'),width = 1024, height = 768)
dev.off()

save.image("env_6_multisp_avg.RData")

