######### 4: WPGMA ######### (or any of the best clustering method) ####
### !!! If you don't want to run the test methods, skkp to "run pvclust by using WPGMA"

#Inputs needed:
#pca18 (in 4_PCA environment) 

session4<-"your_session4_directory_here"
load(paste0(session4,"/env_4_PCA_all_data.RData")) 

session5<-"your_session_directory_here"
setwd(session5)

dir.create("clust_tests")
dir.create("pvclust")

library(ade4)
library(pvclust)

####TEST SEVERAL CLUSTERING METHODS RUN TIME 1-2 HOURS ####

coord_inds<- pca18$li[, 2:3]
cors_coph<-c()
as.matrix(dist.quant (coord_inds, method=1,diag=T,upper=T))->dist_coord_ints_23
methods=c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median","centroid")
dev.off()
for (m in methods) {
  pvclust(dist_coord_ints_23,method.hclust=m,method.dist="euclidean",nboot=1000)->pvclust_a23
  width = 60
  height = 30
  dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
  plot(pvclust_a23)
  pvrect(pvclust_a23,alpha=0.95)
  dev.print(pdf, paste0('./clust_tests/a23_', m,'.pdf'))
  dev.print(svg, paste0('./clust_tests/a23_', m,'.svg'))
  dev.print(png, paste0('./clust_tests/a23_', m,'.png'),width=width+1, height = height+1, unit = "cm",res=300) 
  dev.off()
  
  #based on https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre
  cophenetic(pvclust_a23$hclust)->coph
  cor(as.dist(dist_coord_ints_23),coph)->cor_coph
  
  append(cors_coph,paste0("a23_",m,"=",cor_coph))->cors_coph
  cors_cophfile<-file("./clust_tests/cors_coph.txt")
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
  width = 60
  height = 30
  dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
  plot(pvclust_a13)
  pvrect(pvclust_a13,alpha=0.95)
  dev.print(pdf, paste0('./clust_tests/a13_', m,'.pdf'))
  dev.print(svg, paste0('./clust_tests/a13_', m,'.svg'))
  dev.print(png, paste0('./clust_tests/a13_', m,'.png'),width=width+1, height = height+1, unit = "cm",res=300) 
  dev.off()
  
  #based on https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre
  cophenetic(pvclust_a13$hclust)->coph
  cor(as.dist(dist_coord_ints_13),coph)->cor_coph
  
  append(cors_coph,paste0("a13_",m,"=",cor_coph))->cors_coph
  cors_cophfile<-file("./clust_tests/cors_coph.txt")
  writeLines(cors_coph, cors_cophfile)
  close(cors_cophfile)
}

### RUN PVCLUST BY USING WPGMA (MCQUITTY) ON AXES 2 AND 3 ####
library(ade4)
library(pvclust)

coord_inds<- pca18$li[, 2:3]
cors_coph<-c()
as.matrix(dist.quant (coord_inds, method=1,diag=T,upper=T))->dist_coord_ints_23


pvclust(dist_coord_ints_23,method.hclust="mcquitty",method.dist="euclidean",nboot=1000)->pvclust_a23
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
pvclust(dist_coord_ints_13,method.hclust="mcquitty",method.dist="euclidean",nboot=10)->pvclust_a13

dev.off()
width = 60
height = 30
dev.new(noRStudioGD = T,width = width, height = height, unit = "cm")
plot(pvclust_a13)
pvrect(pvclust_a13,alpha=0.95)

dev.print(pdf, paste0('./pvclust/a13_mcquitty.pdf'))
dev.print(svg, paste0('./pvclust/a13_mcquitty.svg'))
dev.print(png, paste0('./pvclust/a13_mcquitty.png'),width=width+1, height = height+1, unit = "cm",res=300)

cophenetic(pvclust_a13$hclust)->coph
cor(as.dist(dist_coord_ints_13),coph)->cor_coph

append(cors_coph,paste0("a13_mcquitty=",cor_coph))->cors_coph
cors_cophfile<-file("./pvclust/cors_coph.txt")
writeLines(cors_coph, cors_cophfile)
close(cors_cophfile)

save.image("env_5_WPGMA.RData")
