####1: CREATE DIST TABLES FROM DNA AND LOCS####
#This script computes distance matrices and other useful tables for
#each species.

#Inputs needed: 
#1)Maroni_loc.txt   Tab-delimited text file with five columns: Species, Lat, Lon, field.num, Exact_Site
#2)all_seq.fas      Fasta file containing all aligned sequences. Name of sequences must be coded as Genus_species|Sample_code

#Useful tip if all_seq.fas comes fom BOLD: use GBOL.*?[|] to remove the GBOL characters
#Check carefully that species names match between Maroni loc and all_seq

####MUST BE DEFINED BY USER #####
session1<-"G:/DATA/WORK/P4 MARONI BARCODE/BARCODE_MARONI/201209_back_Raph1/R_ANALYSES_SUPP/Maroni_Barcode-master/simulations"
setwd(session1)

####ALL PACKAGES##### not including their dependencies
library(ape) 
#"as.DNAbin","read.dna","dist.dna"

dir.create("outputs")
dir.create("outputs/matdist")
dir.create("outputs/att_tables")

#dfloc: df with 5 c: Species, Lat, Lon, field num, exact site, ordered by field.num
dfloc<-read.table("Maroni_loc.txt",header = T)
dfloc<- dfloc[order(dfloc$field.num),] 
dfloc_subset<-dfloc[c(which(dfloc$Species=="Cteniloricaria_platystoma"),which(dfloc$Species=="Metaloricaria_paucidens"),
                which(dfloc$Species=="Myloplus_ternetzi")),]            
dfloc_subset<- dfloc_subset[order(dfloc_subset$field.num),] 

#dfdna: input=dna sequences; output=df with 2 c: 1)Species|field.numb 2)seq as character strings
#Do not forget to change spaces into underscores before
#AND remove BOLD identifiers with GBOL.*?[|]
read.dna("all_seq.fas",format="fasta",as.character=T)->myDNA #heavy
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

#dfloc_subsetdna: merge df_split and dfloc_subset to get df with 5 c: sp, seq, field.num, Lat, Lon, exact site
dfloc_subsetdna<-data.frame(cbind(df_split[,1],dfdna$seq,df_split[,2]),stringsAsFactors = FALSE)
dfloc_subsetdna<- dfloc_subsetdna[order(dfloc_subsetdna$X3),]

dfloc_subsetdna<-dfloc_subsetdna[c(which(dfloc_subsetdna$X1=="Cteniloricaria_platystoma"),which(dfloc_subsetdna$X1=="Metaloricaria_paucidens"),
                                   which(dfloc_subsetdna$X1=="Myloplus_ternetzi")),] 

dfloc_subsetdna<- dfloc_subsetdna[order(dfloc_subsetdna$X3),] 


dfloc_subsetdna<-cbind(dfloc_subsetdna,dfloc_subset$Lat,dfloc_subset$Lon,dfloc_subset$Exact_Site)
colnames(dfloc_subsetdna)<-c("sp","seq","field.num","Lat","Lon","Exact_Site")

######
#add simulations with less localities:
#Cteniloricaria platystoma: Aweim, Esperance
#Metaloricaria paucidens, remove: Antecume, Taluen, Cayode, Papa??chton
#Myloplus ternetzi, remove: Papa??chton, Taluen, Cayode, Pidima, Cayod


center_locs<-c(grep("Antecume",dfloc_subsetdna$Exact_Site),
grep("Aweim",dfloc_subsetdna$Exact_Site),
grep("Cayod",dfloc_subsetdna$Exact_Site),
grep("Esperance",dfloc_subsetdna$Exact_Site),
grep("Papa",dfloc_subsetdna$Exact_Site),
grep("Pidima",dfloc_subsetdna$Exact_Site),
grep("Taluen",dfloc_subsetdna$Exact_Site))

dfloc_subsetdna.border<-dfloc_subsetdna[-center_locs,]

#I want to add these points for this species bc they areat the border of the distribution
rbind(dfloc_subsetdna.border,dfloc_subsetdna[which(dfloc_subsetdna$sp=="Cteniloricaria_platystoma"&
                        dfloc_subsetdna$Exact_Site=="Vicinity_of_Antecume_Pata"),])->dfloc_subsetdna.border

paste0(dfloc_subsetdna.border$sp,"_bordersim")->dfloc_subsetdna.border$sp

rbind(dfloc_subsetdna,dfloc_subsetdna.border)->dfloc_subsetdna
dfloc_subsetdna<- dfloc_subsetdna[order(dfloc_subsetdna$field.num),] 


#list_dfloc_subsetdna_sp: list of n=number of specimens elements. Each element is a data frame with seq, coord, field num for each species.
#Because of lapply, outputs as many identical df for each species as there are specimens in each species
lapply(dfloc_subsetdna$sp,function(i) assign(paste("df", i, sep = "."), dfloc_subsetdna[dfloc_subsetdna$sp==i,]))->list_dfloc_subsetdna_sp

#list_unique_dfloc_subsetdna_sp:list of n=number of species elements. We just removed the identical dfs
unique(list_dfloc_subsetdna_sp)->list_unique_dfloc_subsetdna_sp

#Change equal latitudes with +0.00001 for each sp
for (j in 1:length(list_unique_dfloc_subsetdna_sp)) {
  k=0.00001
  for (i in 1:length(list_unique_dfloc_subsetdna_sp[[j]]$Lat)) {
    if (list_unique_dfloc_subsetdna_sp[[j]]$Lat[i] %in% list_unique_dfloc_subsetdna_sp[[j]]$Lat[-i]) {
      list_unique_dfloc_subsetdna_sp[[j]]$Lat[i]<-list_unique_dfloc_subsetdna_sp[[j]]$Lat[i]+k
      k=k+0.00001
    }}}
rm(i,j,k)

#list_seqmats: list of n= nb of sp matrices. Each matrix contains dna sequences with nb of rows= nb of specimens/sp
lapply(list_unique_dfloc_subsetdna_sp,function(i) t(sapply(strsplit(i[,2],""), tolower)))->list_seqmats
#list_sp_field.nums (formerly "a"): list of n= nb of sp vectors. Each vector contains the field.num of each specimen
lapply(list_unique_dfloc_subsetdna_sp,function(i) print(i[,3]))->list_sp_field.nums

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

unique(dfloc_subsetdna$sp)->v_species

#save genetic distance as matrices:
for (i in 1:length(list_mat_dist)) {
  write.table(as.matrix(list_mat_dist[[i]]), paste0("outputs/matdist/matdist.",v_species[i],".txt"), sep="\t",row.names = T,quote=F)
}
rm(i)

#save attribute tables: (dont worry too much about 'numcode' (formerly Popcode. Just a quick way to index unique pairs of coords for later))
list()->list_att_tables
for(i in 1:length(list_unique_dfloc_subsetdna_sp)) {
  list_unique_dfloc_subsetdna_sp[[i]]->att
  att[,-c(1,2)]->att
  rownames(att)->att[,5]
  c("Sample_ID","Latitude","Longitude","Exact_Site","numcode")->colnames(att)
  att<-att[,c(5,1:4)]
  list_att_tables[[i]]<-att
}
rm(i,att)
names(list_att_tables)<-v_species

for (i in 1:length(list_att_tables)) {
  write.table(list_att_tables[[i]], paste0("outputs/att_tables/att.",v_species[i],".txt"), sep="\t",row.names = F,quote=F)
}
rm(i)

#no more needed
rm(df_split,dfdna,dfloc_subset,dfloc_subsetdna,list_dfloc_subsetdna_sp,list_seqmats,list_seqmats2,list_sp_field.nums,list_split,list_unique_dfloc_subsetdna_sp,myDNA,session1)

#needed for next steps:
#rm(list_att_tables,list_mat_dist,v_species)

save.image("env_1_Dist_tables.RData")
