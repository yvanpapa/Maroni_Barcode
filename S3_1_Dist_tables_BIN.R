####1: CREATE DIST TABLES FROM DNA AND LOCS####
#This script computes distance matrices and other useful tables for
#each species.

#Inputs needed: 
#1)Maroni_loc.txt   Tab-delimited text file with six columns: Species, Lat, Lon, field.num, Exact_Site, BIN
#2)all_seq.fas      Fasta file containing all aligned sequences. Name of sequences must be coded as Genus_species|Sample_code

#Useful tip if all_seq.fas comes fom BOLD: use GBOL.*?[|] to remove the GBOL characters
#Check carefully that species names match between Maroni loc and all_seq

####MUST BE DEFINED BY USER #####
session1<-"your_session_directory_here"
setwd(session1)

####ALL PACKAGES##### not including their dependencies
library(ape) 
#"as.DNAbin","read.dna","dist.dna"

dir.create("outputs")
dir.create("outputs/matdist")
dir.create("outputs/att_tables")

#dfloc: df with 6 c: Species, Lat, Lon, field num, exact site, BIN, ordered by field.num
dfloc<-read.table("Maroni_loc.txt",header = T)
dfloc<- dfloc[order(dfloc$field.num),] 

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

#dflocdna: merge df_split and dfloc to get df with 6 c: sp, seq, field.num, Lat, Lon, exact site, BIN
dflocdna<-data.frame(cbind(df_split[,1],dfdna$seq,df_split[,2]),stringsAsFactors = FALSE)
dflocdna<- dflocdna[order(dflocdna$X3),]
dflocdna<-cbind(dflocdna,dfloc$Lat,dfloc$Lon,dfloc$Exact_Site,dfloc$BIN)
colnames(dflocdna)<-c("sp","seq","field.num","Lat","Lon","Exact_Site","BIN")

#list_dflocdna_bin: list of n=number of specimens elements. Each element is a data frame with seq, coord, field num for each bin.
#Because of lapply, outputs as many identical df for each bin as there are specimens in each bin
lapply(dflocdna$BIN,function(i) assign(paste("df", i, sep = "."), dflocdna[dflocdna$BIN==i,]))->list_dflocdna_BIN

#list_unique_dflocdna_BIN:list of n=number of species elements. We just removed the identical dfs
unique(list_dflocdna_BIN)->list_unique_dflocdna_BIN

#Change equal latitudes with +0.00001 for each BIN
for (j in 1:length(list_unique_dflocdna_BIN)) {
  k=0.00001
  for (i in 1:length(list_unique_dflocdna_BIN[[j]]$Lat)) {
    if (list_unique_dflocdna_BIN[[j]]$Lat[i] %in% list_unique_dflocdna_BIN[[j]]$Lat[-i]) {
      list_unique_dflocdna_BIN[[j]]$Lat[i]<-list_unique_dflocdna_BIN[[j]]$Lat[i]+k
      k=k+0.00001
    }}}
rm(i,j,k)

#list_seqmats: list of n= nb of BIN matrices. Each matrix contains dna sequences with nb of rows= nb of specimens/BIN
lapply(list_unique_dflocdna_BIN,function(i) t(sapply(strsplit(i[,2],""), tolower)))->list_seqmats
#list_BIN_field.nums (formerly "a"): list of n= nb of BIN vectors. Each vector contains the field.num of each specimen
lapply(list_unique_dflocdna_BIN,function(i) print(i[,3]))->list_bin_field.nums

#list_seqmats2: attribute field.num to each sequence
list()->list_seqmats2
for (i in 1:length(list_bin_field.nums)) {
  list_bin_field.nums[[i]]->names
  list_seqmats[[i]]->m
  names->rownames(m)
  list_seqmats2[[i]]<-m
}
rm(i,m,names)

#apply dist.dna in each group of sequences in list_seqmats 2 -> list_mat_dist
lapply(list_seqmats2,function(i) dist.dna(as.DNAbin(i)))->list_mat_dist

unique(dflocdna$BIN)->v_BINs
gsub(":", "_", v_BINs)->v_BINs #need to replace : by _ or you can't save the file

#save genetic distance as matrices:
for (i in 1:length(list_mat_dist)) {
  write.table(as.matrix(list_mat_dist[[i]]), paste0("outputs/matdist/matdist.",v_BINs[i],".txt"), sep="\t",row.names = T,quote=F)
}
rm(i)

#save attribute tables: (dont worry too much about 'numcode' (formerly Popcode. Just a quick way to index unique pairs of coords for later))
list()->list_att_tables
for(i in 1:length(list_unique_dflocdna_BIN)) {
  list_unique_dflocdna_BIN[[i]]->att
  att[,-c(2,7)]->att
  rownames(att)->att[,6]
  c("sp","Sample_ID","Latitude","Longitude","Exact_Site","numcode")->colnames(att)
  att<-att[,c(6,2:5,1)]
  list_att_tables[[i]]<-att
}
rm(i,att)
names(list_att_tables)<-v_BINs

for (i in 1:length(list_att_tables)) {
  write.table(list_att_tables[[i]], paste0("outputs/att_tables/att.",v_BINs[i],".txt"), sep="\t",row.names = F,quote=F)
}
rm(i)

#no more needed
rm(df_split,dfdna,dfloc,dflocdna,list_dflocdna_BIN,list_seqmats,list_seqmats2,list_bin_field.nums,list_split,list_unique_dflocdna_BIN,myDNA,session1)

#needed for next steps:
#rm(list_att_tables,list_mat_dist,v_BINs)

save.image("env_1_Dist_tables.RData")
