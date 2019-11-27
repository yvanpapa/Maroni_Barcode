#PLEASE REFER TO README.txt TO KNOW REQUIRED INPUTS

session<-"your_session_directory_here"
session<-"E:/BARCODE_MARONI/SCRIPTS_AND_DATA/FIGURE2"
setwd(session)

library(reshape2)
#melt
library(ggplot2)
#ggplot



#How many sp by genus?
taxonomy<-read.table("Taxonomy_spreadsheet.txt", header=T, row.names=NULL,sep="\t")
gensp<-as.data.frame(taxonomy[,8:9])  #only keep Genus and Species column
as.data.frame(table(unique(gensp)$Genus))->genspcount

length(which(genspcount$Freq>1))
length(which(genspcount$Freq>1))/124

#How many individuals by species
as.data.frame(table(taxonomy$Species))->spid
min(spid$Freq)
max(spid$Freq)
mean(spid$Freq)
length(which(spid$Freq>2))
length(which(spid$Freq>2))/length(spid$Freq)
length(which(spid$Freq==1))

#How many sp by order
tax<-as.data.frame(taxonomy[,c(4,5,8,9)])
as.data.frame(table(unique(tax[1:2])$Order))->FamilyFreq
as.data.frame(table(unique(tax[1:3])$Order))->GenFreq
as.data.frame(table(unique(tax[1:4])$Order))->SpFreq
Ordercounts<-as.data.frame(FamilyFreq)
Ordercounts$Genus<-GenFreq$Freq
Ordercounts$Species<-SpFreq$Freq
colnames(Ordercounts)<-c("Order","Family","Genus","Species")
data.m <- melt(Ordercounts, id.vars='Order')

ggplot(data.m, aes(Order, value)) +   
geom_bar(aes(fill = variable),position = "dodge", stat="identity")+
scale_fill_grey()+
guides(fill=guide_legend(title=NULL))+
theme_bw()+
theme(text = element_text(size=18),axis.text.x = element_text(angle=45, hjust=1,vjust=1),axis.title.y=element_blank())

save.image(paste0(session,"/Figure2_data.RData"))
