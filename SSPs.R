library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(gplots)
library(tidyr)
library(Hmisc)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(pivottabler)
library(VennDiagram)

setwd("~/Documents/colleen_brachy/Data/SSPs/")
all_genotype_sig<-read.csv("~/Documents/colleen_brachy/Data/all_sig_by_genotype.csv")
all_timepoint_sig<-read.csv("~/Documents/colleen_brachy/Data/all_sig_by_timepoint.csv")

files <- dir(pattern =".fa.fai")
species<-substr(files, 1, nchar(files)-7)
all<-list()
list<-list()
for(file in species){
  length<-read.delim(paste(file, ".fa.fai", sep=""), header=F)
  colnames(length)<-c("Protein", "Length", "Offset", "Linebases", "LineWidth")
  signalp<-read.delim(paste(file, "_summary.signalp5", sep=""),  header=FALSE, comment.char="#")
  colnames(signalp)<-c("Protein", "Prediction", "SP", "Other", "CS.position")
  join<-merge(length, signalp)
  small<-sum(join$Length<=250)
  signal<-sum(!(join$Prediction=="OTHER"))
  both<-dim(join[(join$Length <= 250 & !(join$Prediction == "OTHER")),])[1]
  table<-cbind(file, small, signal, both)
  list[[length(list)+1]] = table
  join$Chromosome<-paste(file)
  write.csv(join, file=paste(file, "all.csv", sep=""))
  ssps<-join[(join$Length<=250),]
  ssps<-ssps[!(ssps$Prediction=="OTHER"),]
  all[[length(all)+1]]=ssps
  write.csv(ssps, file=paste(file, "ssps.csv", sep=""))
}
Together<-as.data.frame(do.call("rbind", list))
all_putative_SSPs<-as.data.frame(do.call("rbind", all))
write.csv(Together, file="table_ssps.csv")
write.csv(all_SSPs, file="all_chromosomes_ssps.csv")

zymo_tmm <- read.csv("~/Documents/colleen_brachy/Data/SSPs/zymo_tmm.csv")
Brachy_tmm <- read.csv("~/Documents/colleen_brachy/Data/SSPs/Brachy_tmm.csv")
zymo<-merge(all_SSPs[(all_SSPs$Chromosome=="Zymoseptoria_tritici.MG2.pep.all" ),], zymo_tmm, by="Protein")
brachy<-merge(all_SSPs[(all_SSPs$Chromosome=="Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all" ),], Brachy_tmm, by="Protein")

brachy_ssps<-brachy[(brachy$PredHel == 0),]
brachy_ssps<-brachy_ssps[(brachy_ssps$GPI_anchor=="No"),]

zymo_ssps<-zymo[(zymo$PredHel == 0),]
zymo_ssps<-zymo_ssps[(zymo_ssps$GPI_anchor=="No"),]

all_ssps<-rbind(brachy_ssps, zymo_ssps)
all_ssps$Protein<-gsub("Mycgr3P", "Mycgr3T", all_ssps$Protein)

ssps_deg_isolate<-subset(all_genotype_sig, all_genotype_sig$row %in% all_ssps$Protein)
ssps_deg_time<-subset(all_timepoint_sig, all_timepoint_sig$row %in% all_ssps$Protein)
write.csv(ssps_deg_isolate, file="~/Documents/colleen_brachy/Data/SSPs/DESSPs_isolate.csv")
write.csv(ssps_deg_time, file="~/Documents/colleen_brachy/Data/SSPs/DESSPs_time.csv")


