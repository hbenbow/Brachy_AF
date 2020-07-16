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

allowWGCNAThreads()

# Import kallisto files with txi 
# ==================================================================================
# Use these steps to import kallisto files
samples<-dir("~/Documents/colleen_brachy/aligned/") # where the directory 'samples'
# contains the kallisto output directories - 1 per sample.
files <- file.path(samples, "abundance.h5")
setwd("~/Documents/colleen_brachy/aligned/")
names(files) <- paste0(samples)
all(file.exists(files)) # need to navigate into samples directory for this to work!!
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
colData<-read.csv("~/Documents/colleen_brachy/Data/colData.csv")
setwd("~/Documents/colleen_brachy/Data/")
colData$Timepoint<-as.factor(colData$Timepoint)
colData$Rep<-as.factor(colData$Rep)
# ==================================================================================


# if already have a txi object, load it with the metadata (colData)
# load("~/Documents/bmc/Data/txi_and_colData.RData")

# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))
colData$Factor<-paste(colData$Zymo_isolate, colData$Timepoint, sep="_")
colData$Factor<-as.factor(colData$Factor)
# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, 
                                colData, ~  Factor)
dds <- DESeq(dds)

# First lets compare each timepoint to T0. Here, we may have to ignore brachy reads and only look for DE zymo
# genes, as there is no way of telling if the brachy genes are DE because of treatment, or normal development
RES1<-as.data.frame(results(dds, contrast=c("Factor", "553.11_4", "553.11_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES1$Timepoint<-4
RES1$Isolate<-"553.11"
RES2<-as.data.frame(results(dds, contrast=c("Factor", "553.11_9", "553.11_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES2$Timepoint<-9
RES2$Isolate<-"553.11"
RES3<-as.data.frame(results(dds, contrast=c("Factor", "553.11_21", "553.11_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES3$Timepoint<-21
RES3$Isolate<-"553.11"
RES4<-as.data.frame(results(dds, contrast=c("Factor", "560.11_4", "560.11_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES4$Timepoint<-4
RES4$Isolate<-"560.11"
RES5<-as.data.frame(results(dds, contrast=c("Factor", "560.11_9", "560.11_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES5$Timepoint<-9
RES5$Isolate<-"560.11"
RES6<-as.data.frame(results(dds, contrast=c("Factor", "560.11_21", "560.11_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES6$Timepoint<-21
RES6$Isolate<-"560.11"
RES7<-as.data.frame(results(dds, contrast=c("Factor", "IPO323_4", "IPO323_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES7$Timepoint<-4
RES7$Isolate<-"IPO323"
RES8<-as.data.frame(results(dds, contrast=c("Factor", "IPO323_9", "IPO323_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES8$Timepoint<-9
RES8$Isolate<-"IPO323"
RES9<-as.data.frame(results(dds, contrast=c("Factor", "IPO323_21", "IPO323_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES9$Timepoint<-21
RES9$Isolate<-"IPO323"

all_by_timepoint<-rbind(RES1, RES2, RES3, RES4, RES5, RES6, RES7, RES8, RES9)
all_timepoint_sig<-all_by_timepoint[(all_by_timepoint$padj < 0.05),]
all_timepoint_sig<-na.omit(all_timepoint_sig)
rm(RES1, RES2, RES3, RES4, RES5, RES6, RES7, RES8, RES9)


# now lets compare isolates between timepoints.
RES1<-as.data.frame(results(dds, contrast=c("Factor", "553.11_0", "IPO323_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES1$Timepoint<-0
RES1$Comparison<-"553.11 v IPO323"
RES2<-as.data.frame(results(dds, contrast=c("Factor", "553.11_4", "IPO323_4"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES2$Timepoint<-4
RES2$Comparison<-"553.11 v IPO323"
RES3<-as.data.frame(results(dds, contrast=c("Factor", "553.11_9", "IPO323_9"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES3$Timepoint<-9
RES3$Comparison<-"553.11 v IPO323"
RES4<-as.data.frame(results(dds, contrast=c("Factor", "553.11_21", "IPO323_21"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES4$Timepoint<-21
RES4$Comparison<-"553.11 v IPO323"
RES5<-as.data.frame(results(dds, contrast=c("Factor", "560.11_0", "IPO323_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES5$Timepoint<-0
RES5$Comparison<-"560.11 v IPO323"
RES6<-as.data.frame(results(dds, contrast=c("Factor", "560.11_4", "IPO323_4"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES6$Timepoint<-4
RES6$Comparison<-"560.11 v IPO323"
RES7<-as.data.frame(results(dds, contrast=c("Factor", "560.11_9", "IPO323_9"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES7$Timepoint<-9
RES7$Comparison<-"560.11 v IPO323"
RES8<-as.data.frame(results(dds, contrast=c("Factor", "560.11_21", "IPO323_21"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES8$Timepoint<-21
RES8$Comparison<-"560.11 v IPO323"
RES9<-as.data.frame(results(dds, contrast=c("Factor", "560.11_0", "553.11_0"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES9$Timepoint<-0
RES9$Comparison<-"560.11 v 553.11"
RES10<-as.data.frame(results(dds, contrast=c("Factor", "560.11_4", "553.11_4"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES10$Timepoint<-4
RES10$Comparison<-"560.11 v 553.11"
RES11<-as.data.frame(results(dds, contrast=c("Factor", "560.11_9", "553.11_9"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES11$Timepoint<-9
RES11$Comparison<-"560.11 v 553.11"
RES12<-as.data.frame(results(dds, contrast=c("Factor", "560.11_21", "553.11_21"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES12$Timepoint<-21
RES12$Comparison<-"560.11 v 553.11"

all_genotype<-rbind(RES1, RES2, RES3, RES4, RES5, RES6, RES7, RES8, RES9, RES10, RES11, RES12)
all_genotype_sig<-all_genotype[(all_genotype$padj < 0.05),]
all_genotype_sig<-na.omit(all_genotype_sig)

rm(RES1, RES2, RES3, RES4, RES5, RES6, RES7, RES8, RES9, RES10, RES11, RES12)

Zymo_genes <- read.delim("~/Documents/colleen_brachy/Data/Zymo_genes.txt")
Brachy_genes <- read.delim("~/Documents/colleen_brachy/Data/Brachy_genes.txt")
Zymo_genes$Species<-'Zt'
Brachy_genes$Species<-"Bd"
genes<-rbind(Zymo_genes, Brachy_genes)
write.csv(genes, file="~/Documents/colleen_brachy/Data/all_genes.csv")
all_genotype_sig<-merge(all_genotype_sig, genes, by.x="row", by.y="Transcript")
all_timepoint_sig<-merge(all_timepoint_sig, genes, by.x="row", by.y="Transcript")
write.csv(all_genotype_sig, file="~/Documents/colleen_brachy/Data/all_sig_by_genotype.csv")
write.csv(all_timepoint_sig, file="~/Documents/colleen_brachy/Data/all_sig_by_timepoint.csv")


