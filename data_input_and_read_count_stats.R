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
colData$Factor<-paste(colData$Zymo_isolate, colData$Timepoint)
# ==================================================================================
# read count stats chunk starts here

expressed_genes<-txi.kallisto.tsv$counts
expressed_genes<-as.data.frame(expressed_genes)
expressed_genes$GeneID<-row.names(expressed_genes)
expressed_genes<-expressed_genes[,c(37, 1:36)]
expressed_genes_long<-expressed_genes %>% gather(Sample, TPM, 2:37)
all_genes<-merge(expressed_genes_long, colData, by="Sample")
sub<-all_genes[,c(9, 2, 3, 5)]
rep_wise<-spread(sub, key = Rep, value=TPM)
rep_wise$Sum<-rep_wise$`1` + rep_wise$`2` + rep_wise$`3`
rep_wise$test1<-ifelse(rep_wise$`1`>0.5, 1,0)
rep_wise$test2<-ifelse(rep_wise$`2`>0.5, 1,0)
rep_wise$test3<-ifelse(rep_wise$`3`>0.5, 1,0)
write.csv(rep_wise, file='~/Documents/colleen_brachy/Data/repwise.csv', row.names=F)

rep_wise$sum2<-rep_wise$test1 + rep_wise$test2 + rep_wise$test3
expressed<-rep_wise[(rep_wise$sum2 >=2),] 

expressed_brachy<-subset(expressed, expressed$GeneID %in% Brachy_genes$Transcript)
expressed_brachy$Species<-'B. distachyon'
expressed_Zymo<-subset(expressed, expressed$GeneID %in% Zymo_genes$Transcript)
expressed_Zymo$Species<-("Z. tritici")
 kable(table(expressed$Factor), 
      caption="Number of expressed genes per sample (averaged across reps)",
      "pandoc", align="c", col.names=c("Sample", "Number of expressed genes"))

# check correlation of reps
cor<-as.matrix(rep_wise[,c(3,4,5)])
cor<-rcorr(cor)
corrplot(cor$r, type="lower", order="original",p.mat = cor$P, 
         sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 2, 
         tl.srt = 0, tl.offset = 1, method="color", addCoef.col = "white")

rep_wise$sum2<-rep_wise$test1 + rep_wise$test2 + rep_wise$test3
expressed<-rep_wise[(rep_wise$sum2 >=2),] 
length(unique(expressed$GeneID))
unique_gene<-as.data.frame(unique(expressed$GeneID))
colnames(unique_gene)<-"GeneID"


# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, 
                                colData, ~ Timepoint+ Rep + Zymo_isolate)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Timepoint", "Zymo_isolate", "Rep"), returnData=TRUE)
pcaData$Sample<-paste(pcaData$Timepoint, pcaData$Zymo_isolate)
# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(color=Rep), size=7, alpha=0.8) +
  theme_classic() +
  # scale_color_manual(values=c("slategray3", "dodgerblue", "rosybrown3", "orangered"), 
  #                    labels=c("Longbow + Tween20", expression(paste("Longbow + ",italic("Z. tritici"))), "Stigg + Tween20", expression(paste("Stigg + ",italic("Z. tritici")))))+
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black")) +
  xlab("Principle component 1")+
  ylab("Principle component 2")

Zymo_genes <- read.delim("~/Documents/colleen_brachy/Data/Zymo_genes.txt")
Brachy_genes <- read.delim("~/Documents/colleen_brachy/Data/Brachy_genes.txt")

zymo_dds<-dds[as.character(Zymo_genes$Transcript),]
brachy_dds<-dds[as.character(Brachy_genes$Transcript),]
zymo_vsd <- varianceStabilizingTransformation(zymo_dds, blind=FALSE)
brachy_vsd <- varianceStabilizingTransformation(brachy_dds, blind=FALSE)


#  Zymo pca
zymo_pcaData <- plotPCA(zymo_vsd, intgroup=c("Timepoint", "Zymo_isolate", "Rep"), returnData=TRUE)
zymo_pcaData$Sample<-paste(zymo_pcaData$Timepoint, zymo_pcaData$Zymo_isolate)
write.csv(zymo_pcaData, file="~/Documents/colleen_brachy/Data/zymo_PCA.csv", row.names=F)
# plot PC1 vs PC2
ggplot(zymo_pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(color=Rep, shape=Timepoint), size=7, alpha=0.8) +
  theme_classic() +
  # scale_color_manual(values=c("slategray3", "dodgerblue", "rosybrown3", "orangered"), 
  #                    labels=c("Longbow + Tween20", expression(paste("Longbow + ",italic("Z. tritici"))), "Stigg + Tween20", expression(paste("Stigg + ",italic("Z. tritici")))))+
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black")) +
  xlab("Principle component 1")+
  ylab("Principle component 2")

#  Zymo pca
brachy_pcaData <- plotPCA(brachy_vsd, intgroup=c("Timepoint", "Zymo_isolate", "Rep"), returnData=TRUE)
brachy_pcaData$Sample<-paste(brachy_pcaData$Timepoint, brachy_pcaData$Zymo_isolate)
write.csv(brachy_pcaData, file="~/Documents/colleen_brachy/Data/brachy_PCA.csv", row.names=F)

# plot PC1 vs PC2
ggplot(brachy_pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(color=Rep), size=7, alpha=0.8) +
  theme_classic() +
  # scale_color_manual(values=c("slategray3", "dodgerblue", "rosybrown3", "orangered"), 
  #                    labels=c("Longbow + Tween20", expression(paste("Longbow + ",italic("Z. tritici"))), "Stigg + Tween20", expression(paste("Stigg + ",italic("Z. tritici")))))+
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black")) +
  xlab("Principle component 1")+
  ylab("Principle component 2")
