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
setwd("~/Documents/colleen_brachy/Data/")
genes<-read.csv("~/Documents/colleen_brachy/Data/all_genes.csv")
all_genotype_sig<-read.csv("~/Documents/colleen_brachy/Data/all_sig_by_genotype.csv")
all_timepoint_sig<-read.csv("~/Documents/colleen_brachy/Data/all_sig_by_timepoint.csv")

# make some frequency tables, from which we can do venn diagrams (childish but useful)
for(t in unique(all_genotype_sig$Timepoint)){
  data<-all_genotype_sig[(all_genotype_sig$Timepoint == t),]
  tab<- spread(as.data.frame(
    table(data$row, data$Comparison)), key="Var2", value="Freq")
  assign(paste("table", t, sep="_"), tab)
  write.csv(tab, file=paste("~/Documents/colleen_brachy/Data/subset", t, ".csv", sep=""))
}

{data = table_554.11
name<-"table"
  likes <- function(animals) {
  ppl <- data
  names(ppl) <- colnames(data)
  for (i in 1:length(animals)) {
    ppl <- subset(ppl, ppl[animals[i]] == T)
  }
  nrow(ppl)
}

  plotAnimals <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(likes(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(likes(a[1]), likes(a[2]), likes(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[1:2]), 
                            likes(a[2:3]), likes(a[c(1, 3)]), likes(a), euler.d = TRUE, ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]), 
                          likes(a[1:2]), likes(a[c(1, 3)]), likes(a[c(1, 4)]), likes(a[2:3]), 
                          likes(a[c(2, 4)]), likes(a[3:4]), likes(a[1:3]), likes(a[c(1, 2, 
                                                                                     4)]), likes(a[c(1, 3, 4)]), likes(a[2:4]), euler.d = TRUE, likes(a), ...)
  }
  if (!exists(out)) 
    out <- Oops
  return(out)
}
  
  pdf(paste("~/Documents/colleen_brachy/Data/", name, ".pdf", sep=""))
  plotAnimals(c("553.11 v IPO323" ,"560.11 v 553.11","560.11 v IPO323"),
            category = c("553.11 v IPO323" ,"560.11 v 553.11","560.11 v IPO323"), 
            alpha=0.5, scaled=F, 
            fill=c("lightgoldenrod1", "darkolivegreen1", "wheat3"),
            cat.fontfamily =rep("sans", 3),  cex = rep(1.5, 7),
            cat.dist =c(0.1, 0.1, 0.05),
            fontfamily = rep("sans", 7), cat.cex=rep(1.5, 3), margin=.15)
  dev.off()
}
dev.off()

for(i in unique(all_timepoint_sig$Isolate)){
  data<-all_timepoint_sig[(all_timepoint_sig$Isolate == i),]
  tab<- spread(as.data.frame(
    table(data$row, data$Timepoint)), key="Var2", value="Freq")
  assign(paste("table", i, sep="_"), tab)
  write.csv(tab, file=paste("~/Documents/colleen_brachy/Data/subset", i, ".csv", sep=""))

  }
{data = table_IPO323
  name<-"tableIPO323"
  likes <- function(animals) {
    ppl <- data
    names(ppl) <- colnames(data)
    for (i in 1:length(animals)) {
      ppl <- subset(ppl, ppl[animals[i]] == T)
    }
    nrow(ppl)
  }
  
  plotAnimals <- function(a, ...) {
    grid.newpage()
    if (length(a) == 1) {
      out <- draw.single.venn(likes(a), ...)
    }
    if (length(a) == 2) {
      out <- draw.pairwise.venn(likes(a[1]), likes(a[2]), likes(a[1:2]), ...)
    }
    if (length(a) == 3) {
      out <- draw.triple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[1:2]), 
                              likes(a[2:3]), likes(a[c(1, 3)]), likes(a), euler.d = TRUE, ...)
    }
    if (length(a) == 4) {
      out <- draw.quad.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]), 
                            likes(a[1:2]), likes(a[c(1, 3)]), likes(a[c(1, 4)]), likes(a[2:3]), 
                            likes(a[c(2, 4)]), likes(a[3:4]), likes(a[1:3]), likes(a[c(1, 2, 
                                                                                       4)]), likes(a[c(1, 3, 4)]), likes(a[2:4]), euler.d = TRUE, likes(a), ...)
    }
    if (!exists(out)) 
      out <- Oops
    return(out)
  }
  
  pdf(paste("~/Documents/colleen_brachy/Data/", name, ".pdf", sep=""))
  plotAnimals(c("4","9","21" ),
              category = c("4 DPI","9 DPI","21 DPI" ), 
              alpha=0.5, scaled=F, 
              fill=c("lightgoldenrod1", "darkolivegreen1", "wheat3"),
              cat.fontfamily =rep("sans", 3),  cex = rep(1.5, 7),
              cat.dist =c(0.1, 0.1, 0.05),
              fontfamily = rep("sans", 7), cat.cex=rep(1.5, 3), margin=.15)
  dev.off()
}
dev.off()


for(i in unique(all_timepoint_sig$Isolate)){
  data<-all_timepoint_sig[(all_timepoint_sig$Isolate == i),]
  tab<- as.data.frame(
    table(data$Timepoint, data$Species))
  tab$Isolate<-paste(i)
  assign(paste("species", i, sep="_"), tab)
  # write.csv(tab, file=paste("~/Documents/colleen_brachy/Data/subset", i, ".csv", sep=""))
}

for(t in unique(all_genotype_sig$Timepoint)){
  data<-all_genotype_sig[(all_genotype_sig$Timepoint == t),]
  tab<- as.data.frame(
    table(data$Comparison, data$Species))
  tab$Timepoint<-paste(t)
  assign(paste("species", t, sep="_"), tab)
  # write.csv(tab, file=paste("~/Documents/colleen_brachy/Data/subset", t, ".csv", sep=""))
}

species_iso<-rbind(species_553.11, species_560.11, species_IPO323)
colnames(species_iso)<-c("Timepoint", "Species", "DEGs", "Isolate")
species_time<-rbind(species_0, species_21, species_4, species_9)
colnames(species_time)<-c("Comparison", "Species", "DEGs", "Timepoint")
write.csv(species_iso, file="~/Documents/colleen_brachy/Data/DEGs_isolate_by_species.csv")
write.csv(species_time, file="~/Documents/colleen_brachy/Data/DEGs_timepoint_by_species.csv")


a<-ggplot(species_iso, aes(x=as.factor(Timepoint), y=DEGs)) +
  geom_col(aes(fill=Species)) + facet_wrap(~Isolate, ncol=1) +
  theme_bw() +
  scale_fill_manual(values=c("grey65", "grey50"), 
                    labels=c("Bd", "Zt"))+
  theme(text = element_text(size=15, colour="black"),
        axis.text.x = element_text(colour="black"), legend.position = "none") +
  xlab("Timepoint (DPI)")+
  ylab("Number of differentially expressed genes")


b<-ggplot(species_time, aes(x=as.factor(Timepoint), y=DEGs)) +
  geom_col(aes(fill=Species))  + facet_wrap(~Comparison, ncol=1) +
  theme_bw() +
  scale_fill_manual(values=c("grey65", "grey50"),
                    labels=c("Bd", "Zt")) +
  theme(text = element_text(size=15, colour="black"),
        axis.text.x = element_text(colour="black")) +
  xlab("Timepoint (DPI)")+
  ylab("Number of differentially expressed genes")


plot_grid(a,b, labels=c("A", "B"), ncol = 2, nrow = 1)
