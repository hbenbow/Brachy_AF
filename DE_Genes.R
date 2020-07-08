library(DESeq2)
library(ggplot2)
library(pivottabler)
library(VennDiagram)
genes<-read.csv("~/Documents/colleen_brachy/Data/all_genes.csv")
all_genotype_sig<-read.csv("~/Documents/colleen_brachy/Data/all_sig_by_genotype.csv")
all_timepoint_sig<-read.csv("~/Documents/colleen_brachy/Data/all_sig_by_timepoint.csv")

# make some frequency tables, from which we can do venn diagrams (childish but useful)
for(t in unique(all_genotype_sig$Timepoint)){
  data<-all_genotype_sig[(all_genotype_sig$Timepoint == t),]
  tab<- spread(as.data.frame(
    table(data$row, data$Comparison)), key="Var2", value="Freq")
  assign(paste("table", t, sep="_"), tab)
}
  
{data = table_0
name<-"table0"
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
