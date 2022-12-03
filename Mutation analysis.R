### Read in TCGA SKCM CNV data (retrieve from cBioPortal) ###

mutation.data<-read.table("./R/Data/CNV_MC", sep="\t", header=T, row.names=1)



### log2 + 1 transform the gene expression data for consistency ###

mutation.data$IL15.mRNA.expression.RNA.Seq.V2.RSEM <- (log2(mutation.data$IL15.mRNA.expression.RNA.Seq.V2.RSEM) + 1)



### Remove gene expression values that are below zero or infinite ###

for(i in 1:nrow(mutation.data)){
    if(is.finite(mutation.data[i,2]) == FALSE){
        mutation.data<-mutation.data[-i,]
    }
}



### Wilcoxon test to compare expression levels in each CNV group ###

library(ggpubr); library(ggplot2)

compare_means(IL15.mRNA.expression.RNA.Seq.V2.RSEM ~ IL15.Putative.copy.number.alterations.from.GISTIC , data = mutation.data)



### Later on we want to include these Wilcoxon tests in a boxplot ###

my_comparisons <- list( c("Shallow Deletion", "Diploid"), c("Diploid", "Gain"), c("Diploid", "Amplification"))



### Boxplot of IL-15 expression levels in each CNV group ###

p <- ggplot(
    mutation.data, 
    aes(x=IL15.Putative.copy.number.alterations.from.GISTIC, 
    y=IL15.mRNA.expression.RNA.Seq.V2.RSEM)) +

geom_boxplot(
    alpha = 1.0, outlier.shape = NA, 
    aes(fill=IL15.Putative.copy.number.alterations.from.GISTIC), 
    width = 0.5)

p + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    axis.text.x = element_text(face="plain", color="black", size=8), 
    axis.text.y = element_text(face="plain", color="black", size=8), 
    axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
    legend.position = "none") +

stat_summary(
    geom = "errorbar", 
    width = 0.2,
    fun.min = function(z) quantile(z, 0), 
    fun = mean,
    fun.max = function(z) quantile(z, 1)) +

stat_compare_means(comparisons = my_comparisons, label = "p.format") +

scale_y_continuous(breaks=seq(0,10,2.5)) +

labs(
    x = "IL-15 copy number alteration", 
    y = "IL-15 expression") +

scale_x_discrete(
    labels=c("Diploid" = "Diploid\n(n=221)", "Amplification" = "Amplification\n(n=3)","Gain" = "Gain\n(n=48)", "Shallow Deletion" = "Shallow deletion\n(n=90)"),
    limits = c("Diploid", "Amplification", "Gain", "Shallow Deletion"))