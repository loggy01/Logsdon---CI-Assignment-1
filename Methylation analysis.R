### Read in 450k annotation from Illumina, TCGA SKCM methylation data, and TCGA SKCM RNAseq data

meth.annot<-readRDS("./R/Data/annot450k.rds")
meth.data<-read.table("./R/Data/HumanMethylation450_MC", sep="\t", header=T, row.names=1)

rna.data.meth <- read.table("./R/Data/HiSeqV2_MC",sep="\t",head=TRUE,row.names=1)
rna.data.meth <- as.matrix(rna.data.meth)



### Ensure same samples for RNAseq and methylation data ###

rna.data.meth<-rna.data.meth[,which(is.element(colnames(rna.data.meth),colnames(meth.data)))] 
meth.data2<-meth.data[,which(is.element(colnames(meth.data),colnames(rna.data.meth)))] 


### Order RNAseq and methylation data the same as each other ###

rna.data.meth<-as.matrix(rna.data.meth[,order(colnames(rna.data.meth))])
meth.data2<-as.matrix(meth.data2[,order(colnames(meth.data2))])



### Extract IL15 methylation annotation and apply it to our methylation data ###

meth.IL15<-rownames(meth.annot[which(meth.annot$UCSC_RefGene_Name=="IL15"),])
meth.data.IL15<-meth.data2[meth.IL15,]



### Remove rows with > 50% NAs ###

na.count<-apply(meth.data.IL15,1,function(x) sum(as.numeric(is.na(x))))
exclude<-as.numeric(na.count>0.5*ncol(meth.data.IL15))
meth.data.IL15<-meth.data.IL15[which(exclude==0),]


### Construct a array to store IL-15 gene sites where methylation is correlated with expression ###

results<-array(NA,c(nrow(meth.data.IL15),3))
rownames(results)<-rownames(meth.data.IL15)
colnames(results)<-c("Correlation","P value","Adjusted P value")



### Populate above table ###

for (i in 1:nrow(meth.data.IL15)) {
results[i,1]<-cor.test(as.numeric(rna.data.meth["IL15",]),as.numeric(meth.data.IL15[i,]), use="c")$est
results[i,2]<-cor.test(as.numeric(rna.data.meth["IL15",]),as.numeric(meth.data.IL15[i,]), use="c")$p.value
}



### Order results by adjusted p value ###

results <- as.data.frame(results)
results[,3]<-p.adjust(results[,2],method="fdr")
results<-results[order(results[,3], decreasing=F),]



### Construct a dataframe storing IL-15 expression and methylation data, grouped by IL-15 expression level, for later use ###

expression.IL15 <- as.matrix(rna.data.meth["IL15",])
methylation.IL15 <- as.matrix(apply(meth.data.IL15,2,median,na.rm=T))
IL15.high3 <- as.factor(rna.data.meth["IL15",]>median(rna.data.meth["IL15",]))
methylation.scatter <- data.frame(methylation.IL15, expression.IL15,IL15.high3)
levels(methylation.scatter$IL15.high3) <- c("Low", "High")



### This will be used to colour a following scatterplot by IL-15 expression level ###

colors <- c("Low" = "red", "High" = "blue")



### Wilcox test to compare IL-15 methylation levels in high and low IL-15 expression groups ###

low.group <- subset(methylation.scatter, IL15.high3 == "Low")
high.group <- subset(methylation.scatter, IL15.high3 == "High")
wilcox.test(high.group$methylation.IL15, low.group$methylation.IL15)



### Plot scatterplot of IL-15 expression level vs. IL-15 methylation level ###

library(ggpubr); library(ggplot2)

p1 <- ggplot(
    methylation.scatter, 
    aes(x=methylation.IL15, y=expression.IL15)) + 
    
geom_point(aes(colour=IL15.high3), size = 0.5)

p1 + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.position = c(0.80, 0.875), 
    legend.text = element_text(colour="black", size=11, face="plain"),
    panel.background = element_blank(),
    axis.text.x = element_text(face="plain", color="black", size=8),
    axis.text.y = element_text(face="plain", color="black", size=8), 
    axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +

scale_color_manual(values = colors) +

coord_cartesian(
    xlim =c(0,1), 
    ylim = c(0,10)) +

scale_y_continuous(breaks=seq(0,10,2.5)) +

labs(color = "IL-15 group") +

xlab("IL-15 methylation") + 
ylab("IL-15 expression") +

geom_vline(
    xintercept = 0.069, 
    size = 0.5, 
    linetype = "dashed") + 

geom_hline(
    yintercept = 5.0, 
    size = 0.5, 
    linetype = "dashed") +

annotate("text", x=0.825, y=7, label= "p = 8.67e-05") 

