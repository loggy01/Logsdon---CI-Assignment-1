### Read in TCGA SKCM miRNA and RNAseq data

miRNA.data <- read.table("./R/Data/miRNA_MC", header = TRUE, sep = "\t", row.names = 1)

rna.data.miRNA <- read.table("./R/Data/HiSeqV2_MC",sep="\t",head=TRUE,row.names=1)
rna.data.miRNA <- as.matrix(rna.data.miRNA)



### Ensure same samples for RNAseq and miRNA data ###

rna.data.miRNA<-rna.data.miRNA[,which(is.element(colnames(rna.data.miRNA),colnames(miRNA.data)))] 
miRNA.data2<-miRNA.data[,which(is.element(colnames(miRNA.data),colnames(rna.data.miRNA)))] 



### Order RNAseq and methylation data the same as each other ###

rna.data.miRNA<-as.matrix(rna.data.miRNA[,order(colnames(rna.data.miRNA))])
miRNA.data2<-as.matrix(miRNA.data2[,order(colnames(miRNA.data2))])



### Remove rows with > 50% NAs ###

na.count<-apply(miRNA.data2,1,function(x) sum(as.numeric(is.na(x))))
exclude<-as.numeric(na.count>0.5*ncol(miRNA.data2))
miRNA.data2<-miRNA.data2[which(exclude==0),]



### Construct a array to store miRNAs with expression correlated to IL-15 expression ###

results<-array(NA,c(nrow(miRNA.data2),3))
rownames(results)<-rownames(miRNA.data2)
colnames(results)<-c("Correlation","P value","Adjusted P value")



### Populate above table ###

for (i in 1:nrow(miRNA.data2)) {
results[i,1]<-cor.test(as.numeric(rna.data.miRNA["IL15",]),as.numeric(miRNA.data2[i,]), use="c")$est
results[i,2]<-cor.test(as.numeric(rna.data.miRNA["IL15",]),as.numeric(miRNA.data2[i,]), use="c")$p.value
}



### Order results by adjusted p value. I selected MIMAT0000691","MIMAT0000243","MIMAT0000073","MIMAT0000074 ###

results <- as.data.frame(results)
results[,3]<-p.adjust(results[,2],method="fdr")
results<-results[order(results[,1], decreasing=F),]



### Subset information on chosen miRNAs ###

miRNA.hits <- results[c("MIMAT0000691","MIMAT0000243","MIMAT0000073","MIMAT0000074"),]
miRNA.IL15 <- miRNA.data2[rownames(miRNA.hits),]



### Remove NAs ###

col.remove <- c()

for(i in 1:nrow(miRNA.IL15)){
    for(j in 1:ncol(miRNA.IL15)){
        if(is.na(miRNA.IL15[i,j]) == TRUE){
            col.remove<-c(col.remove, j)
            next 
        }
    }
}

col.remove <- colnames(miRNA.IL15[,col.remove])

miRNA.IL15 <- miRNA.IL15[,which(!is.element(colnames(miRNA.IL15),col.remove))]



### Make sure RNAseq data and new miRNA are ordered the same ###

rna.data.miRNA2<-rna.data.miRNA[,which(is.element(colnames(rna.data.miRNA),colnames(miRNA.IL15)))] 

rna.data.miRNA2<-as.matrix(rna.data.miRNA2[,order(colnames(rna.data.miRNA2))])
miRNA.IL15<-as.matrix(miRNA.IL15[,order(colnames(miRNA.IL15))])



### Extract information on miRNAs for IL-15 high samples only ###

IL15.high <- as.matrix((rna.data.miRNA2["IL15",]>median(rna.data.miRNA2["IL15",])))
IL15.high.samples <- matrix(which(IL15.high=="TRUE"))
IL15.high.samples <- colnames(rna.data.miRNA2[,IL15.high.samples])
miRNA.IL15.high <- miRNA.IL15[,IL15.high.samples]



### Store above information in format that is accessible for downstream use ###

library(tidyverse)

miRNA.IL15.high <- as_tibble(t(miRNA.IL15.high))
miRNA.IL15.high <- miRNA.IL15.high %>% 
    pivot_longer(
        cols = c("MIMAT0000691","MIMAT0000243","MIMAT0000073","MIMAT0000074"),
        names_to = "miRNA",
        values_to = "miRNA_expression"
    ) 



### Repeat last two steps for IL-15 low samples ###

IL15.low.samples <- matrix(which(IL15.high=="FALSE"))
IL15.low.samples <- colnames(rna.data.miRNA2[,IL15.low.samples])
miRNA.IL15.low <- miRNA.IL15[,IL15.low.samples]

miRNA.IL15.low <- as_tibble(t(miRNA.IL15.low))
miRNA.IL15.low <- miRNA.IL15.low %>% 
    pivot_longer(
        cols = c("MIMAT0000691","MIMAT0000243","MIMAT0000073","MIMAT0000074"),
        names_to = "miRNA",
        values_to = "miRNA_expression"
    ) 


### Store the above two tibbles in a list ###

l <- list(miRNA.IL15.high,miRNA.IL15.low)
names(l) <- c("High","Low")



### This will allow you to make a single boxplot from two tibbles ###

library(data.table)

miRNA.IL15.all<- rbindlist(l, id="id")



### This will be used to colour the following boxplot by IL-15 expression level ###

colors <- c("Low" = "red", "High" = "blue")



### Boxplot of miRNA expression vs miRNA gene startified by IL-15 expression level ###

library(ggpubr); library(ggplot2)

p <- ggplot(
    miRNA.IL15.all, 
    aes(x=miRNA, y=miRNA_expression, 
    fill=interaction(id))) + 

geom_boxplot(
    alpha = 1, 
    width = 0.5, 
    position = position_dodge(width = 0.8),
    outlier.shape = NA) 
  
p + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.text = element_text(colour="black", size=11, face="plain"),
    panel.background = element_blank(),
    axis.text.x = element_text(face="plain", color="black", size=8), 
    axis.text.y = element_text(face="plain", color="black", size=8), 
    axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +

scale_fill_manual(values = colors) +

labs(
    y="miRNA expression", 
    fill = "IL-15 group") +

coord_cartesian(ylim = c(0,22)) +

stat_summary(
    geom = "errorbar", width = 0.2,
    fun.min = function(z) quantile(z, 0), 
    fun = mean,
    fun.max = function(z) quantile(z, 1), position = position_dodge(width = 0.8)) +

stat_compare_means(
    label = "p.format", 
    label.y = 20) +

scale_y_continuous(breaks=seq(0,20,5)) +

scale_x_discrete(labels=c("MIMAT0000073" = "miR-19a-3p", "MIMAT0000074" = "miR-19b-3p","MIMAT0000243" = "miR-148a-3p", "MIMAT0000691" = "miR-130b-3p"))