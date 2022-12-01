### Read in TCGA SKCM RNAseq data

rna.data.interactions<- read.table("./R/Data/HiSeqV2_MC",sep="\t",head=TRUE,row.names=1)
rna.data.interactions <- as.matrix(rna.data)



### Identify the genes whose expression are most closely associated with IL-15.  ###

library(limma)

cor.design <- cbind(intercept=1,IL15=rna.data.interactions["IL15",]) 
cor.fit <- lmFit(rna.data.interactions[,],design=cor.design) 
cor.fit <- eBayes(cor.fit) 

IL15.100cor.genes <- rownames(topTable(cor.fit,coef=2,number=100))
write.table(IL15.100cor.genes,file="~/VScode/R/Results/MC IL15 correlated100genes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE) # analyse on https://david.ncifcrf.gov/tools.jsp and https://string-db.org/



### After analysing the above results on STRING, read in proteins that are predicted to directly interact with IL-15 (CD45 and CD69) ###

IL15.network.nodes<-read.table("./R/Data/IL15_MC_interactions_names", sep="\t", header=T)



### Extract RNAseq data for CD45 and CD69 ###

rna.data.interactions<-rna.data[which(is.element(rownames(rna.data),IL15.network.nodes[,1])),] 



### Extract information on CD45 and CD69 for IL-15 high samples only ###

IL15.high <- as.matrix((rna.data["IL15",]>median(rna.data["IL15",])))
IL15.high.samples <- matrix(which(IL15.high=="TRUE"))
IL15.high.samples <- colnames(rna.data.interactions[,IL15.high.samples])
interactions.IL15.high <- rna.data.interactions[,IL15.high.samples]



### Store above information in format that is accessible for downstream use ###

library(tidyverse)

interactions.IL15.high <- as_tibble(t(interactions.IL15.high))
interactions.IL15.high <- interactions.IL15.high %>% pivot_longer(
    cols = c("PTPRC","CD69"),
    names_to = "Gene",
    values_to = "Gene_expression"
) 



### Repeat last two steps for IL-15 low samples ###

IL15.low.samples <- matrix(which(IL15.high=="FALSE"))
IL15.low.samples <- colnames(rna.data.interactions[,IL15.low.samples])
interactions.IL15.low <- rna.data.interactions[,IL15.low.samples]

interactions.IL15.low <- as_tibble(t(interactions.IL15.low))
interactions.IL15.low <- interactions.IL15.low %>% 
    pivot_longer(
        cols = c("PTPRC","CD69"),
        names_to = "Gene",
        values_to = "Gene_expression"
    ) 


### Store the above two tibbles in a list ###

l <- list(interactions.IL15.high,interactions.IL15.low)
names(l) <- c("High","Low")


### This will allow you to make a single boxplot from two tibbles ###

library(data.table)

interactions.IL15.all<- rbindlist(l, id="id")



### This will be used to colour the following boxplot by IL-15 expression level ###

colors <- c("Low" = "red", "High" = "blue")



### Boxplot of gene expression vs CD45 and CD69 gene startified by IL-15 expression level ###

library(ggpubr); library(ggplot2)

p <- ggplot(
    interactions.IL15.all, 
    aes(x=Gene, 
    y=Gene_expression, 
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
    y="Gene expression", 
    fill = "IL-15 group") +

 coord_cartesian(ylim = c(0,16)) +
 
 stat_summary(
    geom = "errorbar", width = 0.2,
    fun.min = function(z) quantile(z, 0), 
    fun = mean,
    fun.max = function(z) quantile(z, 1), position = position_dodge(width = 0.8)) +

stat_compare_means(
    label = "p.format", 
    label.y = 16) +

scale_x_discrete(
    labels=c("PTPRC" = "CD45","CD69" = "CD69"), 
    limits = c("PTPRC", "CD69")) 