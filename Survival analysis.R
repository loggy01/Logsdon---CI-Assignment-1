### Read in RNAseq, survival, and clinical TCGA datasets for SKCM ###

rna.data <- read.table("./R/Data/HiSeqV2_MC",sep="\t",head=TRUE,row.names=1)
rna.data <- as.matrix(rna.data)

surv.data <- read.table("./R/Data/survival_MC", sep="\t", header=T,row.names=1)  
rownames(surv.data)<-gsub(rownames(surv.data), pattern="-", replace=".")

clin.data <- read.table("./R/Data/clinical_MC", sep="\t", header=T, row.names=1)
rownames(clin.data)<-gsub(rownames(clin.data), pattern="-", replace=".")



### Generate survival object from TCGA SKCM patients ###

library(survival)

os.time <- surv.data[colnames(rna.data),"OS.time"]/365.25
os.event <- as.numeric(surv.data[colnames(rna.data),"OS"])
MC.os <- Surv(os.time,os.event)



### Observe the clinical features for SKCM to select the most relevant ###

clin.data<-clin.data[colnames(rna.data),]



### Univariate Cox regression analyses of age, disease stage, radiotherapy, ulceration, tumour size respectively###

age <- as.numeric(clin.data$age_at_initial_pathologic_diagnosis)
age[is.na(age)] <- 0
age.high<-rep(0,nrow(clin.data))

for (i in 1:nrow(clin.data)) {
    if (age[i] >= 65) {
        age.high[i] = 1
    }
}

age.uni.coxph <- coxph(MC.os ~ age.high) # old age: > 65



x1<-grep("III",clin.data$pathologic_stage)
x2<-grep("IV",clin.data$pathologic_stage)
stage.high<-rep(0,nrow(clin.data))
stage.high[c(x1,x2)]<-1

stage.uni.coxph <- coxph(MC.os ~ stage.high) # high stage: III and IV



x3<-grep("YES",clin.data$radiation_therapy)
radiation.positive<-rep(0,nrow(clin.data))
radiation.positive[x3]<-1

radiation.uni.coxph <- coxph(MC.os ~ radiation.positive)



x4<-grep("YES",clin.data$melanoma_ulceration_indicator)
ulceration.positive<-rep(0,nrow(clin.data))
ulceration.positive[x4]<-1

ulceration.uni.coxph <- coxph(MC.os ~ ulceration.positive)



x5<-grep("T0",clin.data$pathologic_T)
x6<-grep("T1",clin.data$pathologic_T)
x7<-grep("T2",clin.data$pathologic_T)
tumour.small<-rep(0,nrow(clin.data))
tumour.small[c(x5,x6,x7)]<-1

tumour.uni.coxph <- coxph(MC.os ~ tumour.small) # small tumour: T0, T1, T2



### Construct a dataframe to store genes determined to be most associated with survival based on the above clinical features ###

results.multivariate<-array(NA, c(nrow(rna.data),4)) 
colnames(results.multivariate)<-c("HR","LCI","UCI","PVAL") 
rownames(results.multivariate)<-rownames(rna.data)
results.multivariate<-as.data.frame(results.multivariate)



### Filling in the above table by generating multivariate Cox regression models for each gene with the above clinical features ###

for(i in 1:nrow(rna.data)){
 coxphmodel <- coxph(MC.os ~ rna.data[i,]+age.high+stage.high+radiation.positive+ulceration.positive+tumour.small)
 
 results.multivariate$HR[i] <- summary(coxphmodel)$coef[1,2]
 results.multivariate$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
 results.multivariate$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
 results.multivariate$PVAL[i] <- summary(coxphmodel)$coef[1,5]
}



### Order genes by adjusted p value and select a gene with an interesting hazard ratio. I have selected IL-15 ###

results.multivariate<-as.data.frame(results.multivariate)
results.multivariate$FDR<-p.adjust(results.multivariate$PVAL,method="fdr")
results.multivariate<-results.multivariate[order(results.multivariate$FDR, decreasing=F),]

results.multivariate[1:10,]



### Univariate Cox regression analysis of IL-15 expression level ###

IL15.high <- as.numeric(rna.data["IL15",]>median(rna.data["IL15",]))

IL15.uni.coxph <- coxph(MC.os ~ IL15.high)



### Multivariate Cox regression analysis with IL-15 expression level and previously determined clinical features ###

coxphmodel <- coxph(MC.os ~ IL15.high+age.high+stage.high+radiation.positive+ulceration.positive+tumour.small)

test.coxphmodel <- cox.zph(coxphmodel)



### Log rank test for TCGA SKCM patient survival startified by IL-15 expression level ###

surv.diff <- survdiff(MC.os ~ IL15.high)

surv.stats <- survfit(MC.os ~ IL15.high)



### Construct data table for univariate and multivariate Cox regression analyses of age, disease stage, radiotherapy, ulceration, tumour size, and IL-15 expression in the TCGA SKCM cohort ###

library(dplyr)

unicox.plot.data <- tibble::tibble(mean  = c(1.86,1.74,0.52,2.35,0.56,0.51),
                            lower = c(1.39,1.32,0.31,1.76,0.42,0.39),
                            upper = c(2.48,2.31,0.86,3.15,0.76,0.68),
                            prognostic = c("Old age","High stage","Radiotherapy received","Ulceration present","Small tumour size","High IL-15 expression"),
                            HR = c("1.86(1.39-2.48)", "1.74(1.32-2.31)", "0.52(0.31-0.86)", "2.35(1.76-3.15)", "0.56(0.42-0.76)", "0.51(0.39-0.68)"),
                            Pval = c("2.94e-05", "9.78e-05", "1.12e-02", "9.88e-09", "1.62e-04", "1.79e-06"))

multicox.plot.data <- tibble::tibble(mean  = c(1.47,1.79,0.54,1.77,0.73,0.54),
                            lower = c(1.09,1.35,0.32,1.30,0.53,0.41),
                            upper = c(1.99,2.38,0.91,2.40,1.01,0.71),
                            prognostic = c("Old age","High stage","Radiotherapy received","Ulceration present","Small tumour size","High IL-15 expression"),
                            HR = c("1.47(1.09-1.99)","1.79(1.35-2.38)","0.54(0.32-0.91)","1.77(1.30-2.40)","0.73(0.53-1.01)","0.54(0.41-0.71)"),
                            Pval = c("1.07e-02","5.72e-05","2.01e-02","2.77e-04","5.43e-02","1.57e-05"))



### Function for generating a forest plot ###

library(forestplot)

forest.plotter <- function (cox.data,colour) {
    cox.data %>% 

        forestplot(labeltext = c(prognostic,HR,Pval),
        xticks = c(0, 1, 2, 3),
        lwd.x = 1.5,
        zero = 1,
        lwd.zero = 1.5,
        ci.vertices = TRUE,
        ci.vertices.height = 0.05,
        lwd.ci = 1.5,
        boxsize = 0.3,
        xlab = "Hazard ratio") %>%

            fp_set_style(box = colour,
            line = "blue",
            zero = "black",
            txt_gp = fpTxtGp(ticks=gpar(cex=0.65), xlab=gpar(cex=0.9), label = gpar(cex = 0.9))) %>% 

                fp_add_header(prognostic = c("", ""),
                HR = c("", "Hazard ratio"),
                Pval = c("", "P value"))
}



### Generate forest plots for univariate and multivariate Cox regression analyses of age, disease stage, radiotherapy, ulceration, tumour size, and IL-15 expression in the TCGA SKCM cohort ###

forest.unicox.plot<-forest.plotter(unicox.plot.data,"green")
forest.multicox.plot<-forest.plotter(multicox.plot.data,"orange")



### Generate Kaplan-Meier survival curves for TCGA SKCM patients stratified by IL-15 expression level ###

library(survminer) 

ggsurvplot(
    survfit(MC.os ~ IL15.high), 
    data = clin.data, 
    palette = c("red", "blue"),
    conf.int = TRUE,
    ggtheme = theme_survminer(font.main = 11,font.x = 11, font.y = 11, font.legend = 11,font.tickslab = c(8, "black")),
    size = 0.7,
    surv.median.line = "hv",
    legend = c(0.9,0.8),
    legend.title = "IL-15 expression",
    legend.labs = c("Low", "High"),
    xlab = "Time (years)",
    ylab = "Survival probability",
    xlim = c(0, 32),
    break.time.by = 2,
    pval = "p = 1e-06",
    pval.size = 3.75,
    pval.coord = c(0, 0.1),
    risk.table = TRUE,
    fontsize = 3,
    risk.table.col = "strata",
    tables.height = 0.35,
)