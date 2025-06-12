#Sparse Projection to Latent Structure - Discriminant Analysis for Metabolites
library(mixOmics)

#PARAMETERS
abund_df<-read.csv("../data/metabolite_abundance.csv",header=TRUE,row.names=1)
abund_df[abund_df$Treatment=="No","Treatment"] <- "T1"
abund_df[abund_df$Treatment=="Yes","Treatment"] <- "T2"
abund_df$Treatment <- as.factor(abund_df$Treatment)
# meta_df <- data.frame(c(rep("T1",10), rep("T2",10)))
# rownames(meta_df) <- paste0("MP", 1:20)


#Perform PLS-DA following steps from http://mixomics.org/mixmc/case-study-hmp-bodysites-repeated-measures/
#Step 1: 
abund_df.plsda<-plsda(X=abund_df[-1], Y=abund_df$Treatment, ncomp = 10)
abund_df.plsda.perf<-perf(abund_df.plsda,validation="Mfold",folds=5,progressBar=TRUE,auc=TRUE,nrepeat=10)

# Balanced Error Rate plot
pdf("Step.1-PLS-DA_performance.pdf")
plot(abund_df.plsda.perf, overlay = 'measure', sd = TRUE)
dev.off()

#Step 2: Plot PLS-DA
pdf("Step.2-PLS-DA_plotIndiv.pdf")
plotIndiv(abund_df.plsda,comp=1:2,group=abund_df$Treatment,ind.names=FALSE,ellipse=TRUE,legend=TRUE,title="PLS-DA comp 1-2")
dev.off()

#Step 3: tuning PLS-DA
abund_df.plsda.tune<-tune.splsda(X=abund_df[-1], 
                                    Y=abund_df$Treatment,
                                    ncomp=3,
                                    test.keepX =c(seq(20,200,2)),
                                    validation="Mfold",
                                    folds=5,
                                    progressBar=TRUE,
                                    dist="max.dist",
                                    measure="overall",
                                    nrepeat=10)

pdf("Step.3-sPLS-DA_tuning.pdf")
plot(abund_df.plsda.tune)
dev.off()

#Step 4:
#Choose optimal number of variables to select on tuning_components comps:
select.keepX = abund_df.plsda.tune$choice.keepX[1:3]
#Now we run sPLS-DA multilevel analysis on the selected variables
abund_df.splsda<-splsda(X=abund_df[-1], 
                        Y=abund_df$Treatment, ncomp = 3,keepX = select.keepX)
pdf("Step.4-sPLS-DA_plotIndiv.pdf")
plotIndiv(abund_df.splsda,comp=1:2,group=abund_df$Treatment,ind.names=FALSE,ellipse=TRUE,legend=TRUE,title="sPLS-DA comp 1-2")
dev.off()

#Step 5:
abund_df.splsda.perf<-perf(abund_df.splsda,validation="Mfold",folds=5,progressBar=FALSE,auc=TRUE,nrepeat=10)
pdf("Step.5-sPLS-DA_performance.pdf")
plot(abund_df.splsda.perf, overlay = 'measure', sd = TRUE)
dev.off()

#Step 6:
selectVar(abund_df.splsda, comp = 1)$value
selected.features.comp1 = selectVar(abund_df.splsda, comp = 1)$name
selected.features.comp2 = selectVar(abund_df.splsda, comp = 2)$name
selected.features.comp3 = selectVar(abund_df.splsda, comp = 3)$name
#Step 7:
pdf("Step.6-Plot_loadings.pdf")
plotLoadings(abund_df.splsda,comp=1,contrib="max",method="mean")
dev.off()

#Step 8:
#A heatmap will also help understanding the selected features. 
pdf("Step.8-sPLS-DA_heatmap.pdf",height=10,width=14)
cim(abund_df.splsda,row.sideColors = c(rep("#FFA8BB",10), rep("#5EF1F2",10)),margins=c(20,6))
dev.off()



