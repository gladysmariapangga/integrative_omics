#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for Data Integration Analysis for Biomarker Discovery using latent variable approaches for 'Omic studies (DIABLO) 

library(vegan)
library(ggplot2)
library(ape)
library(phangorn)
library(stringr)
library(reshape2)
library(mixOmics)

#PARAMETERS ###########################
abund_df<-read.csv("../data/taxonomy_abundance_normalised.csv",header=TRUE,row.names=1)

met_df<-read.csv("../data/metabolite_abundance.csv",header=T,row.names=1,check.names=FALSE)
met_df[met_df$Treatment=="No","Treatment"] <- "T1"
met_df[met_df$Treatment=="Yes","Treatment"] <- "T2"
met_df$Treatment <- as.factor(met_df$Treatment)
rownames(met_df) <- c("MP01", "MP02", "MP03", "MP04", "MP05",  
                                 "MP06", "MP07", "MP08",  "MP09",  "MP10",
                                 "MP11", "MP12", "MP13", "MP14", "MP15",
                                 "MP16", "MP17","MP18","MP19","MP20")


#For DIABLO, we need to combine our data
data=list(metagenome=abund_df,metabolome=met_df[-1])
#Now make a design matrix where all blocks (datasets) are connected with a link of 0.1
design=matrix(0.1,ncol=2,nrow=2,dimnames=list(names(data),names(data)))


#Perform DIABLO following steps from http://mixomics.org
#Step 1:
data.diablo=block.splsda(X=data,Y=met_df$Treatment,ncomp=10,design=design)
data.diablo.perf<-perf(data.diablo,validation="loo",folds=5,progressBar=TRUE,auc=TRUE,nrepeat=10)

#The plot indicates a decrease in the classification error rate. The BER stands for Balanced Error Rate 
#and should be considered when we have an unbalanced number of samples per group. 
pdf(paste("Step.1-DIABLO_performance_",label,".pdf",sep=""), height=4,width=10)
plot(data.diablo.perf, overlay = 'measure', sd = TRUE)
dev.off()

#Step 2:
# From the performance plot, we can look at both overall and
# balanced error rate (BER) to see which distance matrix gives the best accuracy
# you may have to change it later in the paramters
# This will give us the number of components to choose

data.diablo.tune=tune.block.splsda(X=data,
                                   Y=met_df$Treatment,
                                   ncomp=2,
                                   design=design,
                                   test.keepX=list(metagenome=c(seq(10,50,5)),metabolome=c(seq(50,500,50))),
                                   validation="loo",
                                   folds=5,
                                   measure="overall",
                                   nrepeat=10,
                                   dist="centroids.dist")

list.keepX=data.diablo.tune$choice.keepX

#Now run the final DIABLO modle
data.diablo=block.splsda(X=data,Y=met_df$Treatment,ncomp=ncomp,keepX=list.keepX,design=design)
Y=met_df$Treatment
for(i in 1:ncomp){
  pdf(paste("Step.2.",i,"-DIABLO_correlation.comp",i,"_",label,".pdf",sep=""),width=8, height=8)
  plotDiablo(data.diablo,ncomp=i)
  dev.off()
}

#Step 3: Now use plotIndiv 
pdf(paste("Step.3-DIABLO_plotIndiv_",label,".pdf",sep=""), width=10, height=4.5)
plotIndiv(data.diablo,group=met_df$Treatment,col=colours[1:length(levels(met_df$Treatment))],ind.names=FALSE,ellipse=draw_ellipse,legend=TRUE)
dev.off()

#Step 5: Get significant OTUs
write.csv(data.diablo$loadings[[1]][abs(rowSums(data.diablo$loadings[[1]]))>0,],paste("Step.5-DIABLO_significant_OTUs_",paste(levels(met_df$Treatment),collapse="_"),"_",label,".csv",sep=""))
write.csv(data.diablo$loadings[[2]][abs(rowSums(data.diablo$loadings[[1]]))>0,],paste("Step.5-DIABLO_significant_Mtbs_",paste(levels(met_df$Treatment),collapse="_"),"_",label,".csv",sep=""))


# metabolite_loadings_df <- as.data.frame(data.diablo$loadings$metabolome)
# metabolite_loadings_df$average <- rowMeans(metabolite_loadings_df[1:3])
# metabolite_loadings_df$Sign <- ifelse(metabolite_loadings_df$average==0,NA,ifelse(metabolite_loadings_df$average>0,"T1","T2"))
# write.csv(metabolite_loadings_df,"Step.5-DIABLO_metabolomeloadings.csv")

#Step 6: Now use plotVar 
#The correlation circle plot highlights the contribution of each selected
#variable to each component. Important variables should be close to the large circle
#plotVar displays the variables from all blocks, selected on component 1 and 2
#Clusters of points indicate a strong correlation between variables.
pdf(paste("Step.6-DIABLO_plotVar_",label,".pdf",sep=""))
plotVar(data.diablo, var.names = plot_names_correlation_circle, pch=rep(15,length(data)),cex=rep(size_features_correlation_circle,length(data)),style = 'graphics', legend = TRUE)
dev.off()

#Step 7: Now use circosPlot
#The circos plot represents the correlations between variables of different
#types, represented on the side quadrants. Several display options
#are possible, to show within and between connexions between blocks, 
#expression levels of each variable according to each class (argument line=TRUE)
#The circos plot is built based on similarity matrix, which was exteded
#to the case of multiple 

pdf(paste("Step.7-DIABLO_circosPlot",label,".pdf",sep=""))
circosPlot(data.diablo,cutoff=0.7,
           line=circos_plot_show_expression,
           color.Y=colours[1:length(levels(met_df$Treatment))],
           color.cor=c("red","blue"),
           showIntraLinks=circos_plot_show_intralinks_correlation,
           color.blocks = alpha(c("green", "darkgreen"), .2),
           size.variables = 0.2)
dev.off()



#color.Y is T1 and T2
#var.adj up and down adjustment
#labels are metagenome/metabolome
#variables are the metabolites, mags

####################Get correlation heatmap####################
#Reference: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cp<-circosPlot(data.diablo,cutoff=circos_plot_correlation_cutoff,
               line=circos_plot_show_expression,
               color.Y=colours[1:length(levels(met_df$Treatment))],
               #Temporary bugfix is to get rid of colors all together
               color.blocks=rep(NA,2),
               color.cor=c("red","blue"),
               showIntraLinks=circos_plot_show_intralinks_correlation,
               size.labels=circos_plot_label_size,
               size.variables=circos_plot_variable_size)
cormat<-cp

colnames(cormat)<-as.character(sapply(colnames(cormat),function(x){if(grepl("OTU_",x)){gsub(";+$","",paste(x,as.character(paste(OTU_taxonomy[x,],collapse=";")),sep=":"))}else{print(x)}}))
rownames(cormat)<-as.character(sapply(rownames(cormat),function(x){if(grepl("OTU_",x)){gsub(";+$","",paste(x,as.character(paste(OTU_taxonomy[x,],collapse=";")),sep=":"))}else{print(x)}}))

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#Now only retain the strongest correlations
in_interval <- function(x, interval){
  stopifnot(length(interval) == 2L)
  interval[1] < x & x < interval[2]
}
melted_cormat<-melted_cormat[!in_interval(melted_cormat$value, c(-0.6,0.6)),]

#Get rid of same names
melted_cormat<-melted_cormat[!sapply(seq(1:nrow(melted_cormat)),function(x){if(as.character(melted_cormat[x,1])==as.character(melted_cormat[x,2])) TRUE else FALSE}),]
# write.csv(melted_cormat, "melted_cormat.csv", row.names=TRUE)

melted_cormat$Var1 <- as.character(melted_cormat$Var1)
melted_cormat$Var2 <- as.character(melted_cormat$Var2)
# #genome x metabolome
Ind<-sapply(rownames(melted_cormat),function(x){!sum(grepl("^C",melted_cormat[x,c("Var1","Var2")]))==2})

melted_cormat2 <- melted_cormat[Ind,]
Ind2<-sapply(rownames(melted_cormat),function(x){sum(grepl("^G",melted_cormat[x,c("Var1","Var2")]))==1})
melted_cormat2 <- melted_cormat[Ind2,]
write.csv(melted_cormat2, "melted_cormat2.csv", row.names=TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90,vjust=1,hjust=1, size = 2))+
  theme(axis.text.y = element_text(vjust = 1,
                                  size = 2, hjust = 1))
  coord_fixed()

#ggheatmap<-ggheatmap + 
#  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) 
ggheatmap<-ggheatmap +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.border = element_blank(),
    # panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
# Print the heatmap
pdf(paste("Step.7-DIABLO_CorrelationHeatmap_",label,".pdf",sep=""),height=20,width=20)
print(ggheatmap)
dev.off()
write.csv(melted_cormat,paste("Step.7-DIABLO_CorrelationHeatmap_",label,".csv",sep=""))



#Step 8: 
#The plotLoadings shows that all features selected on the first component
for (i in 1:ncomp){
  pdf(paste("Step.8.",i,"-DIABLO_loading.comp",i,"_",label,".pdf",sep=""),height=30,width=23)
  plotLoadings(data.diablo,comp=i,legend.color=colours[1:length(levels(met_df$Treatment))],contrib=loading_contrib,method=loading_method,title=paste('DIABLO Loadings (comp',i,')',sep=""))
  dev.off()
}

#Step 9:
#A heatmap will now specifically represent the multi-'omics molecular signature expression
#for each sample
pdf(paste("Step.9-DIABLO_heatmap_",label,".pdf",sep=""),height=23,width=35)
cimDiablo(data.diablo, color.Y=colours[1:length(levels(met_df$Treatment))],margins=c(80,40))
dev.off()


#Step 10
#Get all pairwise-correlations and save them in the form of a CSV file##################

tables_list<-names(data)
discriminant_variables<-NULL
for(i in tables_list){
  discriminant_variables[[i]]=rownames(data.diablo$loadings[[i]][abs(rowSums(data.diablo$loadings[[i]]))>0,])
}

#Get pair-wise combination of the tables
S<-combn(tables_list,2)
pairwise_correlations<-NULL
for(i in 1:ncol(S)){
  tmp<-expand.grid(discriminant_variables[[S[1,i]]],discriminant_variables[[S[2,i]]])
  if(is.null(pairwise_correlations)){pairwise_correlations=tmp} else {pairwise_correlations=rbind(pairwise_correlations,tmp)}
}

pairwise_correlations$cor=-1
pairwise_correlations$p.value=-1
for(i in 1:nrow(pairwise_correlations)){
  var1=as.character(pairwise_correlations[i,1])
  table1=names(discriminant_variables)[sapply(names(discriminant_variables),function(x){var1 %in% discriminant_variables[[x]]})]
  var2=as.character(pairwise_correlations[i,2])
  table2=names(discriminant_variables)[sapply(names(discriminant_variables),function(x){var2 %in% discriminant_variables[[x]]})]
  #Fixed bug here
  correlation=cor.test(data[[table1]][,var1],data[[table2]][,var2])
  pairwise_correlations[i,"p.value"]=correlation$p.value
  pairwise_correlations[i,"cor"]=correlation$estimate
}
#Now adjust for multiple comparisons
pairwise_correlations$adjp.value<-p.adjust(pairwise_correlations$p.value,method="BH")
pairwise_correlations<-pairwise_correlations[pairwise_correlations$adjp.value<=0.05,]
write.csv(pairwise_correlations,paste("Step.10-DIABLO_significant_correlations_",label,".csv",sep=""))
#/Get all pairwise-correlations and save them in the form of a CSV file##################

save.image(paste("Step.11-DIABLO_DATA_",label,".RData",sep=""))

#for performance of model
perf.diablo <- perf(data.diablo, validation = "Mfold", folds = 10, nrepeat = 10)

# Initialize vectors to store selected features
combined_metagenome_features <- c()
combined_metabolome_features <- c()

# Loop through each repetition
for (i in 1:10) {  # Change 10 to the actual number of repetitions
  # Extract metagenome features for the current repetition and component
  metagenome_features <- perf.diablo$features$stable[[paste0("nrep", i)]]$metagenome
  
  # Extract metabolome features for the current repetition and component
  metabolome_features <- perf.diablo$features$stable[[paste0("nrep", i)]]$metabolome
  
  # Combine features from this repetition into the overall list
  combined_metagenome_features <- c(combined_metagenome_features, metagenome_features)
  combined_metabolome_features <- c(combined_metabolome_features, metabolome_features)
}

library(data.table)
library(dplyr)
metagenome_data <- Map(as.data.frame, combined_metagenome_features) 
metagenome_df <- rbindlist(metagenome_data)
metagenome_stability <- metagenome_df %>%
  group_by(Var1) %>%
  summarise(Count = n(),
            Average = mean(Freq))

write.csv(metagenome_stability,"Step.12-metagenome_stability.csv")

metabolome_data <- Map(as.data.frame, combined_metabolome_features) 
metabolome_df <- rbindlist(metabolome_data)
metabolome_stability <- metabolome_df %>%
  group_by(Var1) %>%
  summarise(Count = n(),
            Average = mean(Freq))

write.csv(metabolome_stability,"Step.12-metabolome_stability.csv")

comp1 <- plotLoadings(data.diablo,comp=1,legend.color=colours[1:length(levels(met_df$Treatment))],contrib=loading_contrib,method=loading_method,title=paste('DIABLO Loadings (comp',1,')',sep=""))
comp2 <- plotLoadings(data.diablo,comp=2,legend.color=colours[1:length(levels(met_df$Treatment))],contrib=loading_contrib,method=loading_method,title=paste('DIABLO Loadings (comp',2,')',sep=""))
comp3 <- plotLoadings(data.diablo,comp=3,legend.color=colours[1:length(levels(met_df$Treatment))],contrib=loading_contrib,method=loading_method,title=paste('DIABLO Loadings (comp',3,')',sep=""))

comp1$metagenome$Comp <- 1
comp2$metagenome$Comp <- 2
comp3$metagenome$Comp <- 3

comp1$metabolome$Comp <- 1
comp2$metabolome$Comp <- 2
comp3$metabolome$Comp <- 3

plotload_metagenome <- rbind(comp1$metagenome, comp2$metagenome, comp3$metagenome)
plotload_metagenome$ID <- rownames(plotload_metagenome)
plotload_metabolome <- rbind(comp1$metabolome, comp2$metabolome, comp3$metabolome)
plotload_metabolome$ID <- rownames(plotload_metabolome)

write.xlsx(plotload_metagenome,"plot_loadings_metagenome.xlsx", row.Names=1)
write.xlsx(plotload_metabolome,"plot_loadings_metabolome.xlsx", row.Names=1)

