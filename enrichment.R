#bubble plots

library(openxlsx)

#load from metaboanalyst output files
posdf <- read.csv("metaboanalyst/Positive/mummicho_pathway_enrichment_integ.csv")
negdf <- read.csv("metaboanalyst/Negative/mummicho_pathway_enrichment_integ.csv")

#load from metaboanalyst output files
posdf2 <- read.csv("metaboanalyst/Positive/mummichog_pathway_enrichment_integ.csv")
negdf2 <- read.csv("metaboanalyst/Negative/mummichog_pathway_enrichment_integ.csv")

#Merge them, get it from df2 add it to df2
posdf$cpd.hits <- posdf2$cpd.hits[match(posdf$X,posdf2$X)]
negdf$cpd.hits <- negdf2$cpd.hits[match(negdf$X,negdf2$X)]

#Importing files with t values stat
pos_stat <- read.csv("metaboanalyst/Positive/data_original.csv")
posdf3 <- read.csv("metaboanalyst/Positive/mummichog_matched_compound_all.csv")
colnames(posdf3) <- c("Mass","Compound","MatchedForm","MassDiff")
neg_stat <- read.csv("metaboanalyst/Negative/data_original.csv")
negdf3 <- read.csv("metaboanalyst/Negative/mummichog_matched_compound_all.csv")
colnames(negdf3) <- c("Mass","Compound","MatchedForm","MassDiff")

### Positive ion mode - get which compounds are elevated in T1 or T2

#For loop to get compound hits per pathway in positive mode, only including top 4
df0 <- data.frame()
for (x in 1:4) {
  Compound <- unlist(strsplit(posdf$cpd.hits[x], ";"))
  num <- as.numeric(length(unlist(strsplit(posdf$cpd.hits[x], ";"))))
  Pathway <- rep(posdf[x,"X"],num)
  Mode <- rep("Pos",num)
  tmp <- data.frame(Compound,Pathway,Mode)
  df0 <- rbind(df0,tmp)
}

#Match compound id, by getting query mass, then which matches to another file with tval, pval
df0_merged <- merge(df0, posdf3, by = "Compound")
df0_merged$Stat <- pos_stat$t.score[match(df0_merged$Mass,pos_stat$m.z)]
df0_merged$Pval <- pos_stat$p.value[match(df0_merged$Mass,pos_stat$m.z)]
df0_merged$Treatment <- ifelse(df0_merged$Stat>0,"T1","T2") 

### Negative ion mode - get which compounds are elevated in T1 or T2

#negative mode, including top 8 so 1:8
df1 <- data.frame()
for (x in 1:8) {
  Compound <- unlist(strsplit(negdf$cpd.hits[x], ";"))
  num <- as.numeric(length(unlist(strsplit(negdf$cpd.hits[x], ";"))))
  Pathway <- rep(negdf[x,"X"],num)
  Mode <- rep("Neg",num)
  tmp <- data.frame(Compound,Pathway,Mode)
  df1 <- rbind(df1,tmp)
}
df1_merged <- merge(df1, negdf3, by = "Compound")
df1_merged$Stat <- neg_stat$t.score[match(df1_merged$Mass,neg_stat$m.z)]
df1_merged$Pval <- neg_stat$p.value[match(df1_merged$Mass,neg_stat$m.z)]
df1_merged$Treatment <- ifelse(df1_merged$Stat>0,"T1","T2")


#Bind both pos and neg compound hits
df_all <- rbind(df0_merged,df1_merged)

# write.xlsx(df_all,"EnrichmentKEGGcompounds_stat.xlsx") # to export if needed

#Getting number of those significant elevated per treatment and per pathway
# Below is getting all per compound, a compound can have multiple matches
T1 <- c()
T2 <- c()
dfm <- rbind(posdf,negdf)
Pathways <- c(posdf$X[1:4],negdf$X[1:8])
dfm <- dfm[dfm$Combined_Pvals<0.251,] # using 0.251 to only get top 10 pathways

# To get how many compounds per pathway were elevated according to treatment
for (x in Pathways) {
  a <- nrow(df_all[df_all$Pathway==x & df_all$Pval < 0.05 & df_all$Treatment =="T1",])
  T1 <- c(T1,a)
  b <- nrow(df_all[df_all$Pathway==x & df_all$Pval < 0.05 & df_all$Treatment =="T2",])
  T2 <- c(T2,b)
}
dfm$Path <- Pathways
dfm$T1 <- T1
dfm$T2 <- T2
# write.xlsx(dfm,"Enrichment_table_duplicatedcompounds.xlsx")
