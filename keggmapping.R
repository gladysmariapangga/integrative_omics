library(pathview)
#bubble plots

library(openxlsx)

#Load the following dataframes from MetaboAnalyst results
#These show the top 5 pathways
posdf <- read.csv("metaboanalyst/Positive/mummicho_pathway_enrichment_integ.csv")
posdf$Top <- c(rep("Top",3),rep("No",nrow(posdf)-3))
posdf <- posdf[posdf$Top=="Top",]
negdf <- read.csv("metaboanalyst/Negative/mummicho_pathway_enrichment_integ.csv")
negdf$Top <- c(rep("Top",8),rep("No",nrow(negdf)-8))
negdf <- negdf[negdf$Top=="Top",]

#Import these datafiles, because they should the compound hits
posdf2 <- read.csv("metaboanalyst/Positive/mummichog_pathway_enrichment_integ.csv")
negdf2 <- read.csv("metaboanalyst/Negative/mummichog_pathway_enrichment_integ.csv")

#Merge them, get it from df2 add it to df2
posdf$cpd.hits <- posdf2$cpd.hits[match(posdf$X,posdf2$X)]
negdf$cpd.hits <- negdf2$cpd.hits[match(negdf$X,negdf2$X)]


#Getting list of compounds in enriched pathways
pos1<- unlist(strsplit(posdf$cpd.hits[1], ";")) #"Phenylalanine metabolism"  
pos2<- unlist(strsplit(posdf$cpd.hits[2], ";")) #"Tryptophan metabolism"   
pos3<- unlist(strsplit(posdf$cpd.hits[3], ";")) # "Phenylalanine, tyrosine and tryptophan biosynthesis" 

neg1<- unlist(strsplit(negdf$cpd.hits[1], ";")) # "Glutathione metabolism"  
neg2<- unlist(strsplit(negdf$cpd.hits[2], ";")) # "Glycerophospholipid metabolism" 
neg3<- unlist(strsplit(negdf$cpd.hits[3], ";")) # "beta-Alanine metabolism"
neg5<- unlist(strsplit(negdf$cpd.hits[5], ";")) # "N-Glycan biosynthesis"
neg8<- unlist(strsplit(negdf$cpd.hits[8], ";")) # "Thiamine metabolism" 

#POSITIVE
#phenylalalnine 00360
pathview(gene.data = NULL, cpd.data = pos1, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00360", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga")

#tryptophan
pathview(gene.data = NULL, cpd.data = pos2, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00380", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga") 

#Phenylalanine, tyrosine and tryptophan biosynthesis map00400
pathview(gene.data = NULL, cpd.data = pos3, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00400", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga")



#NEGATIVE
#glutathione map00480
pathview(gene.data = NULL, cpd.data = neg1, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00480", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga")

#Glycerophospholipid metabolism map00564
pathview(gene.data = NULL, cpd.data = neg2, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00564", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga")

#beta-Alanine metabolism map00410	
pathview(gene.data = NULL, cpd.data = neg3, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00410", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga")

#N-Glycan biosynthesis map00510
pathview(gene.data = NULL, cpd.data = neg5, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00510", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga")


#Thiamine metabolism 00730	
pathview(gene.data = NULL, cpd.data = neg8, cpd.idtype = "kegg", 
         map.cpdname = FALSE,
         pathway.id = "gga00730", 
         expand.node = TRUE,
         kegg.native = FALSE,        # Use KEGG native pathway graph
         out.suffix = "compound_map",
         high = list(cpd = "lightblue"),
         species = "gga")



