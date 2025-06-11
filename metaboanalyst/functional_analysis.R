library(MetaboAnalystR)
library(readxl)
library(openxlsx)

mSet<-InitDataObjects("mass_all", "mummichog", FALSE)
mSet<-SetPeakFormat(mSet, "rmp")
mSet<-UpdateInstrumentParameters(mSet, 5.0, "negative", "no", 0.02);
mSet<-Read.PeakListData(mSet, "../neg_summary.txt");
mSet<-SanityCheckMummichogData(mSet)
mSet<-SetPeakEnrichMethod(mSet, "integ", "v2")
mSet<-SetMummichogPval(mSet, 0.05)
mSet<-PerformPSEA(mSet, "gga_kegg", "current", 3 , 100)
mSet<-PlotPSEAIntegPaths(mSet, "plot2", "png", 600, width=NA, labels=10,scale=FALSE)
mSet<-SaveTransformedData(mSet)