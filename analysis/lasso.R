
library(glmnet)

#===========Parameters=============#
#Load the whole meta table, assure you have the treatment column T1 = No, T2 = Yes
meta_table<-read.csv("metabolite_abundance.csv",row.names=1,header=TRUE,check.names=FALSE)
#List the metabolites in that table
mets <- colnames(meta_table)[1:(ncol(meta_table)-1)]
metabolites_list<-mets
dependent_variable<-"Treatment"

#Get rid of any samples for which dependent variable is missing
meta_table[,dependent_variable][meta_table[,dependent_variable]==""]<-NA
meta_table<-meta_table[complete.cases(meta_table[,dependent_variable,drop=FALSE]),,drop=FALSE]
#Make sure the presence/absence column is coded as Yes/No - No (T1). Yes (T2)
meta_table[,dependent_variable]<-as.factor(meta_table[,dependent_variable])
meta_table[,dependent_variable]<-relevel(meta_table[,dependent_variable],"No" )
#===========Hypothesis space========

#In the first instance, shorten the metabolites name, which we can recover later
a<-meta_table[,metabolites_list]
#Make a mapping table
mapping<-data.frame(row.names=colnames(a),paste0("M",seq(1,ncol(a))))
#Use the new mappings
colnames(a)<-mapping[colnames(a),]

#Reference: https://www.statology.org/lasso-regression-in-r/

#use below for no normalisation #BGI already normalised this
data_table <- as.matrix(meta_table)
cv_model <- cv.glmnet(data_table, meta_table[,dependent_variable], alpha = 1,family="binomial")

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min

#find coefficients of best model
best_model <- glmnet(data_table, meta_table[,dependent_variable], alpha = 1, lambda = best_lambda,family="binomial")
coef(best_model)
write.csv(as.matrix(coef(best_model)),"lasso.csv")


