library(dplyr)
library(readxl)
library(openxlsx)

df0 <- read.csv("metabolite_abundance.csv", row.names=1,header=TRUE)

df0[df0$Treatment=="No","Treatment"] <- "T1"
df0[df0$Treatment=="Yes","Treatment"] <- "T2"

# Kolmogorov-Smirnov Test for each column
ks_test_results <- apply(df0, 2, function(column) ks.test(column, "pnorm", mean(column), sd(column))$p.value)

# Print p-values for each column
ks_test_results_df <- data.frame(ks_test_results)

#For normally distributed
t_test_mets <- ks_test_results_df %>%
  filter(ks_test_results>0.05)
t_test_mets <- rownames(t_test_mets)
t_test_df <- df0[,t_test_mets]
t_test_df$Treatment <- c(rep("T1",10),rep("T2",10))
t_test_df <- data.frame(t_test_df[ncol(t_test_df)],t_test_df[1:ncol(t_test_df)-1])

# Separate groups
group1_t_test_df <- t_test_df[t_test_df$Treatment == "T1", ]
group2_t_test_df <- t_test_df[t_test_df$Treatment == "T2", ]
result_t_test_df <- apply(t_test_df[, -1], 2, function(x) {
  t.test(x ~ t_test_df$Treatment)
})

p_values <- sapply(result_t_test_df, function(x) x$p.value)
t_scores <- sapply(result_t_test_df, function(x) x$statistic)

summary_t_test_df <- data.frame(
  Metabolite = colnames(t_test_df[, -1]),
  p_value = p_values,
  t_score = t_scores
)

#For non-normal, wilcox instead of t-test

w_test_mets <- ks_test_results_df %>%
  filter(ks_test_results<0.05)
w_test_mets <- rownames(w_test_mets)
w_test_df <- df0[,w_test_mets]
w_test_df$Treatment <- c(rep("T1",10),rep("T2",10))
w_test_df <- data.frame(w_test_df[ncol(w_test_df)],w_test_df[1:ncol(w_test_df)-1])
# Separate groups
group1_w_test_df <- w_test_df[w_test_df$Treatment == "T1", ]
group2_w_test_df <- w_test_df[w_test_df$Treatment == "T2", ]

result_w_test_df <- apply(w_test_df[, -1], 2, function(x) {
  wilcox.test(x ~ w_test_df$Treatment)
})
p_values2 <- sapply(result_w_test_df, function(x) x$p.value)
t_scores2 <- sapply(result_w_test_df, function(x) x$statistic)

# Create a summary table
summary_w_test_df <- data.frame(
  Metabolite = colnames(w_test_df[, -1]),
  p_value = p_values2,
  t_score = t_scores2
)
all_summary <- rbind(summary_t_test_df, summary_w_test_df)

# Adding metabolite information
df1 <- as.data.frame(read_excel("global_metabolite_list.xlsx"))
all_summary$mz <- df1$`m/z`[match(all_summary$Metabolite,df1$ID)]
all_summary$rtime <- df1$RT[match(all_summary$Metabolite,df1$ID)]
all_summary$Name <- df1$Name[match(all_summary$Metabolite,df1$ID)]
all_summary$ID <- df1$ID[match(all_summary$Metabolite,df1$ID)]
all_summary$ION<- df1$Ion.Mode[match(all_summary$Metabolite,df1$ID)]

#for table files

neg_summary <- all_summary[all_summary$ION=="NEG",]
neg_summary <- data.frame(neg_summary[c(4:5,2:3)])
pos_summary <- all_summary[all_summary$ION=="POS",]
pos_summary <- data.frame(pos_summary[c(4:5,2:3)])
# write.table(neg_summary,"neg_summary.txt",row.names = FALSE)
# write.table(pos_summary,"pos_summary.txt",row.names = FALSE)

#when t scores > 0, T1, <0 T2
