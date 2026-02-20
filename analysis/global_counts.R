library(readxl)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(ggpubr)

# Importing metabolite classification file
df0 <- read_excel("global_metabolite_list.xlsx")


# Selecting columns to include
classification <- df0[,c("ID","ID_name","Name","Super.Class","Sub.Class","Class","Final.Class",
                         "Super.Pathway","Sub.Pathway","Pathway")]

# Renaming super classes
classification[classification$Super.Class=="NA","Super.Class"] <- "Unclassified"
classification[classification$Super.Class=="FA  Fatty acyls","Super.Class"] <- "Fatty acyls"
classification[classification$Super.Class=="ST  Sterol lipids","Super.Class"] <- "Sterol lipids"
classification[classification$Super.Class=="PK  Polyketides","Super.Class"] <- "Polyketides"
classification[classification$Super.Class=="GL  Glycerolipids","Super.Class"] <- "Glycerolipids"
classification[classification$Super.Class=="GP  Glycerophospholipids","Super.Class"] <- "Glycerophospholipids"
classification[classification$Super.Class=="PR  Prenol lipids","Super.Class"] <- "Prenol lipids"
classification[classification$Super.Class=="SP  Sphingolipids","Super.Class"] <- "Sphingolipids"
classification[classification$Super.Class=="Others","Super.Class"] <- "Unclassified"

# NAs as unclassified
classification[classification$Super.Pathway=="NA","Super.Pathway"] <- "Unclassified"

# Grouping by superclass
Superclass_tally <- classification %>% group_by(Super.Class) %>% tally() %>%
  arrange(n)

# Getting top superclasses
Superclass_tally$Superclass <- ifelse(Superclass_tally$n >= 13, Superclass_tally$Super.Class, "Others")
# Getting counts for 'others;
Others_sum <- sum(Superclass_tally[which(Superclass_tally$Superclass=="Others"),2])

# Plotting superclass distribution
Superclass_df <- Superclass_tally[Superclass_tally$n>=13,]
Superclass_df <- rbind(Superclass_df, list("Others", Others_sum, "Others"))
Superclass_df$per <- round((Superclass_df$n)/7554*100,2)
Superclass_df$n2 <- paste(Superclass_df$n, " (",Superclass_df$per,"%)", sep="")

p <- ggplot(Superclass_df, aes(x = n, y = reorder(Superclass, n))) +
  geom_bar(stat="identity", color="#292929",fill='steelblue') +
  geom_text(aes(label = n2), nudge_x = 250, size = 5) +
  theme_bw() +
  labs(x="Number of metabolites",
       y="Superclass") +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title.x = element_text(size=18, face = "bold", vjust = -5),
        axis.title.y = element_text(size=18, face = "bold", vjust = +5),
        title = element_text(face="bold", size = 20)) +
  theme(plot.margin = margin(1.2,1.2,1.2,1.2, "cm")); p



# Plotting superpathway distribution
Superpathway_tally <- classification %>% group_by(Super.Pathway) %>% tally() %>%
  arrange(n)
Superpathway_tally$per <- round((Superpathway_tally$n)/7554*100,2)
Superpathway_tally$n2 <- paste(Superpathway_tally$n, " (",Superpathway_tally$per,"%)", sep="")


p2 <- ggplot(Superpathway_tally, aes(x = n, y = reorder(Super.Pathway, n))) +
  geom_bar(stat="identity", color="#292929",fill='steelblue') +
  geom_text(aes(label = n2), nudge_x = 450, size = 5) +
  theme_bw() +
  labs(x="Number of metabolites",
       y="Superpathway") +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title.x = element_text(size=18, face = "bold", vjust = -5),
        axis.title.y = element_text(size=18, face = "bold", vjust = +5),
        title = element_text(face="bold", size = 20)) +
  theme(plot.margin = margin(1.2,1.2,1.2,1.2, "cm")); p2


# Merging plots
plots <- ggarrange(p, p2, ncol = 2, nrow = 1,
                   font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                   labels = c("(a)","(b)"))

ggsave("Globaloverview.png",height=10,width=30, units="in", dpi = 600)
print(plots)
dev.off()

