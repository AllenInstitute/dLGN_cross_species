setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")


library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(feather)

source("summarySE.R")
load("collected_data_20191220.rda")
cluster_palette <- read.csv("cluster_palette.csv")


#########
## Macaque
#########

load("macaque_allcl.rda")
allcl$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Macaque"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Macaque"])]

mac.meta <- metalist$macaque
mac.meta <- tibble::rownames_to_column(mac.meta, var="sample_id")

mac.anno <- left_join(allcl, mac.meta, by=c("sample_name"="sample_id"))

## gene stats ##
mac.stats <- summarySE(data = mac.anno,
                         measurevar = "Genes.Detected",
                         groupvars = "cluster_id")

mac.gene.stats <- mac.anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>% 
  left_join(mac.stats)

## read stats ##
mac.stats <- summarySE(data = mac.anno,
                         measurevar = "total_reads",
                         groupvars = "cluster_id")

mac.read.stats <- mac.anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>% 
  left_join(mac.stats)






#########
## Human 
#########

load("human_allcl.rda")
allcl$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Human"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Human"])]

hu.meta <- metalist$human

hu.anno <- left_join(allcl, hu.meta, by=c("sample_name"="sample_id"))
colnames(hu.anno)[2:4] <- c("cluster_id","cluster_label", "cluster_color")
hu.anno$Genes.Detected_label <- as.numeric(as.character(hu.anno$Genes.Detected_label))
hu.anno$total_reads_label <- as.numeric(as.character(hu.anno$total_reads_label))

## gene stats ##
hu.stats <- summarySE(data = hu.anno,
                         measurevar = "Genes.Detected_label",
                         groupvars = "cluster_id")

hu.gene.stats <- hu.anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>% 
  left_join(hu.stats)

## read stats ##
hu.stats <- summarySE(data = hu.anno,
                         measurevar = "total_reads_label",
                         groupvars = "cluster_id")

hu.read.stats <- hu.anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>% 
  left_join(hu.stats)






#########
## Mouse 
#########

load("mouse_allcl.rda")
allcl$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Mouse"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Mouse"])]

mo.meta <- metalist$mouse

mo.meta <- tibble::rownames_to_column(mo.meta, var="sample_id")
mo.anno <- left_join(allcl, mo.meta, by=c("sample_name"="sample_id"))


## gene stats ##
mo.stats <- summarySE(data = mo.anno,
                         measurevar = "Genes.Detected",
                         groupvars = "cluster_id")

mo.gene.stats <- mo.anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>% 
  left_join(mo.stats)


## read stats ##
mo.stats <- summarySE(data = mo.anno,
                         measurevar = "total_reads",
                         groupvars = "cluster_id")

mo.read.stats <- mo.anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>% 
  left_join(mo.stats)





#########
## plotting
#########

## ma ##


macp1 <-QC_plot(anno=mac.anno, groupvar="cluster", plot.var="Genes.Detected",stats=mac.gene.stats)
macp1 <- macp1 + scale_y_continuous("Genes Detected (thousands)",
                         limits = c(-3000, 14000), 
                         breaks = seq(0, 14000, 2000),
                         labels = seq(0,14,2)) 


macp2 <- QC_plot(anno=mac.anno, groupvar="cluster", plot.var="total_reads",stats=mac.read.stats)
macp2 <- macp2 + scale_y_continuous("Genes Detected (thousands)",
                         limits = c(-1000000, 7000000), 
                         breaks = seq(0, 7000000, 1000000),
                         labels = seq(0,7,1)) 
macp2[["layers"]][[4]][["mapping"]][["y"]] <- -500000
macp2[["layers"]][[5]][["mapping"]][["y"]] <- -400000



## hu ##

hup1 <- QC_plot(anno=hu.anno, groupvar="cluster", plot.var="Genes.Detected_label",stats=hu.gene.stats)
hup1 <- hup1 + scale_y_continuous("Genes Detected (thousands)",
                         limits = c(-3000, 14000), 
                         breaks = seq(0, 14000, 2000),
                         labels = seq(0,14,2)) 

hup2 <- QC_plot(anno=hu.anno, groupvar="cluster", plot.var="total_reads_label",stats=hu.read.stats)
hup2 <- hup2 + scale_y_continuous("Genes Detected (thousands)",
                         limits = c(-1000000, 7000000), 
                         breaks = seq(0, 7000000, 1000000),
                         labels = seq(0,7,1)) 
hup2[["layers"]][[4]][["mapping"]][["y"]] <- -500000
hup2[["layers"]][[5]][["mapping"]][["y"]] <- -400000

## mo ##
mop1 <- QC_plot(anno=mo.anno, groupvar="cluster", plot.var="Genes.Detected",stats=mo.gene.stats)
mop1 <- mop1 + scale_y_continuous("Genes Detected (thousands)",
                         limits = c(-3000, 14000), 
                         breaks = seq(0, 14000, 2000),
                         labels = seq(0,14,2)) 

mop2 <- QC_plot(anno=mo.anno, groupvar="cluster", plot.var="total_reads",stats=mo.read.stats)
mop2 <- mop2 + scale_y_continuous("Genes Detected (thousands)",
                         limits = c(-1000000, 7000000), 
                         breaks = seq(0, 7000000, 1000000),
                         labels = seq(0,7,1)) 
mop2[["layers"]][[4]][["mapping"]][["y"]] <- -500000
mop2[["layers"]][[5]][["mapping"]][["y"]] <- -400000


layout_matrix <- rbind(c(1,2,3),c(4,5,6))  

ggsave(file="QC_jitters.pdf",gridExtra::marrangeGrob(list(macp1, hup1, mop1, macp2, hup2, mop2),nrow = 2, ncol=3, layout_matrix=layout_matrix), height=8,width=6)



  
ggsave(file.path(dest.d, paste0(prefix,".genes_cluster.alpha.png")),gridExtra::marrangeGrob(g,nrow = 6, ncol=5, layout_matrix=layout_matrix),height=15, width=15)  
  