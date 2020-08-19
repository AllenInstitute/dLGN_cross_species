

########################
## Annotations (+old) for all samples
########################

setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")

library(dplyr)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)
source("doublet.finder.r")
require(Matrix)
source("hicat_modules.r")
require(scrattch.vis)
library(Seurat)
library(grid)
library(gridExtra)
library(scales)
library(reshape2)

load("collected_data_20191220.rda")

st=format(Sys.time(), "%Y%m%d_%H%M_")


#### setup metadata
cluster_palette <- read.csv("cluster_palette.csv")
roi_palette <- read.csv("roi_palette.csv",stringsAsFactors = FALSE)

load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Macaque"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Macaque"])]
allcl_mac$species="Macaque"
allcl_mac$species_color="#939598"
allcl_mac$dissection_roi <- roi_palette$dissection_roi[match(allcl_mac$roi, roi_palette$roi)]
allcl_mac$roi_color <- roi_palette$roi_color[match(allcl_mac$roi, roi_palette$roi)]

load("human_allcl_20200126.rda")
allcl_hum=allcl
allcl_hum$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Human"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Human"])]
allcl_hum$species="Human"
allcl_hum$species_color="#000000"
allcl_hum$dissection_roi <- roi_palette$dissection_roi[match(allcl_hum$roi, roi_palette$roi)]
allcl_hum$roi_color <- roi_palette$roi_color[match(allcl_hum$roi, roi_palette$roi)]

load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Mouse"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Mouse"])]
allcl_mus$species="Mouse"
allcl_mus$species_color="#FF73D0"
allcl_mus$dissection_roi <- roi_palette$dissection_roi[match(allcl_mus$roi, roi_palette$roi)]
allcl_mus$roi_color <- roi_palette$roi_color[match(allcl_mus$roi, roi_palette$roi)]


#### read in historic cluster calls

hu.old <- read.csv("human_LGN_2018-06-14_samples-columns.csv")
ma.old <- read.csv("macaque_LGN_2018-06-14_samples-columns.csv")
load("20191206_LGN_Mmul_10_star_samp.dat.comb.Rdata")
mo.old <- read.csv("mouse_LGd_2018-06-14_samples-columns.csv")

hu.old$sample_name <- hu.old$seq_name
anno.hum <- left_join(allcl_hum, hu.old)

samp.dat$sample_name <- rownames(samp.dat)
anno.mac <- left_join(allcl_mac, samp.dat)
ma.old$sample_name <- ma.old$seq_name
anno.mac <- left_join(anno.mac, select(ma.old, sample_name, class, cluster))

mo.old$sample_name <- mo.old$seq_name
anno.mus <- left_join(allcl_mus, mo.old)

save(anno.hum, anno.mac, anno.mus, file= "anno_species.rda")

colnames(anno.mac)[colnames(anno.mac)=="cell_prep_type"] <- "sample_type"
anno.mac$organism <- "macaca nemestrina"

colnames(anno.hum)[which(names(anno.hum)=="genes_detected_cpm_criterion")] <- "Genes.Detected.CPM"
colnames(anno.mus)[which(names(anno.mus)=="genes_detected_cpm_criterion")] <- "Genes.Detected.CPM"

intersect(colnames(anno.mus), colnames(anno.mac))

sel.cols <- c("sample_name",                  
			"cluster_id",                   
			"cluster_label" ,               
			"cluster_color" ,                                                    
			"species" ,                     
			"species_color"   , 
			"organism", 
			"donor" ,           
			"dissection_roi"  ,             
			"roi"   , 
			"roi_color" ,
			"class",
			"cluster",
			"total_reads",
			"Genes.Detected.CPM"
			      )

anno.mac <- anno.mac[,sel.cols]
anno.hum <- anno.hum[,sel.cols]
anno.mus <- anno.mus[,sel.cols]

anno.mus$donor <- as.character(anno.mus$donor)

anno.df <- anno.mac %>% full_join(anno.hum) %>% full_join(anno.mus)




########################
## Donor signatures 
########################


# Macaque (snRNA-seq)
load("20191206_LGN_Mmul_10_star_cpm.Rdata")
load("anno_species.rda")

rownames(anno.mac) <- anno.mac$sample_name

samp.rh.nuc.subset <- anno.mac[anno.mac$cluster_label %in% c("MP1", "MP2"),]

keep.rh.nuc.samp <- anno.mac$sample_name[anno.mac$cluster_label %in% c("MP1", "MP2")]
expr.rh.nuc.subset <- log2(cpmR[, keep.rh.nuc.samp] + 1)


#### Calculate macaque magno vs. parvo DEX genes (snRNA-seq)
dex.genes.rh.nuc <- list()


# MC/PC DEX
design <- model.matrix(~ 0 + samp.rh.nuc.subset$cluster_label)
colnames(design) <- c("MP1","MP2")

fit <- lmFit(expr.rh.nuc.subset, design)
for (contrast1 in c("MP1-MP2")) {
  cont.matrix <- makeContrasts(contrasts=contrast1, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
  top1 <- topTable(fit2, number = Inf, p.value = 1, 
                   adjust = "BH", sort.by = "none")
  top1$gene <- row.names(top1)
  dex.genes.rh.nuc[[contrast1]] <- top1
}


# MC/PC DEX (by donor)
for (donor1 in unique(samp.rh.nuc.subset$donor)) {
  keep.rh.nuc.donor <- which(samp.rh.nuc.subset$donor == donor1)
  expr.rh.nuc.donor <- expr.rh.nuc.subset[, keep.rh.nuc.donor]
  samp.rh.nuc.donor <- droplevels(samp.rh.nuc.subset[keep.rh.nuc.donor, ])
  
  # LGN DEX
  design <- model.matrix(~ 0 + samp.rh.nuc.donor$cluster_label)
  colnames(design) <- c("MP1","MP2")
  fit <- lmFit(expr.rh.nuc.donor, design)
  
  cont.matrix <- makeContrasts(contrasts="MP1-MP2", levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
  top1 <- topTable(fit2, number = Inf, p.value = 1, 
                   adjust = "BH", sort.by = "none")
  top1$gene <- row.names(top1)
  colnames(top1) <- paste0(colnames(top1), "_", donor1)
  dex.genes.rh.nuc[["MP1-MP2"]] <- cbind(dex.genes.rh.nuc[["MP1-MP2"]], top1)
}

dex.genes.mac <-  dex.genes.rh.nuc[["MP1-MP2"]]
write.csv(dex.genes.mac, file="deg_mac_mp1_mp2_donor.csv")


# Human (snRNA-seq)

# Load human data
expr.h <- as.data.frame(read_feather("human_data.feather"))
row.names(expr.h) <- expr.h$sample_id
expr.h <- t(expr.h[, ! grepl("sample_id", colnames(expr.h))])
load("anno_species.rda")

rownames(anno.hum) <- anno.hum$sample_name

samp.hum.subset <- anno.hum[anno.hum$cluster_label == "MP1",]
expr.h.subset <- logCPM(expr.h[,samp.hum.subset$sample_name])


#### Calculate human magno vs. parvo DEX genes (snRNA-seq)
dex.genes.hum <- list()


# MC/PC DEX
design <- model.matrix(~ 0 + samp.hum.subset$dissection_roi)
colnames(design) <- c("KC","MC","PC")

fit <- lmFit(expr.h.subset, design)
for (contrast1 in c("MC-PC")) {
  cont.matrix <- makeContrasts(contrasts=contrast1, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
  top1 <- topTable(fit2, number = Inf, p.value = 1, 
                   adjust = "BH", sort.by = "none")
  top1$gene <- row.names(top1)
  dex.genes.hum[[contrast1]] <- top1
}


# MC/PC DEX (by donor)
for (donor1 in unique(samp.hum.subset$donor)) {
  keep.hum.donor <- which(samp.hum.subset$donor == donor1)
  expr.hum.donor <- expr.h.subset[, keep.hum.donor]
  samp.hum.donor <- droplevels(samp.hum.subset[keep.hum.donor, ])
  
  # LGN DEX
  design <- model.matrix(~ 0 + samp.hum.donor$dissection_roi)
  colnames(design) <- c("KC","MC","PC")
  fit <- lmFit(expr.hum.donor, design)
  
  cont.matrix <- makeContrasts(contrasts="MC-PC", levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
  top1 <- topTable(fit2, number = Inf, p.value = 1, 
                   adjust = "BH", sort.by = "none")
  top1$gene <- row.names(top1)
  colnames(top1) <- paste0(colnames(top1), "_donor.", donor1)
  dex.genes.hum[["MC-PC"]] <- cbind(dex.genes.hum[["MC-PC"]], top1)
}

dex.genes.hum <-  dex.genes.hum[["MC-PC"]]
write.csv(dex.genes.hum, file="deg_hum_mp1_mcpc_donor.csv")
