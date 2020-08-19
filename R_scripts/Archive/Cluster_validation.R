
########################
## Setup
########################

setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)
source("doublet.finder.r")
require(Matrix)
source("hicat_modules.r")
require(scrattch.vis)
require(Seurat)
library("RColorBrewer")
library(gplots)

load("collected_data_20191220.rda")

st=format(Sys.time(), "%Y%m%d_%H%M_")

########################
## Assess clustering robustness
########################


getConfusionMatrix <- function(orgCluster,
                               foundCluster,
                               proportions = TRUE) {
    orgCluster <- as.character(orgCluster)
    foundCluster <- as.character(foundCluster)
    #lev <- sort(unique(c(orgCluster, foundCluster)))
    #orgCluster <- factor(orgCluster, levels = lev)
    #foundCluster <- factor(orgCluster, levels = lev)
    confusion <- table(foundCluster, orgCluster)
    if (proportions) {
        cs <- colSums(confusion)
        for (i in 1:dim(confusion)[1])
            confusion[i, ] <- confusion[i, ] / pmax(cs, 1e-08)
    }
    confusion
}

#=================
##=== Macaque ===#
#=================

load("macaque_clusters_subclusters.rda")
m.int <- FindVariableFeatures(object = m.int, selection.method = "vst", nfeatures = 2000)
mac.markers <- VariableFeatures(m.int)

load("macaque_allcl.rda")
cl.consensus <- setNames(allcl$cluster_id, allcl$sample_name)
cl.df.consensus <- allcl %>% group_by(cluster_id, cluster_label) %>% summarise(size=n())

norm.dat<-as.matrix(datlist[["macaque"]][,allcl$sample_name])

cell.cl.map.df = get_core_transition(norm.dat, cl.consensus, mac.markers, n.bin=5, n.iter=100, mc.cores=1)

cell.cl.map.df$org.label <- cl.df.consensus$cluster_label[match(cell.cl.map.df$org.cl,cl.df.consensus$cluster_id)]
cell.cl.map.df$best.label <- cl.df.consensus$cluster_label[match(cell.cl.map.df$best.cl,cl.df.consensus$cluster_id)]

membConfusionProp  <- getConfusionMatrix(cell.cl.map.df$org.label,cell.cl.map.df$best.label,TRUE)
clOrd <- cl.df.consensus$cluster_label # Cluster order

library(paletteer)
colfunc <- paletteer_c("pals::ocean.tempo", n = 6)
library(cartography)
colfunc<-colorRampPalette(carto.pal(pal1 = "harmo.pal" ,n1 = 5))(n=100)

#colfunc <- colorRampPalette(c("#F2EEE7FF","#5CA687FF", "#147C77FF", "#1A4D5FFF", "#151D44FF"))(n=100)
colfunc <- colorRampPalette(c("#FDF1E2","#A77697" ,"#71408C","#1E315B"))(n=100)
pdf("macaque_confusion.pdf", height=10, width=10, useDingbats=FALSE)
heatmap.2(pmin(membConfusionProp,0.25)[clOrd,clOrd],col=colfunc ,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none", main="Confusion Matrix")
dev.off()

#===============
##=== Human ===#
#===============

load("human_clusters_subclusters_20200126.rda")

m.int <- FindVariableFeatures(object = m.int, selection.method = "vst", nfeatures = 2000)
hum.markers <- VariableFeatures(m.int)

load("human_allcl_20200126.rda")
cl.consensus <- setNames(allcl$cluster_id, allcl$sample_name)
cl.df.consensus <- allcl %>% group_by(cluster_id, cluster_label) %>% summarise(size=n())

norm.dat<-as.matrix(datlist[["human"]][,allcl$sample_name])

cell.cl.map.df = get_core_transition(norm.dat, cl.consensus, rownames(norm.dat), n.bin=5, n.iter=100, mc.cores=1)

cell.cl.map.df$org.label <- cl.df.consensus$cluster_label[match(cell.cl.map.df$org.cl,cl.df.consensus$cluster_id)]
cell.cl.map.df$best.label <- cl.df.consensus$cluster_label[match(cell.cl.map.df$best.cl,cl.df.consensus$cluster_id)]

membConfusionProp  <- getConfusionMatrix(cell.cl.map.df$org.label,cell.cl.map.df$best.label,TRUE)
clOrd <- cl.df.consensus$cluster_label # Cluster order

library(paletteer)
colfunc <- paletteer_c("pals::ocean.tempo", n = 6)
library(cartography)
colfunc<-colorRampPalette(carto.pal(pal1 = "harmo.pal" ,n1 = 5))(n=100)

#colfunc <- colorRampPalette(c("#F2EEE7FF","#5CA687FF", "#147C77FF", "#1A4D5FFF", "#151D44FF"))(n=100)
colfunc <- colorRampPalette(c("#FDF1E2","#A77697" ,"#71408C","#1E315B"))(n=100)
pdf("human_confusion.pdf", height=10, width=10, useDingbats=FALSE)
heatmap.2(pmin(membConfusionProp,0.25)[clOrd,clOrd],col=colfunc ,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none", main="Confusion Matrix")
dev.off()

## split per donor
donors <- as.character(unique(allcl$donor))
for(donor in donors) {

	subcl <- allcl[allcl$donor == donor,]
	cl <- setNames(subcl$cluster_id, subcl$sample_name)
	cldf <- subcl %>% group_by(cluster_id, cluster_label) %>% summarise(size=n())
	
	norm.dat<-as.matrix(datlist[["human"]][,subcl$sample_name])

	cell.cl.map.df = get_core_transition(norm.dat, cl, hum.markers, n.bin=5, n.iter=100, mc.cores=1)

	cell.cl.map.df$org.label <- cldf$cluster_label[match(cell.cl.map.df$org.cl,cldf$cluster_id)]
	cell.cl.map.df$best.label <- cldf$cluster_label[match(cell.cl.map.df$best.cl,cldf$cluster_id)]

	membConfusionProp  <- getConfusionMatrix(cell.cl.map.df$org.label,cell.cl.map.df$best.label,TRUE)
	clOrd <- cldf$cluster_label # Cluster order

	colfunc <- colorRampPalette(c("#FDF1E2","#A77697" ,"#71408C","#1E315B"))(n=100)
	pdf(paste0("human",donor,"confusion.pdf"), height=10, width=10, useDingbats=FALSE)
	heatmap.2(pmin(membConfusionProp,0.25)[clOrd,clOrd],col=colfunc ,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none", main="Confusion Matrix")
dev.off()

}

#===============
##=== Mouse ===#
#===============

load("mouse_clusters_subclusters.rda")
m.int <- FindVariableFeatures(object = m.int, selection.method = "vst", nfeatures = 2000)
mo.markers <- VariableFeatures(m.int)

load("mouse_allcl.rda")


cl.consensus <- setNames(allcl$cluster_id, allcl$sample_name)
cl.df.consensus <- allcl %>% group_by(cluster_id, cluster_label) %>% summarise(size=n())

norm.dat<-as.matrix(datlist[["mouse"]][,allcl$sample_name])

cell.cl.map.df = get_core_transition(norm.dat, cl.consensus, mo.markers, n.bin=5, n.iter=100, mc.cores=1)

cell.cl.map.df$org.label <- cl.df.consensus$cluster_label[match(cell.cl.map.df$org.cl,cl.df.consensus$cluster_id)]
cell.cl.map.df$best.label <- cl.df.consensus$cluster_label[match(cell.cl.map.df$best.cl,cl.df.consensus$cluster_id)]

membConfusionProp  <- getConfusionMatrix(cell.cl.map.df$org.label,cell.cl.map.df$best.label,TRUE)
clOrd <- cl.df.consensus$cluster_label # Cluster order

library(paletteer)
colfunc <- paletteer_c("pals::ocean.tempo", n = 6)
library(cartography)
colfunc<-colorRampPalette(carto.pal(pal1 = "harmo.pal" ,n1 = 5))(n=100)

#colfunc <- colorRampPalette(c("#F2EEE7FF","#5CA687FF", "#147C77FF", "#1A4D5FFF", "#151D44FF"))(n=100)
colfunc <- colorRampPalette(c("#FDF1E2","#A77697" ,"#71408C","#1E315B"))(n=100)
pdf("mouse_confusion.pdf", height=10, width=10, useDingbats=FALSE)
heatmap.2(pmin(membConfusionProp,0.25)[clOrd,clOrd],col=colfunc ,Rowv=FALSE,Colv=FALSE,trace="none",dendrogram="none", main="Confusion Matrix")
dev.off()



