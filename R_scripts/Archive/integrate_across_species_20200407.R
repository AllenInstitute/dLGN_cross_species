setwd("C:/Users/menonv/Dropbox (HHMI)/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")
####integrate human and macaque data###
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(scrattch.vis)
library(feather)
require(Matrix)

load("collected_data_20191220.rda")

spec="macaque"
orthos=read.csv("ortholog_table_20191122.csv",as.is=T)
nonlocgenes=orthos$rhesus_symbol[intersect(grep("^LOC",orthos$human_symbol,invert=T),grep("^LOC",orthos$rhesus_symbol))]
keepgen=setdiff(rownames(datlist[[spec]]),setdiff(grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|^TRNA",rownames(datlist[[spec]]),val=T),nonlocgenes))
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
spec="human"
keepgen=grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|[0-9]P$|-AS|-PS|DUX4L",rownames(datlist[[spec]]),val=T,invert=T)
keepdons=which(metalist$human$external_donor_name_label %in% c("H200.1023","H200.1025","H200.1030"))
datlist[[spec]]=datlist[[spec]][keepgen,keepdons]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
spec="mouse"
keepgen=grep("^LOC|^mt-|^Rpl|^Rps",rownames(datlist[[spec]]),val=T,invert=T)
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6


cluster_palette <- read.csv("cluster_palette.csv")

load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Macaque"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Macaque"])]
allcl_mac$species="Macaque"
load("human_allcl_20200126.rda")
allcl_hum=allcl
allcl_hum$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Human"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Human"])]
allcl_hum$species="Human"
load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Mouse"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Mouse"])]
allcl_mus$species="Mouse"

####integrate all neuronal cell types: note that the macaque donors need to be integrated separately###
classkeep=c("GA","GL")
keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Macaque"])
keepcells_mac1=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q18"))]

keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Human"])
keepcells_hum=allcl_hum$sample_name[which(allcl_hum$cluster_label %in% keepclasses)]

keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Mouse"])
keepcells_mus=allcl_mus$sample_name[which(allcl_mus$cluster_label %in% keepclasses)]

datlist2=datlist
for (ii in 1:length(datlist2)) {
  rownames(datlist2[[ii]])=tolower(rownames(datlist2[[ii]]))
}
obj.list = list()
obj.list[["macaque1"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac1]),meta.data=allcl_mac[keepcells_mac1,])
obj.list[["macaque2"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac2]),meta.data=allcl_mac[keepcells_mac2,])
obj.list[["human"]]=CreateSeuratObject(round(datlist2[["human"]][,keepcells_hum]),meta.data=allcl_hum[keepcells_hum,])
obj.list[["mouse"]]=CreateSeuratObject(round(datlist2[["mouse"]][,keepcells_mus]),meta.data=allcl_mus[keepcells_mus,])
for (ii in 1:length(obj.list)) {
  obj.list[[ii]] <- SCTransform(obj.list[[ii]], verbose = FALSE)
}
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

require(ggplot2)
require(cowplot)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
uniqcl=unique(obj.integrated@meta.data$cluster_label)
uniqcl=uniqcl[order(uniqcl)]
colvec=obj.integrated@meta.data$cluster_color[match(uniqcl,obj.integrated@meta.data$cluster_label)]
pdf("human_macaque_mouse_allneurons_integrated_umap_dimensions_20200407.pdf",width=15)
for (dimval in c(5,10,15,20,25,30)) {
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
  p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T,split.by="species",cols=colvec)
  print(p1)
}
dev.off()


#########All Glutamatergic cells###
classkeep=c("GL")
keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Macaque"])
keepcells_mac1=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q18"))]

keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Human"])
keepcells_hum=allcl_hum$sample_name[which(allcl_hum$cluster_label %in% keepclasses)]

keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Mouse"])
keepcells_mus=allcl_mus$sample_name[which(allcl_mus$cluster_label %in% keepclasses)]

datlist2=datlist
for (ii in 1:length(datlist2)) {
  rownames(datlist2[[ii]])=tolower(rownames(datlist2[[ii]]))
}
obj.list = list()
obj.list[["macaque1"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac1]),meta.data=allcl_mac[keepcells_mac1,])
obj.list[["macaque2"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac2]),meta.data=allcl_mac[keepcells_mac2,])
obj.list[["human"]]=CreateSeuratObject(round(datlist2[["human"]][,keepcells_hum]),meta.data=allcl_hum[keepcells_hum,])
obj.list[["mouse"]]=CreateSeuratObject(round(datlist2[["mouse"]][,keepcells_mus]),meta.data=allcl_mus[keepcells_mus,])
for (ii in 1:length(obj.list)) {
  obj.list[[ii]] <- SCTransform(obj.list[[ii]], verbose = FALSE)
}
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

require(ggplot2)
require(cowplot)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
uniqcl=unique(obj.integrated@meta.data$cluster_label)
uniqcl=uniqcl[order(uniqcl)]
colvec=obj.integrated@meta.data$cluster_color[match(uniqcl,obj.integrated@meta.data$cluster_label)]
pdf("human_macaque_mouse_allglut_integrated_umap_dimensions_20200407.pdf",width=15)
for (dimval in c(5,10,15,20,25,30)) {
  #dimval=30
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
  p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T,split.by="species",cols=colvec)
  print(p1)
}
dev.off()


###GABAergic Cells###
classkeep=c("GA")
keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Macaque"])
keepcells_mac1=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q18"))]

keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Human"])
keepcells_hum=allcl_hum$sample_name[which(allcl_hum$cluster_label %in% keepclasses)]

keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Mouse"])
keepcells_mus=allcl_mus$sample_name[which(allcl_mus$cluster_label %in% keepclasses)]

datlist2=datlist
for (ii in 1:length(datlist2)) {
  rownames(datlist2[[ii]])=tolower(rownames(datlist2[[ii]]))
}
obj.list = list()
obj.list[["macaque1"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac1]),meta.data=allcl_mac[keepcells_mac1,])
obj.list[["macaque2"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac2]),meta.data=allcl_mac[keepcells_mac2,])
obj.list[["human"]]=CreateSeuratObject(round(datlist2[["human"]][,keepcells_hum]),meta.data=allcl_hum[keepcells_hum,])
obj.list[["mouse"]]=CreateSeuratObject(round(datlist2[["mouse"]][,keepcells_mus]),meta.data=allcl_mus[keepcells_mus,])
for (ii in 1:length(obj.list)) {
  obj.list[[ii]] <- SCTransform(obj.list[[ii]], verbose = FALSE)
}
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30,k.filter=50)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

require(ggplot2)
require(cowplot)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
uniqcl=unique(obj.integrated@meta.data$cluster_label)
uniqcl=uniqcl[order(uniqcl)]
colvec=obj.integrated@meta.data$cluster_color[match(uniqcl,obj.integrated@meta.data$cluster_label)]
pdf("human_macaque_mouse_allgaba_integrated_umap_dimensions_20200407.pdf",width=15)
for (dimval in c(5,10,15,20,25,30)) {
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
  p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T,split.by="species",cols=colvec)
  print(p1)
}
dev.off()


###Magno/Parvo/Konio and LGd/LP
keepclasses=unique(as.character(grep("K|MP|LP|LGD",cluster_palette$cluster_label,val=T)))
keepcells_mac1=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(which(allcl_mac$cluster_label %in% keepclasses),which(allcl_mac$donor=="Q18"))]

keepcells_hum=allcl_hum$sample_name[which(allcl_hum$cluster_label %in% keepclasses)]

keepcells_mus=allcl_mus$sample_name[which(allcl_mus$cluster_label %in% keepclasses)]

datlist2=datlist
for (ii in 1:length(datlist2)) {
  rownames(datlist2[[ii]])=tolower(rownames(datlist2[[ii]]))
}
obj.list = list()
obj.list[["macaque1"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac1]),meta.data=allcl_mac[keepcells_mac1,])
obj.list[["macaque2"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac2]),meta.data=allcl_mac[keepcells_mac2,])
obj.list[["human"]]=CreateSeuratObject(round(datlist2[["human"]][,keepcells_hum]),meta.data=allcl_hum[keepcells_hum,])
obj.list[["mouse"]]=CreateSeuratObject(round(datlist2[["mouse"]][,keepcells_mus]),meta.data=allcl_mus[keepcells_mus,])
for (ii in 1:length(obj.list)) {
  obj.list[[ii]] <- SCTransform(obj.list[[ii]], verbose = FALSE)
}
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

require(ggplot2)
require(cowplot)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
uniqcl=unique(obj.integrated@meta.data$cluster_label)
uniqcl=uniqcl[order(uniqcl)]
colvec=obj.integrated@meta.data$cluster_color[match(uniqcl,obj.integrated@meta.data$cluster_label)]
pdf("human_macaque_mouse_konioparvomagnolgdlp_integrated_umap_dimensions_20200407.pdf",width=15)
for (dimval in c(5,10,15,20,25,30)) {
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
  p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T,split.by="species",cols=colvec)
  print(p1)
}
dev.off()




