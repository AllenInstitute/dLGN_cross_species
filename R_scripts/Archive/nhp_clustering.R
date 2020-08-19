library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)

load("20191206_LGN_Mmul_10_star_exon.Rdata")
load("20191206_LGN_Mmul_10_star_intron.Rdata")
mac_dat=exon+intron
load("20191206_LGN_Mmul_10_star_samp.dat.comb.Rdata")
mac_dat=sweep(mac_dat,2,colSums(mac_dat),"/")*10^6
mac_dat2=Matrix(as.matrix(mac_dat))
rownames(mac_dat2)=rownames(mac_dat)
colnames(mac_dat2)=colnames(mac_dat)

source("doublet.finder.r")
require(Matrix)
alldoublet <- doubletFinder(mac_dat2, select.genes=rownames(mac_dat2), 
                                               proportion.artificial = 0.2)
#samp.dat$doublet.score.new <- doublet.score[row.names(samp.dat)]
#doublet.df <- data.frame(score=c(doublet.score, samp.dat$doublet.score.new),
#                         cat = rep(c("new", "old"), rep(nrow(samp.dat), 2)))
# g <- ggplot(doublet.df, aes(x=score, color=cat)) + 
#   stat_ecdf(geom="step")
# ggsave(g, file="score.ecdf.pdf", width= 6, height = 6)
# 
# g <- ggplot(data=samp.dat, aes(x = doublet.score.new, y=Genes.Detected.CPM)) + 
#   geom_point(alpha=0.1, size=0.5)
# ggsave(g, file="Gene_vs_doubletScore.pdf", width= 6, height = 6)

save(alldoublet,file="doublet.score_lgn_macaque_all.rda")

keepcells1=intersect(names(alldoublet)[alldoublet<0.2],colnames(mac_dat2))

metadat=data.frame(samp.dat)
metadat=metadat[keepcells1,]
m.int <- CreateSeuratObject(counts = mac_dat2[,keepcells1], meta.data = data.frame(metadat))
#m.int <- AddMetaData(m.int, metadata = cl.lab, col.name = "cluster")
###Integrate over donors####

m.int <- NormalizeData(object = m.int)
m.int <- ScaleData(object = m.int, verbose = FALSE)
m.int <- FindVariableFeatures(object = m.int, selection.method = "vst", 
                              top.genes = 2000)
m.int <- RunPCA(object = m.int, npcs = 30, verbose = FALSE)
set.seed(0)
m.int <- RunTSNE(object = m.int, reduction = "pca", dims = 1:30)
m.int <- RunUMAP(object = m.int, reduction = "pca", dims = 1:30)
m.int <- FindNeighbors(object = m.int, k.param = 8)
m.int <- FindClusters(object = m.int, resolution = 0.2)
DimPlot(m.int,reduction = "umap",label = TRUE)
DotPlot(m.int,features=c("GAD2","SLC17A6","SNAP25","AQP4","AIF1","OLIG2","OPALIN"))
#FindMarkers(m.int,ident.1=1)
table(m.int@active.ident,m.int@meta.data$roi)

#keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% 0:4]
keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% c(0:7)]
datval=mac_dat2[,keep.noglia]
metadat2=metadat[keep.noglia,]
m.int.noglia <- CreateSeuratObject(counts = datval, meta.data = data.frame(metadat2))
#m.int.noglia <- AddMetaData(m.int.noglia, metadata = cl.lab, col.name = "cluster")
m.int.noglia <- NormalizeData(object = m.int.noglia)
m.int.noglia <- ScaleData(object = m.int.noglia, verbose = FALSE)
m.int.noglia <- FindVariableFeatures(object = m.int.noglia, selection.method = "vst", 
                              top.genes = 2000)
m.int.noglia <- RunPCA(object = m.int.noglia, npcs = 30, verbose = FALSE)
m.int.noglia <- RunTSNE(object = m.int.noglia, reduction = "pca", dims = 1:30)
m.int.noglia <- RunUMAP(object = m.int.noglia, reduction = "pca", dims = 1:30)
m.int.noglia <- FindNeighbors(object = m.int.noglia, k.param = 8)
m.int.noglia <- FindClusters(object = m.int.noglia, resolution = 0.2)
DimPlot(m.int.noglia,reduction = "umap",label = TRUE)
DotPlot(m.int.noglia,features=c("GAD2","SLC17A6","SNAP25","PRKCG","SGCD","CASQ2","IL15RA","KCNH3","SYNDIG1L"))


#keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% 0:4]
keep.noglia.glut = colnames(m.int.noglia@assays$RNA@counts)[m.int.noglia@active.ident %in% c(0,3,4,5)]
datval=human_dat2[keep.noglia.glut,]
metadat2=metadat[keep.noglia.glut,]
m.int.noglia.glut <- CreateSeuratObject(counts = t(datval), meta.data = data.frame(metadat2))
#m.int.noglia.glut <- AddMetaData(m.int.noglia.glut, metadata = cl.lab, col.name = "cluster")
m.int.noglia.glut <- NormalizeData(object = m.int.noglia.glut)
m.int.noglia.glut <- ScaleData(object = m.int.noglia.glut, verbose = FALSE)
m.int.noglia.glut <- FindVariableFeatures(object = m.int.noglia.glut, selection.method = "vst", 
                                     top.genes = 2000)
m.int.noglia.glut <- RunPCA(object = m.int.noglia.glut, npcs = 30, verbose = FALSE)
ElbowPlot(m.int.noglia.glut)
m.int.noglia.glut <- RunTSNE(object = m.int.noglia.glut, reduction = "pca", dims = 1:10)
m.int.noglia.glut <- RunUMAP(object = m.int.noglia.glut, reduction = "pca", dims = 1:10)
m.int.noglia.glut <- FindNeighbors(object = m.int.noglia.glut, k.param = 8,dim=1:10)
m.int.noglia.glut <- FindClusters(object = m.int.noglia.glut, resolution = 0.2)
DimPlot(m.int.noglia.glut,reduction = "umap",label = TRUE)
DotPlot(m.int.noglia.glut,features=c("GAD2","SLC17A6","SNAP25","BTNL9","HLA-DRB1"))




pdf("cluster_allcells_nodoublets_20190430.pdf",useDingbats=F)
for (red in c("umap","tsne")) {
  print(DimPlot(object = m.int, group.by = "roi",reduction=red))
  print(DimPlot(object = m.int, group.by = "ident",reduction=red))
  print(DimPlot(object = m.int, group.by = "external_donor_name",reduction=red))
  print(FeaturePlot(object = m.int, features = c("Snap25","Gad1","Slc17a6","Aqp4","Olig1","Ctss"),reduction = red))
}
for (red in c("umap","tsne")) {
  print(DimPlot(object = m.int.noglia, group.by = "roi",reduction=red))
  print(DimPlot(object = m.int.noglia, group.by = "ident",reduction=red))
  print(DimPlot(object = m.int.noglia, group.by = "external_donor_name",reduction=red))
  print(FeaturePlot(object = m.int.noglia, features = c("Snap25","Gad1","Slc32a1","Slc17a6"),reduction=red))
}
dev.off()

save(m.int,m.int.noglia,file="mouse_seurat_20191222.rda")


#########end#############3

