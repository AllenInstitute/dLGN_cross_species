library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)


human_dat=read_feather("human_data.feather")
human_meta=read_feather("human_anno.feather")


source("doublet.finder.r")
require(Matrix)
human_dat2=Matrix(as.matrix(human_dat[,-ncol(human_dat)]))
rownames(human_dat2)=human_dat[,ncol(human_dat)][[1]]
alldoublet <- doubletFinder(t(human_dat2), select.genes=rownames(t(human_dat2)), 
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

save(alldoublet,file="doublet.score_lgn_human_all.rda")

keepcells1=intersect(names(alldoublet)[alldoublet<0.2],rownames(human_dat2))
human_dat2@x=human_dat2@x/rep.int(Matrix::colSums(human_dat2), diff(human_dat2@p))*10^6
                     
metadat=as.matrix(human_meta)
metadat=metadat[match(keepcells1,metadat[,1]),]
rownames(metadat)=metadat[,1]
m.int <- CreateSeuratObject(counts = t(human_dat2[keepcells1,]), meta.data = data.frame(metadat,stringsAsFactors = F))
#m.int <- AddMetaData(m.int, metadata = cl.lab, col.name = "cluster")
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
DotPlot(m.int,features=c("GAD2","SLC17A6","SNAP25","GJA1","AIF1","OLIG2","OPALIN"))
#FindMarkers(m.int,ident.1=1)
table(m.int@active.ident,m.int@meta.data$roi)

#keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% 0:4]
keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% c(0,2,3,5)]
datval=human_dat2[keep.noglia,]
metadat2=metadat[keep.noglia,]
m.int.noglia <- CreateSeuratObject(counts = t(datval), meta.data = data.frame(metadat2))
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
DotPlot(m.int.noglia,features=c("GAD2","SLC17A6","SNAP25","BTNL9","HLA-DRB1"))


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

