######==========================######
### species independent clustering ###
######==========================######


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)
source("R_scripts/doublet.finder.r")
require(Matrix)

load("int_data/collected_data_20200818.rda")



spec="human"


if (spec=="macaque") {
  orthos=read.csv("data/ortholog_table_20191122.csv",as.is=T)
  nonlocgenes=orthos$rhesus_symbol[intersect(grep("^LOC",orthos$human_symbol,invert=T),grep("^LOC",orthos$rhesus_symbol))]
  keepgen=setdiff(rownames(datlist[[spec]]),setdiff(grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|^TRNA",rownames(datlist[[spec]]),val=T),nonlocgenes))
  datlist[[spec]]=datlist[[spec]][keepgen,]
  datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
} else {
  if (spec=="human") {
    keepgen=grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|[0-9]P$|-AS|-PS|DUX4L",rownames(datlist[[spec]]),val=T,invert=T)
    keepdons=which(metalist$human$external_donor_name_label %in% c("H200.1023","H200.1025","H200.1030"))
    datlist[[spec]]=datlist[[spec]][keepgen,keepdons]
    datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
  } else {
    keepgen=grep("^LOC|^mt-|^Rpl|^Rps",rownames(datlist[[spec]]),val=T,invert=T)
  }
}

alldoublet <- doubletFinder(datlist[[spec]], select.genes=keepgen, 
                            proportion.artificial = 0.2)
#save(alldoublet,file="doublet.score_lgn_macaque_all.rda")

keepcells1=intersect(names(alldoublet)[alldoublet<0.2],colnames(datlist[[spec]]))


###step one - isolate neurons###
m.int <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcells1], meta.data = metalist[[spec]][keepcells1,])
m.int <- NormalizeData(object = m.int)
if (spec=="macaque") {
  m.int <- ScaleData(object = m.int, verbose = FALSE,vars.to.regress="donor")
} else {
  m.int <- ScaleData(object = m.int, verbose = FALSE)
}
m.int <- FindVariableFeatures(object = m.int, selection.method = "vst", 
                              nfeatures = 3000) #2000 for mouse and macaque
m.int <- RunPCA(object = m.int, npcs = 30, verbose = FALSE)
set.seed(0)
m.int <- RunTSNE(object = m.int, reduction = "pca", dims = 1:30)
m.int <- RunUMAP(object = m.int, reduction = "pca", dims = 1:30)
m.int <- FindNeighbors(object = m.int, k.param = 8,dims=1:30)
m.int <- FindClusters(object = m.int, resolution = 0.4) ##0.2 for macaque,0.4 for mouse
DimPlot(m.int,reduction = "umap",label = TRUE)
DimPlot(m.int,reduction = "umap",group.by="external_donor_name_label",label = TRUE)
FeaturePlot(m.int,features=c("AIF1","GJA1","FGFR3"))
DotPlot(m.int,features=c("GAD2","SLC17A6","SNAP25","AQP4","AIF1","OLIG2","OPALIN","LAMP5"))
DotPlot(m.int,features=c("Gad2","Slc17a6","Snap25","Aqp4","Ctss","Olig2","Opalin","Lamp5","Pdgfra"))


####no glia clustering###
if (spec=="macaque") {
  keepcells2=names(m.int@active.ident)[m.int@active.ident %in% c(0:9)]
  m.int2=m.int
} else {
  if (spec=="human") {
    #keepcells2=names(m.int@active.ident)[m.int@active.ident %in% c(0,2,3,5,7)]
    keepcells2=names(m.int@active.ident)[m.int@active.ident %in% c(0,2,4,5)]
  }
  if (spec=="mouse") {
    keepcells2=names(m.int@active.ident)[m.int@active.ident %in% c(0,1,2,3,4,5,6,8,9,11,12)]
  }
  m.int2 <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcells2], meta.data = metalist[[spec]][keepcells2,])
  m.int2 <- NormalizeData(object = m.int2)
  m.int2 <- ScaleData(object = m.int2, verbose = FALSE)
  m.int2 <- FindVariableFeatures(object = m.int2, selection.method = "vst", 
                                 nfeatures = 3000) #2000 for mouse and macaque
  m.int2 <- RunPCA(object = m.int2, npcs = 30, verbose = FALSE)
  set.seed(0)
  m.int2 <- RunTSNE(object = m.int2, reduction = "pca", dims = 1:30)
  m.int2 <- RunUMAP(object = m.int2, reduction = "pca", dims = 1:30)
  m.int2 <- FindNeighbors(object = m.int2, k.param = 8,dims=1:30)
  m.int2 <- FindClusters(object = m.int2, resolution = 0.5) #0.5 for mouse and macaque
}
DimPlot(m.int2,reduction = "umap",label = TRUE)
DimPlot(m.int2,reduction = "umap",label = TRUE,group.by="external_donor_name_label")
DotPlot(m.int2,features=c("GAD2","SLC17A6","SNAP25","AQP4","AIF1","OLIG2","OPALIN","BTNL9","FOXP2","PRKCG"))
FeaturePlot(m.int2,features=c("GAD2","SLC17A6","SNAP25","AQP4","AIF1","OLIG2","OPALIN"))
FeaturePlot(m.int2,features=c("Gad2"))
DotPlot(m.int2,features=c("Gad2","Slc17a6","Snap25","Aqp4","Ctss","Olig2","Opalin","Lamp5","Pdgfra"))
ttt=FindMarkers(m.int2,ident.1=4,ident.2=0)

###glutamatergic neurons###
if (spec=="human") {
  keepcells3=names(m.int2@active.ident)[m.int2@active.ident %in% c(0,1)]
}
if (spec=="macaque") {
  keepcells3=names(m.int2@active.ident)[m.int2@active.ident %in% c(0,1,3,4,5,6,8)]
}
if (spec=="mouse") {
  keepcells3=names(m.int2@active.ident)[m.int2@active.ident %in% c(0,1,3,7,12)]
}
m.int3 <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcells3], meta.data = metalist[[spec]][keepcells3,])
m.int3 <- NormalizeData(object = m.int3)
if (spec=="macaque") {
  m.int3 <- ScaleData(object = m.int3, verbose = FALSE,vars.to.regress = "donor")
  dimval=30
} else {
  m.int3 <- ScaleData(object=m.int3,verbose=FALSE)
  if (spec=="mouse") {
    dimval=25
  } else {
    dimval=30
  }
}
m.int3 <- FindVariableFeatures(object = m.int3, selection.method = "vst", 
                               nfeatures = 3000) #2000 for macaque and mouse
m.int3 <- RunPCA(object = m.int3, npcs = 30, verbose = FALSE)
set.seed(0)
m.int3 <- RunTSNE(object = m.int3, reduction = "pca", dims = 1:dimval)
m.int3 <- RunUMAP(object = m.int3, reduction = "pca", dims = 1:dimval)
m.int3 <- FindNeighbors(object = m.int3, k.param = 8,dim=1:dimval)
m.int3 <- FindClusters(object = m.int3, resolution = 0.4,dims=1:dimval) #0.2 macaque, 0.1 for mouse 
DimPlot(m.int3,reduction = "umap",label = TRUE)
DimPlot(m.int3,reduction="umap",group.by="roi")
DotPlot(m.int3,reduction="umap",features=c("Snap25","Slc17a6","Gad2","Lamp5",""))
DimPlot(m.int3,reduction = "umap",label = TRUE,group.by="external_donor_name_label")
DotPlot(m.int3,features=c("SNAP25","AIF1","AQP4"))
DotPlot(m.int3,features=c("GAD2","SLC17A6","SNAP25","PRKCG","SGCD","CASQ2","IL15RA","KCNH3","SYNDIG1L","EBF1","FOXP2","LHX9","CD5","BTNL9"))
#head(FindMarkers(m.int3,ident.1=5),n=20)
table(paste0(m.int3@meta.data$external_donor_name,m.int3@meta.data$roi),m.int3@active.ident)
table(paste0(m.int3@meta.data$external_donor_name_label,m.int3@meta.data$roi_label),m.int3@active.ident)


###sub group: non-konio/lgd
if (spec=="macaque") {
  keepcells4=names(m.int3@active.ident)[m.int3@active.ident %in% c(0,3,4,6)]
}
if (spec=="human") {
  keepcells4=names(m.int2@active.ident)[m.int2@active.ident %in% c(0)]
}
if (spec=="mouse") {
  keepcells4=names(m.int3@active.ident)[m.int3@active.ident %in% c(0,2)]
}
m.int4 <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcells4], meta.data = metalist[[spec]][keepcells4,])
m.int4 <- NormalizeData(object = m.int4)
if (spec=="macaque") {
  m.int4 <- ScaleData(object = m.int4, verbose = FALSE,vars.to.regress = "donor")
} else {
  if (spec %in% c("mouse")) {
    #m.int4 <- ScaleData(object = m.int4, verbose = FALSE,vars.to.regress = "external_donor_name_label")
    m.int4 <- ScaleData(object = m.int4, verbose = FALSE)
  } else {
    m.int4 <- ScaleData(object = m.int4, verbose=FALSE,vars.to.regress = "external_donor_name_label")
  }
}
m.int4 <- FindVariableFeatures(object = m.int4, selection.method = "vst", 
                               nfeatures = 3000)
m.int4 <- RunPCA(object = m.int4, npcs = 30, verbose = FALSE)
set.seed(0)
m.int4 <- RunTSNE(object = m.int4, reduction = "pca", dims = 1:30)
m.int4 <- RunUMAP(object = m.int4, reduction = "pca", dims = 1:30)
m.int4 <- FindNeighbors(object = m.int4, k.param = 8, dims = 1:30)
m.int4 <- FindClusters(object = m.int4, resolution = 0.7) ###0.4 for human, 0.2 for macaque,0.2 for mouse
DimPlot(m.int4,reduction = "umap",label = TRUE)
DimPlot(m.int4,reduction = "umap",group.by="roi_label",cols=rep(c("green","red","blue"),each=2))
DimPlot(m.int4,reduction = "umap",label = TRUE,group.by="external_donor_name_label")
DotPlot(m.int4,features=c("GAD2","SLC17A6","SNAP25","FOXP2","CRH","SYT2","BTNL9","PRKCH"))
FeaturePlot(m.int4,features=c("FOXP2","CRH","BTNL9","PRKCH"))


pdf("mouse_dim_range.pdf") 
for (numgen in c(1000,2000,3000)) {  
  m.int4 <- FindVariableFeatures(object = m.int4, selection.method = "vst", 
                                 top.genes = numgen)
  m.int4 <- RunPCA(object = m.int4, npcs = 30, verbose = FALSE)
  for (dimval in 3:30) {
    set.seed(0)
    m.int4 <- RunUMAP(object = m.int4, reduction = "pca", dims = 1:dimval)
    p1=DimPlot(m.int4,reduction = "umap",group.by="roi")
    print(p1)
  }
  print(c(numgen,dimval))
}
dev.off()



###sub group: konio/LP
if (spec=="macaque") {
  keepcells5=names(m.int3@active.ident)[m.int3@active.ident %in% c(1)]  ###otherwise 2
}
if (spec=="human") {
  keepcells5=names(m.int2@active.ident)[m.int2@active.ident %in% c(1)]  ###otherwise 2
}
if (spec=="mouse") {
  keepcells5=names(m.int3@active.ident)[m.int3@active.ident %in% c(1)]  ###otherwise 2
}
m.int5 <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcells5], meta.data = metalist[[spec]][keepcells5,])
m.int5 <- NormalizeData(object =m.int5)
if (spec=="macaque") {
  m.int5 <- ScaleData(object =m.int5, verbose = FALSE,vars.to.regress = "donor")
} else {
  if (spec %in% c("mouse")) {
    #m.int5 <- ScaleData(object =m.int5, verbose = FALSE,vars.to.regress = "external_donor_name_label")
    m.int5 <- ScaleData(object =m.int5, verbose = FALSE)
  } else {
    m.int5 <- ScaleData(object =m.int5, verbose = FALSE)
  }
}
m.int5 <- FindVariableFeatures(object =m.int5, selection.method = "vst", 
                               top.genes = 2000) #2000 for macaque
m.int5 <- RunPCA(object =m.int5, npcs = 30, verbose = FALSE)
set.seed(0)
m.int5 <- RunTSNE(object =m.int5, reduction = "pca", dims = 1:30)
m.int5 <- RunUMAP(object =m.int5, reduction = "pca", dims = 1:5) #30 for macaque and mouse, 5 for human
m.int5 <- FindNeighbors(object =m.int5, k.param = 8, dims = 1:30)
m.int5 <- FindClusters(object =m.int5, resolution = 0.9) ##0.2 for macaque
DimPlot(m.int5,reduction = "umap",label = TRUE)
DimPlot(m.int5,reduction = "umap",label = TRUE,group.by="external_donor_name_label")
DotPlot(m.int5,features=c("GAD2","SLC17A6","SNAP25","PRKCG","SGCD","CASQ2","IL15RA","KCNH3","SYNDIG1L","EBF1","FOXP2","AQP4","KCNJ2","HAPLN2","CWH43","MLPH"))
DotPlot(m.int5,features=c("GAD2","SLC17A6","SNAP25","PRKCG","SGCD","BCHE","GRB14"))
allmarks=FindAllMarkers(m.int5)
FeaturePlot(m.int5,features=c("PRKCG","BCHE","GRB14"))
table(paste0(m.int5@meta.data$external_donor_name,m.int5@meta.data$roi),m.int5@active.ident)

pdf("human_dim_range.pdf") 
for (numgen in c(1000,2000,3000)) {  
  m.int5 <- FindVariableFeatures(object = m.int5, selection.method = "vst", 
                                 nfeatures=numgen)
  m.int5 <- RunPCA(object = m.int5, npcs = 30, verbose = FALSE)
  for (dimval in 3:30) {
    set.seed(0)
    m.int5 <- RunUMAP(object = m.int5, reduction = "pca", dims = 1:dimval)
    m.int5 <- FindClusters(object =m.int5, resolution = 0.6) ##0.2 for macaque
    p1=DimPlot(m.int5,reduction = "umap",label=T)
    print(p1)
    p1=FeaturePlot(m.int5,features=c("PENK","BCHE"))
    print(p1)
  }
  print(c(numgen,dimval))
}
dev.off()




###subgroup - GABAergic1
if (spec=="macaque") {
  keepcells6=names(m.int2@active.ident)[m.int2@active.ident %in% c(2,7)]
}
if (spec=="human") {
  keepcells6=names(m.int2@active.ident)[m.int2@active.ident %in% c(3)]  ##other one is 2
}

if (spec=="mouse") {
  keepcells6=names(m.int2@active.ident)[m.int2@active.ident %in% c(2,4)]
}


m.int6 <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcells6], meta.data = metalist[[spec]][keepcells6,])
m.int6 <- NormalizeData(object =m.int6)
if (spec=="macaque") {
  m.int6 <- ScaleData(object =m.int6, verbose = FALSE,vars.to.regress = "donor")
} else {
  if (spec %in% c("human","mouse")) {
    #m.int6 <- ScaleData(object =m.int6, verbose = FALSE,vars.to.regress = "external_donor_name_label")
    m.int6 <- ScaleData(object =m.int6, verbose = FALSE)
  }
}
m.int6 <- FindVariableFeatures(object =m.int6, selection.method = "vst", 
                               top.genes = 2000)
m.int6 <- RunPCA(object =m.int6, npcs = 30, verbose = FALSE)
set.seed(0)
ElbowPlot(m.int6)
if (spec=="macaque") {
  m.int6 <- RunTSNE(object =m.int6, reduction = "pca", dims = 1:10) ##15 for macaque
  m.int6 <- RunUMAP(object =m.int6, reduction = "pca", dims = 1:10)
  m.int6 <- FindNeighbors(object =m.int6, k.param = 8, dims = 1:10)
  m.int6 <- FindClusters(object =m.int6, resolution = 0.2)
  DotPlot(m.int6,features=c("GAD2","SLC17A6","SNAP25","LAMP5","TRPC4","KRT80","DSEL","VGF","OXTR","GBP3","XIST","FBN2"))
  table(m.int6@meta.data$external_donor_name,m.int6@active.ident)
  ###merge clusters 4,1
  #  m.int6@active.ident[m.int6@active.ident==4]=1
}
if (spec=="human") {
  dimval=8
  #m.int6 <- RunTSNE(object =m.int6, reduction = "pca", dims = 1:10)
  #m.int6 <- RunUMAP(object =m.int6, reduction = "pca", dims = 1:dimval)
  m.int6 <- FindNeighbors(object =m.int6, k.param = 3, dims = 1:dimval)
  m.int6 <- FindClusters(object =m.int6, resolution = 0.3)
  #DimPlot(m.int6,reduction="umap",label=T)
  #DotPlot(m.int6,features=c("GAD2","SLC17A6","SNAP25","LAMP5","TRPC4","NTRK1","CTXN3","ALCAM"))
  #FeaturePlot(m.int6,features=c("GAD2","LAMP5","TRPC4"))
  print(table(m.int6@meta.data$cluster_label,m.int6@active.ident))   ###separation into clusters 0-2, 3
}
if (spec=="mouse") {
  dimval=20
  m.int6 <- RunTSNE(object =m.int6, reduction = "pca", dims = 1:dimval)
  m.int6 <- RunUMAP(object =m.int6, reduction = "pca", dims = 1:dimval)
  m.int6 <- FindNeighbors(object =m.int6, k.param = 8, dims = 1:dimval)
  m.int6 <- FindClusters(object =m.int6, resolution = 0.2)
  DimPlot(m.int6,reduction="umap",label=T)
  DimPlot(m.int6,reduction="umap",group.by="roi")
  DotPlot(m.int6,features=c("Gad2","Snap25","Slc17a6","Lamp5","Chrna6","Plac9a","Pax7","Dscaml1","Pvalb","Spata20","Rspo2","Dlx5","Cbln2","Shox2","Htr1a","Sox14"))
  print(table(m.int6@meta.data$cluster_label,m.int6@active.ident))
}



###subgroup - GABAergic2
if (spec=="human") {
  keepcells7=names(m.int2@active.ident)[m.int2@active.ident %in% c(2)]  ##other one is 3
}
if (spec=="mouse") {
  keepcells7=names(m.int2@active.ident)[m.int2@active.ident %in% c(5,6,8,9,10,11)]
}
m.int7 <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcells7], meta.data = metalist[[spec]][keepcells7,])
m.int7 <- NormalizeData(object =m.int7)
if (spec %in% c("human","mouse")) {
  #m.int7 <- ScaleData(object =m.int7, verbose = FALSE,vars.to.regress = "external_donor_name_label")
  m.int7 <- ScaleData(object =m.int7, verbose = FALSE)
}
m.int7 <- FindVariableFeatures(object =m.int7, selection.method = "vst", 
                               top.genes = 2000)
m.int7 <- RunPCA(object =m.int7, npcs = 30, verbose = FALSE)
set.seed(0)
if (spec=="human") {
  dimval=6
  #m.int7 <- RunTSNE(object =m.int7, reduction = "pca", dims = 1:10) ##15 for macaque
  #m.int7 <- RunUMAP(object =m.int7, reduction = "pca", dims = 1:dimval)
  m.int7 <- FindNeighbors(object =m.int7, k.param = 3, dims = 1:dimval)
  m.int7 <- FindClusters(object =m.int7, resolution = 0.1)
  #DimPlot(m.int7,reduction="umap",label=T)
  DotPlot(m.int7,features=c("GAD2","SLC17A6","SNAP25","LAMP5","TRPC4","NTRK1","CTXN3","ALCAM"))
  #FeaturePlot(m.int7,features=c("GAD2","LAMP5","TRPC4"))
  print(table(m.int7@meta.data$cluster_label,m.int7@active.ident))
}
if (spec=="mouse") {
  dimval=20
  m.int7 <- RunTSNE(object =m.int7, reduction = "pca", dims = 1:10)
  m.int7 <- RunUMAP(object =m.int7, reduction = "pca", dims = 1:dimval)
  m.int7 <- FindNeighbors(object =m.int7, k.param = 8, dims = 1:dimval)
  m.int7 <- FindClusters(object =m.int7, resolution = 0.2)
  DimPlot(m.int7,reduction="umap",label=T)
  DotPlot(m.int7,features=c("Gad2","Dlx5","Cbln2","Pvalb","Plek2","Htr1a","Spata20","Rspo2","Irx1","Eltd1","Rspo3"))
  DotPlot(m.int7,features=c("Gad2","Snap25","Slc17a6","Lamp5","Chrna6","Plac9a","Pax7","Dscaml1","Pvalb","Spata20","Rspo2","Dlx5","Cbln2","Shox2","Htr1a","Sox14"))
  #FeaturePlot(m.int7,features=c("GAD2","LAMP5","TRPC4"))
  print(table(m.int7@meta.data$cluster_label,m.int7@active.ident))
}
bbb=FindAllMarkers(m.int7)




if (spec=="human") {
  save(m.int,m.int2,m.int3,m.int4,m.int5,m.int6,m.int7,file="human_clusters_subclusters_20200126.rda")
}
if (spec=="macaque") {
  save(m.int,m.int2,m.int3,m.int4,m.int5,m.int6,file="macaque_clusters_subclusters.rda")
}
if (spec=="mouse") {
  save(m.int,m.int2,m.int3,m.int4,m.int5,m.int6,m.int7,file="mouse_clusters_subclusters.rda")
}



####human - donor specific clustering###
load("human_clusters_subclusters.rda")
donorvals=unique(m.int4@meta.data$external_donor_name_label)
table(m.int4@meta.data$external_donor_name_label,m.int4@meta.data$roi_label)
for (ii in 1:length(donorvals)) {
  don=donorvals[ii]
  keepcellsd=names(m.int4@active.ident)[m.int4@meta.data$external_donor_name_label==don]
  m.int_don <- CreateSeuratObject(counts = datlist[[spec]][keepgen,keepcellsd], meta.data = metalist[[spec]][keepcellsd,])
  m.int_don <- NormalizeData(object =m.int_don)
  m.int_don <- ScaleData(object =m.int_don, verbose = FALSE)
  m.int_don <- FindVariableFeatures(object =m.int_don, selection.method = "vst", 
                                    top.genes = 2000)
  m.int_don <- RunPCA(object =m.int_don, npcs = 30, verbose = FALSE)
  set.seed(0)
  pdf(paste0("human_",don,".pdf"),useDingbats=F)
  for (dimval in 3:20) {
    m.int_don <- RunUMAP(object =m.int_don, reduction = "pca", dims = 1:dimval)
    m.int_don <- FindNeighbors(object =m.int_don, k.param = 8, dims = 1:dimval)
    m.int_don <- FindClusters(object =m.int_don, resolution = 0.5)
    p=DimPlot(m.int_don,reduction="umap",label=T)
    print(p)
    p=DimPlot(m.int_don,reduction="umap",label=T,group.by="roi_label")
    print(p)
  }
  p=FeaturePlot(m.int_don,features=c("FOXP2","EBF1"))
  print(p)
  dev.off()
}
