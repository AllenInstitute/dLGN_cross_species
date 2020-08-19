library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)


pdf("mouse_tsne_plots.pdf",useDingbats=F,width=10)
#### Add mouse analysis ####
load("20180105_RSC-005-049_mouse_LGN_star_cpm.Rdata")
lgn_cpmR=cpmR
load("20180105_RSC-005-049_mouse_LGN_star_samp.dat.Rdata")
load("LP_ss2v4_RSC191.RDA")
lp_cpmR=cpmR
rm(cpmR)

source("doublet.finder.r")
require(Matrix)
doublet.score <- doubletFinder(Matrix(lgn_cpmR), select.genes=rownames(lgn_cpmR), 
                                               proportion.artificial = 0.2)
save(doublet.score, file="doublet.score_lgn_mouse.rda")

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

doublet.score.lp <- doubletFinder(Matrix(lp_cpmR), select.genes=rownames(lp_cpmR), 
                               proportion.artificial = 0.2)
#save(doublet.score.lp, file="doublet.score_lp.rda")

samp.dat.lp$doublet.score.new <- doublet.score.lp[row.names(samp.dat.lp)]

alldoublet=doubletFinder(Matrix(cbind(lgn_cpmR,lp_cpmR)), select.genes=rownames(lgn_cpmR), 
                                          proportion.artificial = 0.2)
save(alldoublet,file="doublet.score_lgn_mouse_all.rda")

keepcells1=intersect(intersect(names(alldoublet)[alldoublet<0.2],rownames(samp.dat)),names(doublet.score[doublet.score<0.2]))
keepcells2=intersect(intersect(names(alldoublet)[alldoublet<0.2],rownames(samp.dat.lp)),names(doublet.score.lp)[doublet.score.lp<0.2])

common.cols=intersect(colnames(samp.dat),colnames(samp.dat.lp))
metadat=rbind(samp.dat[keepcells1,common.cols],samp.dat.lp[keepcells2,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
m.int <- CreateSeuratObject(counts = cbind(lgn_cpmR[,keepcells1],lp_cpmR[,keepcells2]), meta.data = metadat)
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
DotPlot(m.int,features=c("Gad2","Slc17a7","Snap25","Gja1","Ctss","Olig2"))
table(m.int@active.ident,m.int@meta.data$roi)

#keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% 0:4]
keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% c(0,1,2,3,5)]
datval=cbind(lgn_cpmR,lp_cpmR)
datval=datval[,keep.noglia]
metadat=rbind(samp.dat[,common.cols],samp.dat.lp[,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
metadat=metadat[match(keep.noglia,metadat$exp_component_name),]
m.int.noglia <- CreateSeuratObject(counts = datval, meta.data = metadat)
m.int.noglia <- AddMetaData(m.int.noglia, metadata = cl.lab, col.name = "cluster")
m.int.noglia <- NormalizeData(object = m.int.noglia)
m.int.noglia <- ScaleData(object = m.int.noglia, verbose = FALSE)
m.int.noglia <- FindVariableFeatures(object = m.int.noglia, selection.method = "vst", 
                              top.genes = 2000)
m.int.noglia <- RunPCA(object = m.int.noglia, npcs = 30, verbose = FALSE)
m.int.noglia <- RunTSNE(object = m.int.noglia, reduction = "pca", dims = 1:30)
m.int.noglia <- RunUMAP(object = m.int.noglia, reduction = "pca", dims = 1:30)
m.int.noglia <- FindNeighbors(object = m.int.noglia, k.param = 8)
m.int.noglia <- FindClusters(object = m.int.noglia, resolution = 0.2)
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


####Donor specific values####
pdf("donor_specific_20190430.pdf",useDingbats = F,width=10)
keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% c(0,1,2,3,5)]
donorkeep=c("207722","209943","211623","211624","213715","213716")
datval=cbind(lgn_cpmR,lp_cpmR)
metadat=rbind(samp.dat[,common.cols],samp.dat.lp[,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
for (ii in donorkeep) {
  keepcols=intersect(keep.noglia,metadat$exp_component_name[metadat$external_donor_name==as.numeric(ii)])
  metadat_d=metadat[match(keepcols,metadat$exp_component_name),]
  datval_d=datval[,keepcols]
  m.int.donor <- CreateSeuratObject(counts = datval_d, meta.data = metadat_d)
  #m.int.donor <- AddMetaData(m.int.donor, metadata = cl.lab, col.name = "cluster")
  m.int.donor <- NormalizeData(object = m.int.donor)
  m.int.donor <- ScaleData(object = m.int.donor, verbose = FALSE)
  m.int.donor <- FindVariableFeatures(object = m.int.donor, selection.method = "vst", 
                                       top.genes = 2000)
  m.int.donor <- RunPCA(object = m.int.donor, npcs = 30, verbose = FALSE)
  m.int.donor <- RunTSNE(object = m.int.donor, reduction = "pca", dims = 1:30,perplexity=10)
  m.int.donor <- RunUMAP(object = m.int.donor, reduction = "pca", dims = 1:30)
  m.int.donor <- FindNeighbors(object = m.int.donor, k.param = 8)
  m.int.donor <- FindClusters(object = m.int.donor, resolution = 0.1)
  for (red in c("umap","tsne")) {
    print(DimPlot(object = m.int.donor, group.by = "roi",reduction=red,main=ii))
    print(FeaturePlot(object = m.int.donor, features = c("Snap25","Gad1","Slc32a1","Slc17a6"),reduction=red))
    print(DimPlot(object = m.int.donor, group.by = "ident",main=ii,reduction=red))
    print(DimPlot(object = m.int.donor, group.by = "external_donor_name",main=ii,reduction=red))
  }
}
dev.off()

####Integrate across donors####
pdf("cross_donor_integration_20190430.pdf",useDingbats = F,width=10)
keep.noglia = colnames(m.int.noglia@assays$RNA@counts)
donorkeep=c("207722","209943","211623","211624","213715","213716")
datval=cbind(lgn_cpmR,lp_cpmR)
metadat=rbind(samp.dat[,common.cols],samp.dat.lp[,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
keepcols=intersect(keep.noglia,metadat$exp_component_name[metadat$external_donor_name %in% as.numeric(donorkeep)])
metadat_d=metadat[match(keepcols,metadat$exp_component_name),]
datval_d=datval[,keepcols]
m.dat <- CreateSeuratObject(counts = datval_d, meta.data = metadat_d)
m.list <- SplitObject(object = m.dat, split.by = "external_donor_name")
for (i in 1:length(x = m.list)) {
  m.list[[i]] <- NormalizeData(object = m.list[[i]], verbose = FALSE)
  m.list[[i]] <- FindVariableFeatures(object = m.list[[i]], 
                                      selection.method = "vst", 
                                      nfeatures = 2000, verbose = FALSE)
}
m.anchors <- FindIntegrationAnchors(object.list = m.list, 
                                    anchor.features = 2000, dims = 1:30,k.filter=30,k.score=20)
m.int <- IntegrateData(anchorset = m.anchors, dims = 1:30, 
                       features.to.integrate = row.names(m.dat@assays$RNA))
DefaultAssay(object = m.int) <- "integrated"

m.int <- ScaleData(object = m.int, verbose = FALSE)
m.int <- RunPCA(object = m.int, npcs = 30, 
                features = m.int@assays$integrated@var.features, 
                verbose = FALSE)
m.int <- RunTSNE(object = m.int, reduction = "pca", dims = 1:30)
m.int <- RunUMAP(object = m.int, reduction = "pca", dims = 1:30)
m.int <- FindNeighbors(object = m.int, k.param = 8)
m.int <- FindClusters(object = m.int, resolution = 0.1)
for (red in c("umap","tsne")) {
  print(DimPlot(object = m.int, group.by = "roi",reduction=red,main=ii))
  print(FeaturePlot(object = m.int, features = c("Snap25","Gad1","Slc32a1","Slc17a6"),reduction=red))
  print(DimPlot(object = m.int, group.by = "ident",main=ii,reduction=red))
  print(DimPlot(object = m.int, group.by = "external_donor_name",main=ii,reduction=red))
}
dev.off()



#####PCA to identify key continua######
keepcells1=intersect(intersect(names(alldoublet)[alldoublet<0.2],rownames(samp.dat)),names(doublet.score[doublet.score<0.2]))
keepcells2=intersect(intersect(names(alldoublet)[alldoublet<0.2],rownames(samp.dat.lp)),names(doublet.score.lp)[doublet.score.lp<0.2])

common.cols=intersect(colnames(samp.dat),colnames(samp.dat.lp))
metadat=rbind(samp.dat[keepcells1,common.cols],samp.dat.lp[keepcells2,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
m.int <- CreateSeuratObject(counts = cbind(lgn_cpmR[,keepcells1],lp_cpmR[,keepcells2]), meta.data = metadat)
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

#keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% 0:4]
keep.noglia = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% c(0,1,2,3,5)]
datval=cbind(lgn_cpmR,lp_cpmR)
datval=datval[,keep.noglia]
metadat=rbind(samp.dat[,common.cols],samp.dat.lp[,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
metadat=metadat[match(keep.noglia,metadat$exp_component_name),]
m.int.noglia <- CreateSeuratObject(counts = datval, meta.data = metadat)
m.int.noglia <- AddMetaData(m.int.noglia, metadata = cl.lab, col.name = "cluster")
m.int.noglia <- NormalizeData(object = m.int.noglia)
m.int.noglia <- ScaleData(object = m.int.noglia, verbose = FALSE)
m.int.noglia <- FindVariableFeatures(object = m.int.noglia, selection.method = "vst", 
                                     top.genes = 2000)
m.int.noglia <- RunPCA(object = m.int.noglia, npcs = 30, verbose = FALSE)
m.int.noglia <- RunTSNE(object = m.int.noglia, reduction = "pca", dims = 1:30)
m.int.noglia <- RunUMAP(object = m.int.noglia, reduction = "pca", dims = 1:30)
m.int.noglia <- FindNeighbors(object = m.int.noglia, k.param = 8)
m.int.noglia <- FindClusters(object = m.int.noglia, resolution = 0.2)

####rerun PCA on large glut cluster###
keep.glut = colnames(m.int@assays$RNA@counts)[m.int@active.ident %in% c(0)]
metadat=rbind(samp.dat[,common.cols],samp.dat.lp[,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
metadat=metadat[match(keep.glut,metadat$exp_component_name),]
metadat=metadat[grep("LGN",metadat$roi),]
datval=cbind(lgn_cpmR,lp_cpmR)
datval=datval[,metadat$exp_component_name]
m.int.glut <- CreateSeuratObject(counts = datval[genecorma], meta.data = metadat)
m.int.glut <- AddMetaData(m.int.glut, metadata = cl.lab, col.name = "cluster")
m.int.glut <- NormalizeData(object = m.int.glut)
m.int.glut <- ScaleData(object = m.int.glut, verbose = FALSE)
m.int.glut <- FindVariableFeatures(object = m.int.glut, selection.method = "vst", nfeatures=5000,
                                     top.genes = 5000)
m.int.glut <- RunPCA(object = m.int.glut, npcs = 10, verbose = FALSE)
m.int.glut <- RunTSNE(object = m.int.glut, reduction = "pca", dims = 1:10)
m.int.glut <- RunUMAP(object = m.int.glut, reduction = "pca", dims = 1:10)
m.int.glut <- FindNeighbors(object = m.int.glut, k.param = 8)
m.int.glut <- FindClusters(object = m.int.glut, resolution = 0.2)
DimPlot(object = m.int.glut, group.by = "roi",reduction=red)
DimPlot(object = m.int.glut, group.by = "ident",reduction=red)

#pctab=m.int.glut@reductions$pca@cell.embeddings
#pcgenes=m.int.glut@reductions$pca@feature.loadings
require(irlba)
vargen=m.int.glut@assays$RNA@var.features
pcvar1=prcomp_irlba(scale(t(log2(datval[vargen,]+1))),n = 20)
pctab=pcvar1$x

keepmetacols=c(1:3,5:18,59,61,65,66,67,68,69,82)
techcor=cor(pctab,metadat[,keepmetacols])
apply(techcor,1,max)

###select first PC with <0.2 correlation : PC 3###
rownames(pcvar1$rotation)=vargen
pcvar1$rotation[order(-abs(pcvar1$rotation[,3]))[1:20],1:3]
apply(pcvar1$rotation,2,function(x){return(sum(abs(x)>abs(x[match("Foxp2",vargen)])))})

###make heatmap###
keepgen=rownames(pcvar1$rotation)[order(-abs(pcvar1$rotation[,3]))[1:100]]
genecor=t(cor(log2(t((as.matrix(datval[grep("Gm|^LOC",rownames(datval),invert=T),])+1))),pctab[,3],method="spearman"))
keepgen=colnames(genecor)[order(-abs(genecor))]
#keepgen=keepgen[keepgen %in% rownames(datval)[rowSums(as.matrix(datval)>0)>0]]
keepgen=keepgen[1:50]

plotmat=datval[keepgen,order(pctab[,3])]
require(gplots)
breakvals=seq(from=-2,to=2,length.out=256)
roicols=c("green","black")[as.numeric(as.factor(metadat$roi[order(pctab[,2])]))]
hclustout=hclust(as.dist(1-cor(scale(t(log2(as.matrix(plotmat)+1))))),method="ward")


pdf("core_shell_pca_genes.pdf",useDingbats=F)
heatmap.2(as.matrix(log2(plotmat+1)),scale="row",trace='none',Colv=F,Rowv=as.dendrogram(hclustout),breaks=breakvals,col=bluered(255),ColSideColors=roicols)
shellhist=hist(which(roicols=="black"),breaks=20,plot=F)
corehist=hist(which(roicols=="green"),breaks=20,plot=F)
plot(shellhist$mids,shellhist$counts/sum(shellhist$counts),col="green",pch='',xaxt='n',ylab='Relative frequency',xlab=NA)
lines(shellhist$mids,shellhist$counts/sum(shellhist$counts),col="black")
lines(corehist$mids,corehist$counts/sum(corehist$counts),col="green")
plotmat2=datval[c("Foxp2","Robo2","Grid2","Ebf1"),order(pctab[,3])]
heatmap.2(as.matrix(log2(plotmat2+1)),scale="row",trace='none',Colv=F,Rowv=F,breaks=breakvals,col=bluered(255),ColSideColors=roicols)
dev.off()



       