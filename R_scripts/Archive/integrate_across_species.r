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
load("human_clusters_subclusters_20200126.rda")

load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl_20200126.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")

####integrate populations####
#grepchar="MP"
grepchar="MP|K"
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
obj.list = list()
obj.list[["macaque1"]]=CreateSeuratObject(round(datlist[["macaque"]][,keepcells_mac1]),meta.data=allcl_mac[keepcells_mac1,])
obj.list[["macaque2"]]=CreateSeuratObject(round(datlist[["macaque"]][,keepcells_mac2]),meta.data=allcl_mac[keepcells_mac2,])
obj.list[["human"]]=CreateSeuratObject(round(datlist[["human"]][,keepcells_hum]),meta.data=allcl_hum[keepcells_hum,])
for (ii in 1:length(obj.list)) {
  #obj.list[[ii]] <- NormalizeData(obj.list[[ii]], verbose = FALSE)
  #obj.list[[ii]] <- FindVariableFeatures(obj.list[[ii]], selection.method = "vst",nfeatures = 3000, verbose = FALSE)
  obj.list[[ii]] <- SCTransform(obj.list[[ii]], verbose = FALSE)
}
#obj.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
#obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = obj.features, 
#                                    verbose = FALSE)




obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)



library(ggplot2)
library(cowplot)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
#for (dimval in 3:30) {
dimval=30
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8,dims=1:dimval)
obj.integrated = FindClusters(obj.integrated,resolution = 0.4)
DimPlot(obj.integrated,reduction="umap",label=T,group.by="cluster_label")
DimPlot(obj.integrated,reduction="umap")
#FeaturePlot(obj.integrated,features=c("FOXP2","CRH","EBF1"))
#table(obj.integrated@active.ident,obj.integrated@meta.data$cluster_label)
p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label")
pdf("human_macaque_magpar_integration.pdf")
print(p1)
dev.off()
save(obj.integrated,file="human_macaque_magpar_integrated.rda")

####plot FOXP2 expression###
coords=obj.integrated@reductions$umap@cell.embeddings
int_keepmac=intersect(rownames(coords)[coords[,1]< -5],allcl_mac$sample_name[allcl_mac$cluster_label %in% c("MP1_M","MP2_M")])
int_keephum=intersect(rownames(coords)[coords[,1]< -5],allcl_hum$sample_name[allcl_hum$cluster_label %in% c("MP1_H")])
pdf("human_macaque_MP_plots.pdf",useDingbats=F,width=24,height=12)
par(mfrow=c(1,2))
plot(coords[int_keepmac,1:2],
     col=allcl_mac$cluster_color[match(int_keepmac,allcl_mac$sample_name)],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
plot(coords[int_keephum,1:2],
     col=allcl_hum$cluster_color[match(int_keephum,allcl_hum$sample_name)],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
for (gene in c("FOXP2","SUSD2","EPHA1","GRID2","EBF1")) {
#gene="FOXP2"
genemac=datlist[["macaque"]][gene,int_keepmac]
genehum=datlist[["human"]][gene,int_keephum]
colpalette=colorRampPalette(c("black","red","yellow"))
colmac1=ceiling(254*log10(genemac+1)/max(log10(genemac+1)))+1
colhum1=ceiling(254*log10(genehum+1)/max(log10(genehum+1)))+1
colmac2=colpalette(256)[colmac1]
colhum2=colpalette(256)[colhum1]
par(mfrow=c(1,2))
plot(coords[int_keepmac[order(colmac1)],1:2],
     col=c(colmac2[order(colmac1)]),pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",main=gene)
plot(coords[int_keephum[order(colhum1)],1:2],
     col=c(colhum2[order(colhum1)]),pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
}
dev.off()

####use pca###
keepnam=intersect(allcl_hum$sample_name,rownames(obj.integrated@reductions$umap@cell.embeddings))
ndim=30
meandist=matrix(0,nrow=length(keepnam),ncol=2)
rownames(meandist)=keepnam
indset1=intersect(allcl_mac$sample_name[allcl_mac$cluster_label=="MP1_M"],rownames(obj.integrated@reductions$umap@cell.embeddings))
indset2=intersect(allcl_mac$sample_name[allcl_mac$cluster_label=="MP2_M"],rownames(obj.integrated@reductions$umap@cell.embeddings))
distmat=as.matrix(dist(obj.integrated@reductions$pca@cell.embeddings[,1:ndim]),"euclidean")
for (ii in keepnam) {
  meandist[ii,1]=mean(distmat[ii,indset1])
  meandist[ii,2]=mean(distmat[ii,indset2])
}
boxplot(datlist[["human"]]["FOXP2",rownames(meandist)]~as.factor(apply(meandist,1,which.max)))
colvec2=ifelse(meandist[,1]>meandist[,2],"red","blue")
#obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8,dims=1:ndim)


keepnam2=intersect(rownames(obj.integrated@reductions$umap@cell.embeddings),allcl_mac$sample_name)
colvec=datlist[["macaque"]]["FOXP2",keepnam2]
colvec=ceiling(colvec/max(colvec)*255)+1
scalevec=colorRampPalette(c("black","red","yellow"))
colvec2=scalevec(256)[colvec]
plot(obj.integrated@reductions$umap@cell.embeddings[keepnam2[order(colvec)],1:2],col=colvec2[order(colvec)],pch=19)
#     col=allcl_mac$cluster_color[match(keepnam2,allcl_mac$sample_name)],pch=19)


###ideas: use DE genes from macaque clusters (versus random genes) ###
require(edgeR)
mp1_nam=allcl_mac$sample_name[allcl_mac$cluster_label=="MP1_M"]
mp2_nam=allcl_mac$sample_name[allcl_mac$cluster_label=="MP2_M"]
classvec=as.factor(rep(c(1,2),times=c(length(mp1_nam),length(mp2_nam))))
startmat_cpm=datlist[["macaque"]][,c(mp1_nam,mp2_nam)]
e_design=model.matrix(~classvec)
y2 = DGEList(counts=datlist[["macaque"]][,c(mp1_nam,mp2_nam)])
y2 = estimateDisp(y2, e_design)
fit = glmQLFit(y2, e_design)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
outval2=topTags(qlf.2vs1,n=nrow(startmat_cpm),p.value=1)
mean1=apply(startmat_cpm[rownames(outval2$table),classvec==1],1,mean)
mean2=apply(startmat_cpm[rownames(outval2$table),classvec==2],1,mean)
frac1=apply(startmat_cpm[rownames(outval2$table),classvec==1]>0,1,mean)
frac2=apply(startmat_cpm[rownames(outval2$table),classvec==2]>0,1,mean)
outval2$table=cbind(outval2$table,mean1,mean2,frac1,frac2)


###heatmap of macaque genes###
keepgendiff=outval2$table[abs(outval2$table$logFC)>1 & outval2$table$FDR<0.01 & (outval2$table$frac1>0.3 | outval2$table$frac2>0.3),]
keepgendiff1=rownames(keepgendiff)[order(keepgendiff$logFC)[1:30]]
keepgendiff2=rownames(keepgendiff)[order(-keepgendiff$logFC)[1:30]]
require(gplots)
expmat=datlist$macaque[c(keepgendiff1,keepgendiff2),c(mp1_nam,mp2_nam)]
breakvals=seq(from=-2,to=2,length.out=256)
colnames(expmat)=NULL
pdf("macaque_mp_heatmap.pdf")
heatmap.2(log10(as.matrix(expmat+1)),scale='row',trace='none',breaks=breakvals,col=bluered(255),ColSideColors = rep(c("blue","red"),times=c(length(mp1_nam),length(mp2_nam))),Rowv=F,Colv=F,dendrogram='none')
dev.off()

mag_nam=intersect(allcl_hum$sample_name[allcl_hum$cluster_label=="MP1_H"],rownames(metalist[["human"]])[grep("MC",metalist$human$roi_label)])
par_nam=intersect(allcl_hum$sample_name[allcl_hum$cluster_label=="MP1_H"],rownames(metalist[["human"]])[grep("PC",metalist$human$roi_label)])
classvec=as.factor(rep(c(1,2),times=c(length(mag_nam),length(par_nam))))
startmat_cpm=datlist[["human"]][,c(mag_nam,par_nam)]
e_design=model.matrix(~classvec)
y2 = DGEList(counts=datlist[["human"]][,c(mag_nam,par_nam)])
y2 = estimateDisp(y2, e_design)
fit = glmQLFit(y2, e_design)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
outval2_mp=topTags(qlf.2vs1,n=nrow(startmat_cpm),p.value=1)
mean1=apply(startmat_cpm[rownames(outval2_mp$table),classvec==1],1,mean)
mean2=apply(startmat_cpm[rownames(outval2_mp$table),classvec==2],1,mean)
frac1=apply(startmat_cpm[rownames(outval2_mp$table),classvec==1]>0,1,mean)
frac2=apply(startmat_cpm[rownames(outval2_mp$table),classvec==2]>0,1,mean)
outval2_mp$table=cbind(outval2_mp$table,mean1,mean2,frac1,frac2)

save(outval2,outval2_mp,file="degenes_magpar.rda")


####cluster###
macgenes=rownames(outval2$table)[outval2$table$FDR<0.05 & abs(outval2$table$logFC)>0]
remcells=c("SM-D9CXJ_S07_E1-50","SM-D9CXL_S88_E1-50","SM-D9E69_S19_E1-50","SM-D9E6E_S56_E1-50","SM-D9EOQ_S48_E1-50","SM-D9EP3_S66_E1-50","SM-D9EPA_S72_E1-50",
           "SM-D9EPI_S57_E1-50","SM-D9EPU_S80_E1-50","SM-D9EPU_S83_E1-50")
keepcells1=setdiff(allcl_hum$sample_name[allcl_hum$cluster_label=="MP1_H"],remcells)
mp.int <- CreateSeuratObject(counts = datlist[["human"]][,keepcells1], meta.data = metalist[["human"]][keepcells1,])
mp.int <- NormalizeData(object = mp.int)
mp.int <- ScaleData(object = mp.int, verbose = FALSE)
#mp.int <- FindVariableFeatures(object = mp.int, selection.method = "vst", 
#                              top.genes = 2000)
mp.int <- RunPCA(object = mp.int, npcs = 30, verbose = FALSE,features=intersect(macgenes,rownames(datlist[["human"]])))
set.seed(0)
pdf("human_mp1.pdf")
for (dimval in 2:30) {
mp.int <- RunUMAP(object = mp.int, reduction = "pca", dims = 1:dimval)
mp.int <- FindNeighbors(object = mp.int, k.param = 8,dims=1:dimval)
mp.int <- FindClusters(object = mp.int, resolution = 0.4) ##0.2 for macaque
p1=DimPlot(mp.int,reduction = "umap",label = TRUE)
print(p1)
#p1=DimPlot(mp.int,reduction = "umap",group.by="roi_label",label = TRUE)
p1=DotPlot(mp.int,features=c("BTNL9","FOXP2","CRH","SUSD2","EPHA7"))
print(p1)
}
dev.off()
table(mp.int@active.ident,mp.int@meta.data$roi_label)








suballcl=allcl_hum[rownames(meandist),]
suballcl$cluster_id[meandist[,1]>meandist[,2]]=2
suballcl$cluster_id[meandist[,1]<meandist[,2]]=1
suballcl$cluster_color[suballcl$cluster_id==2]="#00CCDD"
sampdat=data.frame(sample_name=suballcl$sample_name,t(as.matrix(datlist[["human"]][,suballcl$sample_name])))
group_violin_plot(sampdat,suballcl,grouping="cluster",
                 genes=c("BTNL9","FOXP2","CRH","EBF1","ADGRG6","SERPINI2","KCNB2"),
                 log_scale=F)
fracdiff=data.frame(p1=apply(datlist[["human"]][,suballcl$sample_name[suballcl$cluster_id==1]]>0,1,mean),
                    p2=apply(datlist[["human"]][,suballcl$sample_name[suballcl$cluster_id==2]]>0,1,mean))



#p1 = DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label", label = TRUE, 
#              repel = TRUE)
#coords=obj.integrated@reductions$umap@cell.embeddings
#macplot=coords[intersect(rownames(coords),allcl_mac$sample_name),]
#humplot=coords[intersect(rownames(coords),allcl_hum$sample_name),]
#xvals=range(coords[,1])
#yvals=range(coords[,2])
#humcl2=rownames(humplot)[humplot[,1]>-2.3 & humplot[,2]<0]
#humcl1=setdiff(rownames(humplot),humcl2)
#humcl1=intersect(allcl_hum$sample_name,names(obj.integrated@active.ident)[obj.integrated@active.ident!=2])
#humcl2=intersect(allcl_hum$sample_name,names(obj.integrated@active.ident)[obj.integrated@active.ident==2])
#plot(macplot[,1:2],col=allcl_mac$cluster_color[match(rownames(macplot),allcl_mac$sample_name)],pch=19,xlim=xvals,ylim=yvals,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
#plot(humplot[,1:2],col=allcl_hum$cluster_color[match(rownames(humplot)[1],allcl_hum$sample_name)],pch=ifelse(rownames(humplot) %in% humcl1,19,1),xlim=xvals,ylim=yvals,xlab="UMAP Axis 1",ylab="UMAP Axis 2")

# suballcl=allcl_hum[c(humcl1,humcl2),]
# suballcl$cluster_id[match(humcl2,suballcl$sample_name)]=2
# suballcl$cluster_id[match(humcl1,suballcl$sample_name)]=1
# suballcl$cluster_color[match(humcl2,suballcl$sample_name)]="#00CCDD"
# sampdat=data.frame(sample_name=suballcl$sample_name,t(as.matrix(datlist[["human"]][,suballcl$sample_name])))
# group_violin_plot(sampdat,suballcl,grouping="cluster",
#                  genes=c("SNAP25","GAD2","SLC17A6","LAMP5","NXPH2","TRPC4","ALCAM","KRT80","PVALB","HGF",
#                          "PRKCG","CALB1","CASQ2","GRB14","PENK","BCHE",
#                          "BTNL9","FOXP2","CRH","EBF1"),
#                  log_scale=F)
dev.off()
ttt=FindMarkers(obj.integrated,ident.1=2)

####heatmap plots along UMAP axes####









######integrate mouse, macaque,and human
require(Seurat)
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


load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl_20200126.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")
load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_label=paste0(allcl_mus$cluster_label,"_Ms")


####integrate LGD/LP/K/PM populations####
#grepchar="MP"
grepchar="MP|K|LGD|LP"
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
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



library(ggplot2)
library(cowplot)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
pdf("human_macaque_mouse_magparkonlgdlp_integration_20200126.pdf")
for (dimval in 5:30) {
#dimval=30
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8,dims=1:dimval)
obj.integrated = FindClusters(obj.integrated,resolution = 0.4)
#DimPlot(obj.integrated,reduction="umap",label=T,group.by="cluster_label")
#DimPlot(obj.integrated,reduction="umap")
#FeaturePlot(obj.integrated,features=c("FOXP2","CRH","EBF1"))
#table(obj.integrated@active.ident,obj.integrated@meta.data$cluster_label)
p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T)
print(p1)
}
dev.off()
save(obj.integrated,file="human_macaque_mouse_magparkonlgdlp_integrated_20200126.rda")


###MP/LGD
grepchar="MP|LGD"
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
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



library(ggplot2)
library(cowplot)
DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
pdf("human_macaque_mouse_magparlgd_integration_20200126.pdf")
for (dimval in 5:30) {
  #dimval=30
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
  obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8,dims=1:dimval)
  obj.integrated = FindClusters(obj.integrated,resolution = 0.4)
  #DimPlot(obj.integrated,reduction="umap",label=T,group.by="cluster_label")
  #DimPlot(obj.integrated,reduction="umap")
  #FeaturePlot(obj.integrated,features=c("FOXP2","CRH","EBF1"))
  #table(obj.integrated@active.ident,obj.integrated@meta.data$cluster_label)
  p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T)
  print(p1)
}
dev.off()
save(obj.integrated,file="human_macaque_mouse_magparlgdlp_integrated_20200126.rda")


####integrate all cells####
#grepchar="MP"
grepchar=""
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
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
  obj.list[[ii]] <- SCTransform(obj.list[[ii]], verbose = TRUE)
}
pdf("human_macaque_mouse_allcells_integration_20200126.pdf")
for (anchorfeat in c(2000,3000,5000,8000)) {
  gc()
  obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30,anchor.features=anchorfeat)
  obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)
  DefaultAssay(obj.integrated) <- "integrated"
  obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
  obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
  for (dimval in 5:30) {
    #dimval=30
    obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
    obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8,dims=1:dimval)
    obj.integrated = FindClusters(obj.integrated,resolution = 0.4)
    #DimPlot(obj.integrated,reduction="umap",label=T,group.by="cluster_label")
    #DimPlot(obj.integrated,reduction="umap")
    #FeaturePlot(obj.integrated,features=c("FOXP2","CRH","EBF1"))
    #table(obj.integrated@active.ident,obj.integrated@meta.data$cluster_label)
    p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T)
    print(p1)
  }
}
dev.off()
save(obj.integrated,file="human_macaque_mouse_allcells_integrated_20200126.rda")

####integrate GABAergic cells####
#grepchar="MP"
grepchar="GABA"
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor %in% c("Q17","Q18")))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
datlist2=datlist
for (ii in 1:length(datlist2)) {
  rownames(datlist2[[ii]])=tolower(rownames(datlist2[[ii]]))
}
obj.list = list()
obj.list[["macaque1"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac1]),meta.data=allcl_mac[keepcells_mac1,])
#obj.list[["macaque2"]]=CreateSeuratObject(round(datlist2[["macaque"]][,keepcells_mac2]),meta.data=allcl_mac[keepcells_mac2,])
obj.list[["human"]]=CreateSeuratObject(round(datlist2[["human"]][,keepcells_hum]),meta.data=allcl_hum[keepcells_hum,])
obj.list[["mouse"]]=CreateSeuratObject(round(datlist2[["mouse"]][,keepcells_mus]),meta.data=allcl_mus[keepcells_mus,])
for (ii in 1:length(obj.list)) {
  obj.list[[ii]] <- SCTransform(obj.list[[ii]], verbose = TRUE)
}
pdf("human_macaque_mouse_GABAergic_integration_20200126.pdf")
for (anchorfeat in c(2000,3000,5000,8000)) {
  gc()
  obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30,anchor.features=anchorfeat,k.filter=40)
  obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)
  DefaultAssay(obj.integrated) <- "integrated"
  obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
  obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = TRUE)
  for (dimval in 5:30) {
    #dimval=30
    obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
    obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8,dims=1:dimval)
    obj.integrated = FindClusters(obj.integrated,resolution = 0.4)
    #DimPlot(obj.integrated,reduction="umap",label=T,group.by="cluster_label")
    #DimPlot(obj.integrated,reduction="umap")
    #FeaturePlot(obj.integrated,features=c("FOXP2","CRH","EBF1"))
    #table(obj.integrated@active.ident,obj.integrated@meta.data$cluster_label)
    p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T)
    print(p1)
  }
}
dev.off()
save(obj.integrated,file="human_macaque_mouse_GABAcells_integrated_20200126.rda")




###make integrated plots####
###plot 1: all cells, all species###
load("human_macaque_mouse_allcells_integrated.rda")
load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")
load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_label=paste0(allcl_mus$cluster_label,"_Ms")
dimval=8
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
grepchar=""
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
coords=obj.integrated@reductions$umap@cell.embeddings[,1:2]
colclust=rep("grey",nrow(coords))
names(colclust)=rownames(coords)
mac_cells=intersect(rownames(coords),c(keepcells_mac1,keepcells_mac2))
colclust[mac_cells]=allcl_mac$cluster_color[match(mac_cells,allcl_mac$sample_name)]
hum_cells=intersect(rownames(coords),c(keepcells_hum))
colclust[hum_cells]=allcl_hum$cluster_color[match(hum_cells,allcl_hum$sample_name)]
mouse_cells=intersect(rownames(coords),c(keepcells_mus))
colclust[mouse_cells]=allcl_mus$cluster_color[match(mouse_cells,allcl_mus$sample_name)]
table(colclust)
plot(coords[,1:2],col=colclust,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")

colreg=rep("grey",nrow(coords))
names(colreg)=rownames(coords)
colreg[mac_cells]=allcl_mac$roi[match(mac_cells,allcl_mac$sample_name)]
hum_cells=intersect(rownames(coords),c(keepcells_hum))
colreg[hum_cells]=as.character(metalist$human$roi_label[match(hum_cells,metalist$human$sample_id)])
mouse_cells=intersect(rownames(coords),c(keepcells_mus))
colreg[mouse_cells]=allcl_mus$roi[match(mouse_cells,allcl_mus$sample_name)]
table(colreg)
colreg[grep("KC",colreg)]="#00B914"
colreg[grep("MC",colreg)]="#0000FF"
colreg[grep("PC",colreg)]="#FF0000"
colreg[grep("Magnocellular",colreg)]="#0000FF"
colreg[grep("Parvocellular",colreg)]="#FF0000"
colreg[grep("Pulvinar",colreg)]="#BEBEBE"
colreg[grep("Core",colreg)]="#580077"
colreg[grep("Shell",colreg)]="#D066E0"
colreg[grep("LGv",colreg)]="#8FBE3F"
colreg[grep("LP",colreg)]="#1C75BC"
table(colreg)
plot(coords[,1:2],col=colreg,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")

####by species###
xlimval=range(coords[,1])
ylimval=range(coords[,2])
pdf("human_macaque_mouse_allcells_integration_panels.pdf",height=25,width=12,useDingbats=F)
par(mfrow=c(4,2))
plot(coords[,1:2],col=colclust,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[,1:2],col=colreg,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mac_cells,1:2],col=colclust[mac_cells],pch=19,xlab="UMAP AXis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mac_cells,1:2],col=colreg[mac_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[hum_cells,1:2],col=colclust[hum_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[hum_cells,1:2],col=colreg[hum_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mouse_cells,1:2],col=colclust[mouse_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mouse_cells,1:2],col=colreg[mouse_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
dev.off()


###plot 1: only LGD/MPK, all species###
load("human_macaque_mouse_magparlgd_integrated.rda")
load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")
load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_label=paste0(allcl_mus$cluster_label,"_Ms")
dimval=22
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
grepchar=""
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
coords=obj.integrated@reductions$umap@cell.embeddings[,1:2]
colclust=rep("grey",nrow(coords))
names(colclust)=rownames(coords)
mac_cells=intersect(rownames(coords),c(keepcells_mac1,keepcells_mac2))
colclust[mac_cells]=allcl_mac$cluster_color[match(mac_cells,allcl_mac$sample_name)]
hum_cells=intersect(rownames(coords),c(keepcells_hum))
colclust[hum_cells]=allcl_hum$cluster_color[match(hum_cells,allcl_hum$sample_name)]
mouse_cells=intersect(rownames(coords),c(keepcells_mus))
colclust[mouse_cells]=allcl_mus$cluster_color[match(mouse_cells,allcl_mus$sample_name)]
table(colclust)
plot(coords[,1:2],col=colclust,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")

colreg=rep("grey",nrow(coords))
names(colreg)=rownames(coords)
colreg[mac_cells]=allcl_mac$roi[match(mac_cells,allcl_mac$sample_name)]
hum_cells=intersect(rownames(coords),c(keepcells_hum))
colreg[hum_cells]=as.character(metalist$human$roi_label[match(hum_cells,metalist$human$sample_id)])
mouse_cells=intersect(rownames(coords),c(keepcells_mus))
colreg[mouse_cells]=allcl_mus$roi[match(mouse_cells,allcl_mus$sample_name)]
table(colreg)
colreg[grep("KC",colreg)]="#00B914"
colreg[grep("MC",colreg)]="#0000FF"
colreg[grep("PC",colreg)]="#FF0000"
colreg[grep("Magnocellular",colreg)]="#0000FF"
colreg[grep("Parvocellular",colreg)]="#FF0000"
colreg[grep("Pulvinar",colreg)]="#BEBEBE"
colreg[grep("Core",colreg)]="#580077"
colreg[grep("Shell",colreg)]="#D066E0"
colreg[grep("LGv",colreg)]="#8FBE3F"
colreg[grep("LP",colreg)]="#1C75BC"
table(colreg)
plot(coords[,1:2],col=colreg,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")

####by species###
xlimval=range(coords[,1])
ylimval=range(coords[,2])
pdf("human_macaque_mouse_magparlgd_integration_panels.pdf",height=25,width=12,useDingbats=F)
par(mfrow=c(4,2))
plot(coords[,1:2],col=colclust,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[,1:2],col=colreg,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mac_cells,1:2],col=colclust[mac_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mac_cells,1:2],col=colreg[mac_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[hum_cells,1:2],col=colclust[hum_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[hum_cells,1:2],col=colreg[hum_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mouse_cells,1:2],col=colclust[mouse_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mouse_cells,1:2],col=colreg[mouse_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
dev.off()


###plot 3: only LGD+LP/MPK, all species###
load("human_macaque_mouse_magparlgdlp_integrated.rda")
load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")
load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_label=paste0(allcl_mus$cluster_label,"_Ms")
dimval=22
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
grepchar=""
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
coords=obj.integrated@reductions$umap@cell.embeddings[,1:2]
colclust=rep("grey",nrow(coords))
names(colclust)=rownames(coords)
mac_cells=intersect(rownames(coords),c(keepcells_mac1,keepcells_mac2))
colclust[mac_cells]=allcl_mac$cluster_color[match(mac_cells,allcl_mac$sample_name)]
hum_cells=intersect(rownames(coords),c(keepcells_hum))
colclust[hum_cells]=allcl_hum$cluster_color[match(hum_cells,allcl_hum$sample_name)]
mouse_cells=intersect(rownames(coords),c(keepcells_mus))
colclust[mouse_cells]=allcl_mus$cluster_color[match(mouse_cells,allcl_mus$sample_name)]
table(colclust)
plot(coords[,1:2],col=colclust,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")

colreg=rep("grey",nrow(coords))
names(colreg)=rownames(coords)
colreg[mac_cells]=allcl_mac$roi[match(mac_cells,allcl_mac$sample_name)]
hum_cells=intersect(rownames(coords),c(keepcells_hum))
colreg[hum_cells]=as.character(metalist$human$roi_label[match(hum_cells,metalist$human$sample_id)])
mouse_cells=intersect(rownames(coords),c(keepcells_mus))
colreg[mouse_cells]=allcl_mus$roi[match(mouse_cells,allcl_mus$sample_name)]
table(colreg)
colreg[grep("KC",colreg)]="#00B914"
colreg[grep("MC",colreg)]="#0000FF"
colreg[grep("PC",colreg)]="#FF0000"
colreg[grep("Magnocellular",colreg)]="#0000FF"
colreg[grep("Parvocellular",colreg)]="#FF0000"
colreg[grep("Pulvinar",colreg)]="#BEBEBE"
colreg[grep("Core",colreg)]="#580077"
colreg[grep("Shell",colreg)]="#D066E0"
colreg[grep("LGv",colreg)]="#8FBE3F"
colreg[grep("LP",colreg)]="#1C75BC"
table(colreg)
plot(coords[,1:2],col=colreg,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")

####by species###
xlimval=range(coords[,1])
ylimval=range(coords[,2])
pdf("human_macaque_mouse_magparlgdlp_integration_panels.pdf",height=25,width=12,useDingbats=F)
par(mfrow=c(4,2))
plot(coords[,1:2],col=colclust,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[,1:2],col=colreg,pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mac_cells,1:2],col=colclust[mac_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mac_cells,1:2],col=colreg[mac_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[hum_cells,1:2],col=colclust[hum_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[hum_cells,1:2],col=colreg[hum_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mouse_cells,1:2],col=colclust[mouse_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
plot(coords[mouse_cells,1:2],col=colreg[mouse_cells],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2",xlim=xlimval,ylim=ylimval,cex=0.5)
dev.off()


###plot 3 - use UMAP axes in each direction###
load("collected_data_20191220.rda")
spec="macaque"
orthos=read.csv("ortholog_table_20191122.csv",as.is=T)
nonlocgenes=orthos$rhesus_symbol[intersect(grep("^LOC",orthos$human_symbol,invert=T),grep("^LOC",orthos$rhesus_symbol))]
keepgen=setdiff(rownames(datlist[[spec]]),setdiff(grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|^TRNA",rownames(datlist[[spec]]),val=T),nonlocgenes))
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
spec="human"
keepgen=grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|[0-9]P$",rownames(datlist[[spec]]),val=T,invert=T)
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
spec="mouse"
keepgen=grep("^LOC|^mt-|^Rpl|^Rps",rownames(datlist[[spec]]),val=T,invert=T)
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6

axis_genelists=list()
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP",allcl_mac$cluster_label)] & coords[,2]>-2)]
mac_cells_axis1=coords[keepcells,1]
genvals=t(cor(mac_cells_axis1,as.matrix(t(datlist[["macaque"]][,keepcells])),method="spearman"))
axis_genelists[["macaque_axis1"]]=rownames(genvals)[c(order(-genvals)[1:20],order(genvals)[1:20])]

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP|K",allcl_mac$cluster_label)])]
mac_cells_axis1=coords[keepcells,2]
genvals=t(cor(mac_cells_axis1,as.matrix(t(datlist[["macaque"]][,keepcells])),method="spearman"))
axis_genelists[["macaque_axis2"]]=rownames(genvals)[c(order(-genvals)[1:20],order(genvals)[1:20])]

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP",allcl_hum$cluster_label)] & coords[,2]>-2)]
hum_cells_axis1=coords[keepcells,1]
genvals=t(cor(hum_cells_axis1,as.matrix(t(datlist[["human"]][,keepcells])),method="spearman"))
axis_genelists[["human_axis1"]]=rownames(genvals)[c(order(-genvals)[1:20],order(genvals)[1:20])]

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP|K",allcl_hum$cluster_label)])]
hum_cells_axis1=coords[keepcells,2]
genvals=t(cor(hum_cells_axis1,as.matrix(t(datlist[["human"]][,keepcells])),method="spearman"))
axis_genelists[["human_axis2"]]=rownames(genvals)[c(order(-genvals)[1:20],order(genvals)[1:20])]

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)] & coords[,2]>-2)]
mus_cells_axis1=coords[keepcells,1]
genvals=t(cor(mus_cells_axis1,as.matrix(t(datlist[["mouse"]][,keepcells])),method="spearman"))
axis_genelists[["mouse_axis1"]]=rownames(genvals)[c(order(-genvals)[1:20],order(genvals)[1:20])]

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)])]
mus_cells_axis1=coords[keepcells,2]
genvals=t(cor(mus_cells_axis1,as.matrix(t(datlist[["mouse"]][,keepcells])),method="spearman"))
axis_genelists[["mouse_axis2"]]=rownames(genvals)[c(order(-genvals)[1:20],order(genvals)[1:20])]

###make axis heatmaps###
require(gplots)
colreg=allcl_mac$roi
colreg[grep("Magnocellular",colreg)]="#0000FF"
colreg[grep("Parvocellular",colreg)]="#FF0000"
colreg[grep("Pulvinar",colreg)]="#BEBEBE"
allcl_mac$roi_color=colreg

colreg=as.character(metalist$human$roi_label[match(allcl_hum$sample_name,metalist$human$sample_id)])
colreg[grep("KC",colreg)]="#00B914"
colreg[grep("MC",colreg)]="#0000FF"
colreg[grep("PC",colreg)]="#FF0000"
allcl_hum$roi_color=colreg

colreg=allcl_mus$roi
colreg[grep("Core",colreg)]="#580077"
colreg[grep("Shell",colreg)]="#D066E0"
colreg[grep("LGv",colreg)]="#8FBE3F"
colreg[grep("LP",colreg)]="#1C75BC"
allcl_mus$roi_color=colreg

pdf("gene_heatmaps.pdf",useDingbats=F)
##axis 1 genes
#allaxis1=c(axis_genelists$macaque_axis1,axis_genelists$human_axis1,axis_genelists$mouse_axis1)
#allaxis1=allaxis1[which(!duplicated(tolower(allaxis1)))]
##macaque, axis 1
mac_genes=axis_genelists$macaque_axis1
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP",allcl_mac$cluster_label)] & coords[,2]>-2)]
mac_cells_axis1=coords[keepcells,1]
plotmat=datlist$macaque[match(tolower(mac_genes),tolower(rownames(datlist$macaque))),keepcells[order(mac_cells_axis1)]]
plotmat=plotmat[apply(plotmat,1,sd)>0,]
corval=cor(1:ncol(plotmat),t(as.matrix(plotmat)))
plotmat=plotmat[order(-corval),]
breakvals=seq(from=-2,to=2,length.out=256)
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_mac$cluster_color[match(colnames(plotmat),allcl_mac$sample_name)])
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_mac$roi_color[match(colnames(plotmat),allcl_mac$sample_name)])


hum_genes=axis_genelists$human_axis1
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP",allcl_hum$cluster_label)] & coords[,2]>-2)]
hum_cells_axis1=coords[keepcells,1]
plotmat=datlist$human[match(tolower(hum_genes),tolower(rownames(datlist$human))),keepcells[order(hum_cells_axis1)]]
plotmat=plotmat[apply(plotmat,1,sd)>0,]
corval=cor(1:ncol(plotmat),t(as.matrix(plotmat)))
plotmat=plotmat[order(-corval),]
breakvals=seq(from=-2,to=2,length.out=256)
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_hum$roi_color[match(colnames(plotmat),allcl_hum$sample_name)])

mus_genes=axis_genelists$mouse_axis1
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)] & coords[,2]>-2)]
mus_cells_axis1=coords[keepcells,1]
plotmat=datlist$mouse[match(tolower(mus_genes),tolower(rownames(datlist$mouse))),keepcells[order(mus_cells_axis1)]]
plotmat=plotmat[apply(plotmat,1,sd)>0,]
corval=cor(1:ncol(plotmat),t(as.matrix(plotmat)))
plotmat=plotmat[order(-corval),]
breakvals=seq(from=-2,to=2,length.out=256)
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_mus$roi_color[match(colnames(plotmat),allcl_mus$sample_name)])

mac_genes=axis_genelists$macaque_axis2
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP|K",allcl_mac$cluster_label)])]
mac_cells_axis1=coords[keepcells,2]
plotmat=datlist$macaque[match(tolower(mac_genes),tolower(rownames(datlist$macaque))),keepcells[order(mac_cells_axis1)]]
plotmat=plotmat[apply(plotmat,1,sd)>0,]
corval=cor(1:ncol(plotmat),t(as.matrix(plotmat)))
plotmat=plotmat[order(-corval),]
breakvals=seq(from=-2,to=2,length.out=256)
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_mac$cluster_color[match(colnames(plotmat),allcl_mac$sample_name)])
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_mac$roi_color[match(colnames(plotmat),allcl_mac$sample_name)])

hum_genes=axis_genelists$human_axis2
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP|K",allcl_hum$cluster_label)])]
hum_cells_axis1=coords[keepcells,2]
plotmat=datlist$human[match(tolower(hum_genes),tolower(rownames(datlist$human))),keepcells[order(hum_cells_axis1)]]
plotmat=plotmat[apply(plotmat,1,sd)>0,]
corval=cor(1:ncol(plotmat),t(as.matrix(plotmat)))
plotmat=plotmat[order(-corval),]
breakvals=seq(from=-2,to=2,length.out=256)
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_hum$cluster_color[match(colnames(plotmat),allcl_hum$sample_name)])
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_hum$roi_color[match(colnames(plotmat),allcl_hum$sample_name)])

mus_genes=axis_genelists$mouse_axis2
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)])]
mus_cells_axis1=coords[keepcells,2]
plotmat=datlist$mouse[match(tolower(mus_genes),tolower(rownames(datlist$mouse))),keepcells[order(mus_cells_axis1)]]
plotmat=plotmat[apply(plotmat,1,sd)>0,]
corval=cor(1:ncol(plotmat),t(as.matrix(plotmat)))
plotmat=plotmat[order(-corval),]
breakvals=seq(from=-2,to=2,length.out=256)
heatmap.2(log10(as.matrix(plotmat)+1),scale="row",trace="none",Colv=F,Rowv=F,labCol='',col=bluered,breaks=breakvals,ColSideColors=allcl_mus$roi_color[match(colnames(plotmat),allcl_mus$sample_name)])
dev.off()


####correlation with axes plots###
load("collected_data_20191220.rda")
spec="macaque"
orthos=read.csv("ortholog_table_20191122.csv",as.is=T)
nonlocgenes=orthos$rhesus_symbol[intersect(grep("^LOC",orthos$human_symbol,invert=T),grep("^LOC",orthos$rhesus_symbol))]
keepgen=setdiff(rownames(datlist[[spec]]),setdiff(grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|^TRNA",rownames(datlist[[spec]]),val=T),nonlocgenes))
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
spec="human"
keepgen=grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|[0-9]P$",rownames(datlist[[spec]]),val=T,invert=T)
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
spec="mouse"
keepgen=grep("^LOC|^mt-|^Rpl|^Rps",rownames(datlist[[spec]]),val=T,invert=T)
datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6

load("human_macaque_mouse_magparlgdlp_integrated.rda")
load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")
load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_label=paste0(allcl_mus$cluster_label,"_Ms")
dimval=22
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
grepchar=""
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
coords=obj.integrated@reductions$umap@cell.embeddings[,1:2]

require(Matrix)
axis_genelists=list()
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP",allcl_mac$cluster_label)] & coords[,2]>-2)]
mac_cells_axis1=coords[keepcells,1]
genvals=t(cor(mac_cells_axis1,as.matrix(t(datlist[["macaque"]][,keepcells])),method="spearman"))
axis_genelists[["macaque_axis1"]]=genvals

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP|K",allcl_mac$cluster_label)])]
mac_cells_axis1=coords[keepcells,2]
genvals=t(cor(mac_cells_axis1,as.matrix(t(datlist[["macaque"]][,keepcells])),method="spearman"))
axis_genelists[["macaque_axis2"]]=genvals

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP",allcl_hum$cluster_label)] & coords[,2]>-2)]
hum_cells_axis1=coords[keepcells,1]
genvals=t(cor(hum_cells_axis1,as.matrix(t(datlist[["human"]][,keepcells])),method="spearman"))
axis_genelists[["human_axis1"]]=genvals

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP|K",allcl_hum$cluster_label)])]
hum_cells_axis1=coords[keepcells,2]
genvals=t(cor(hum_cells_axis1,as.matrix(t(datlist[["human"]][,keepcells])),method="spearman"))
axis_genelists[["human_axis2"]]=genvals

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)] & coords[,2]>-2)]
mus_cells_axis1=coords[keepcells,1]
genvals=t(cor(mus_cells_axis1,as.matrix(t(datlist[["mouse"]][,keepcells])),method="spearman"))
axis_genelists[["mouse_axis1"]]=genvals

keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)])]
mus_cells_axis1=coords[keepcells,2]
genvals=t(cor(mus_cells_axis1,as.matrix(t(datlist[["mouse"]][,keepcells])),method="spearman"))
axis_genelists[["mouse_axis2"]]=genvals

require(wordcloud)
pdf("pairwise_genecors_axis.pdf",useDingbats=F)
keepgen=intersect(tolower(rownames(axis_genelists[["macaque_axis1"]])),tolower(rownames(axis_genelists[["human_axis1"]])))
g1=axis_genelists[["macaque_axis1"]][match(keepgen,tolower(rownames(axis_genelists[["macaque_axis1"]]))),1]
g2=axis_genelists[["human_axis1"]][match(keepgen,tolower(rownames(axis_genelists[["human_axis1"]]))),1]
#g1[is.na(g1)]=0;g2[is.na(g2)]=0;ord1=order(-g1);ord2=order(-g2)
#keepplots=intersect(ord1[1:200],ord2[1:200])
#keepplots2=intersect(rev(ord1)[1:200],rev(ord2)[1:200])
#keepplots=c(keepplots,keepplots2)
keepplots=which((g1>0.25 & g2>0.25) |(g1< -0.25 & g2< -0.25) )
rangevals=range(c(g1[keepplots],g2[keepplots]))
textplot(g1[keepplots],g2[keepplots],names(g1)[keepplots],cex=0.8,xlim=rangevals,ylim=rangevals,xlab="Correlation with UMAP Axis 1 (Macaque)",ylab="Correlation with UMAP Axis 1 (Human)",col="red")

keepgen=intersect(tolower(rownames(axis_genelists[["macaque_axis1"]])),tolower(rownames(axis_genelists[["mouse_axis1"]])))
g1=axis_genelists[["macaque_axis1"]][match(keepgen,tolower(rownames(axis_genelists[["macaque_axis1"]]))),1]
g2=axis_genelists[["mouse_axis1"]][match(keepgen,tolower(rownames(axis_genelists[["mouse_axis1"]]))),1]
#g1[is.na(g1)]=0;g2[is.na(g2)]=0;ord1=order(-g1);ord2=order(-g2)
#keepplots=intersect(ord1[1:200],ord2[1:200])
#keepplots2=intersect(rev(ord1)[1:200],rev(ord2)[1:200])
#keepplots=c(keepplots,keepplots2)
keepplots=which((g1>0.183 & g2>0.183) |(g1< -0.183 & g2< -0.183) )
rangevals=range(c(g1[keepplots],g2[keepplots]))
textplot(g1[keepplots],g2[keepplots],names(g1)[keepplots],cex=0.8,xlim=rangevals,ylim=rangevals,xlab="Correlation with UMAP Axis 1 (Macaque)",ylab="Correlation with UMAP Axis 1 (Mouse)",col="blue")

keepgen=intersect(tolower(rownames(axis_genelists[["human_axis1"]])),tolower(rownames(axis_genelists[["mouse_axis1"]])))
g1=axis_genelists[["human_axis1"]][match(keepgen,tolower(rownames(axis_genelists[["human_axis1"]]))),1]
g2=axis_genelists[["mouse_axis1"]][match(keepgen,tolower(rownames(axis_genelists[["mouse_axis1"]]))),1]
#g1[is.na(g1)]=0;g2[is.na(g2)]=0;ord1=order(-g1);ord2=order(-g2)
#keepplots=intersect(ord1[1:200],ord2[1:200])
#keepplots2=intersect(rev(ord1)[1:200],rev(ord2)[1:200])
#keepplots=c(keepplots,keepplots2)
keepplots=which((g1>0.183 & g2>0.183) |(g1< -0.183 & g2< -0.183) )
rangevals=range(c(g1[keepplots],g2[keepplots]))
textplot(g1[keepplots],g2[keepplots],names(g1)[keepplots],cex=0.8,xlim=rangevals,ylim=rangevals,xlab="Correlation with UMAP Axis 1 (Human)",ylab="Correlation with UMAP Axis 1 (Mouse)",col="purple")

keepgen=intersect(tolower(rownames(axis_genelists[["macaque_axis2"]])),tolower(rownames(axis_genelists[["human_axis2"]])))
g1=axis_genelists[["macaque_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["macaque_axis2"]]))),1]
g2=axis_genelists[["human_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["human_axis2"]]))),1]
keepplots=which((g1>0.25 & g2>0.25) |(g1< -0.25 & g2< -0.25) )
rangevals=range(c(g1[keepplots],g2[keepplots]))
textplot(g1[keepplots],g2[keepplots],names(g1)[keepplots],cex=0.8,xlim=rangevals,ylim=rangevals,xlab="Correlation with UMAP Axis 2 (Macaque)",ylab="Correlation with UMAP Axis 1 (Human)",col="red")

keepgen=intersect(tolower(rownames(axis_genelists[["macaque_axis2"]])),tolower(rownames(axis_genelists[["mouse_axis2"]])))
g1=axis_genelists[["macaque_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["macaque_axis2"]]))),1]
g2=axis_genelists[["mouse_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["mouse_axis2"]]))),1]
keepplots=which((g1>0.183 & g2>0.183) |(g1< -0.183 & g2< -0.183) )
rangevals=range(c(g1[keepplots],g2[keepplots]))
textplot(g1[keepplots],g2[keepplots],names(g1)[keepplots],cex=0.8,xlim=rangevals,ylim=rangevals,xlab="Correlation with UMAP Axis 2 (Macaque)",ylab="Correlation with UMAP Axis 2 (Mouse)",col="blue")

keepgen=intersect(tolower(rownames(axis_genelists[["human_axis2"]])),tolower(rownames(axis_genelists[["mouse_axis2"]])))
g1=axis_genelists[["human_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["human_axis2"]]))),1]
g2=axis_genelists[["mouse_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["mouse_axis2"]]))),1]
keepplots=which((g1>0.183 & g2>0.183) |(g1< -0.183 & g2< -0.183) )
rangevals=range(c(g1[keepplots],g2[keepplots]))
textplot(g1[keepplots],g2[keepplots],names(g1)[keepplots],cex=0.8,xlim=rangevals,ylim=rangevals,xlab="Correlation with UMAP Axis 2 (Human)",ylab="Correlation with UMAP Axis 2 (Mouse)",col="purple")
dev.off()


###triple plot###
keepgen=intersect(tolower(rownames(axis_genelists[["macaque_axis2"]])),intersect(tolower(rownames(axis_genelists[["human_axis2"]])),tolower(rownames(axis_genelists[["mouse_axis2"]]))))
g1=axis_genelists[["macaque_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["macaque_axis2"]]))),1]
g2=axis_genelists[["human_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["human_axis2"]]))),1]
g3=axis_genelists[["mouse_axis2"]][match(keepgen,tolower(rownames(axis_genelists[["mouse_axis2"]]))),1]
g1[is.na(g1)]=0;g2[is.na(g2)]=0;g3[is.na(g3)]=0
ord1=1:length(g1);ord1[order(g1)]=1:length(g1);
ord2=1:length(g2);ord2[order(g2)]=1:length(g2);
ord3=1:length(g3);ord3[order(g3)]=1:length(g3);
allord=cbind(ord1,ord2,ord3);rownames(allord)=names(g1)
allgen=cbind(g1,g2,g3);rownames(allgen)=names(g1)
minval=apply(allgen,1,min)
range(minval)
