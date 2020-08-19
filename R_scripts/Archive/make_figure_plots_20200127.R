library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)
source("doublet.finder.r")
require(Matrix)

load("collected_data_20191220.rda")

spec="macaque"

####make plots####
###violin plot###
source("hicat_modules.r")
require(scrattch.vis)
require(Seurat)
if (spec=="macaque") {
  orthos=read.csv("ortholog_table_20191122.csv",as.is=T)
  nonlocgenes=orthos$rhesus_symbol[intersect(grep("^LOC",orthos$human_symbol,invert=T),grep("^LOC",orthos$rhesus_symbol))]
  keepgen=setdiff(rownames(datlist[[spec]]),setdiff(grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|^TRNA",rownames(datlist[[spec]]),val=T),nonlocgenes))
  datlist[[spec]]=datlist[[spec]][keepgen,]
  datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
  load("macaque_clusters_subclusters.rda")
  allcl=data.frame(sample_name=colnames(m.int@assays$RNA),cluster_id=1,cluster_label=0,cluster_color="black",stringsAsFactors = F)
  rownames(allcl)=allcl$sample_name
  allcl[colnames(m.int6@assays$RNA),"cluster_label"]=paste0("GABA",m.int6@active.ident)
  allcl[colnames(m.int5@assays$RNA),"cluster_label"]=paste0("K",m.int5@active.ident)
  allcl[colnames(m.int4@assays$RNA),"cluster_label"]=paste0("MP",m.int4@active.ident)
  allcl[colnames(m.int3@assays$RNA)[m.int3@active.ident %in% c(2,5,7)],"cluster_label"]=paste0("O",m.int3@active.ident[m.int3@active.ident %in% c(2,5,7)])
  table(allcl$cluster_label)
  ##remove outlier clusters - low quality cells###
  allcl=allcl[allcl$cluster_label!="MP2",]
  ###merge clusters with sex-specific genes
  allcl$cluster_label=gsub("GABA0","GABA3",allcl$cluster_label)
  allcl$cluster_label=gsub("GABA5","GABA2",allcl$cluster_label)
  allcl$cluster_label=gsub("K2","K1",allcl$cluster_label)
  allcl$cluster_label=gsub("K0","K2",allcl$cluster_label)
  allcl$cluster_label=gsub("MP3","MP1",allcl$cluster_label)
  allcl$cluster_label=gsub("MP0","MP2",allcl$cluster_label)
  allcl$cluster_label=gsub("O2","O1",allcl$cluster_label)
  allcl$cluster_label=gsub("O5","O2",allcl$cluster_label)
  allcl$cluster_label=gsub("O7","O3",allcl$cluster_label)
  
  
  ###assign order & colors###
  allcl$cluster_id=match(allcl$cluster_label,c("GABA1","GABA2","GABA3","GABA4","K1","K2","MP1","MP2","O1","O2","O3")) 
  colvec=c("#6A3D9A","#F781BF","#E31A1C","#FF7F00","#33A02C","#00DD00","#377EB8","#00CCDD","#DDDDDD","#888888","#444444")
  allcl$cluster_color=colvec[allcl$cluster_id]
  
  #allcl$cluster_id=as.numeric(as.factor(allcl$cluster_label))
  #allcl$cluster_color=rainbow(length(unique(allcl$cluster_label)),start=0,end=4/6)[as.numeric(as.factor(allcl$cluster_label))]
  
  allcl$donor_sample=metalist[[spec]]$external_donor_name[match(allcl$sample_name,rownames(metalist[[spec]]))]
  allcl$donor=metalist[[spec]]$donor[match(allcl$sample_name,rownames(metalist[[spec]]))]
  allcl$roi=metalist[[spec]]$roi[match(allcl$sample_name,rownames(metalist[[spec]]))]
  save(allcl,file="macaque_allcl.rda")
  
  
  sampdat=data.frame(sample_name=allcl$sample_name,t(as.matrix(datlist[[spec]][,allcl$sample_name])))
  pdf("macaque_violin_plot.pdf",height=10,width=3)
  group_violin_plot(sampdat,allcl,grouping="cluster",
                    genes=c("SNAP25","GAD2","SLC17A6","DLX1","NPY","TNFAIP8L3","LAMP5","KRT80","DSEL","SLITRK2","COL2A1","NXPH2","PRKCG",
                            "G0S2","PENK","SGCD","DAB2","LVRN","CASQ2",
                            "IL15RA","FOXP2","EBF1","FGF11",
                            "GRIK3","PLPP4","RUFY4","PTPN20","GPR37","KCNJ2"),
                    log_scale=F)
  dev.off()
  
  
  table(paste(allcl$donor,allcl$roi),allcl$cluster_label)
  ####UMAP plots###
  ###new overall umap plot###
  plot_objects=list()
  subcells=list()
  subcells[["all"]]=1:nrow(allcl)
  subcells[["mp"]]=grep("^MP",allcl$cluster_label)
  subdat=datlist[[spec]][,allcl$sample_name]
  subdat=sweep(subdat,2,Matrix::colSums(subdat),"/")*10^6
  for (nam in names(subcells)) {
    keepcells=subcells[[nam]]
    m.overall <- CreateSeuratObject(counts = subdat[,allcl$sample_name[keepcells]], meta.data = allcl[keepcells,])
    m.overall <- NormalizeData(object =m.overall)
    m.overall <- ScaleData(object =m.overall, verbose = TRUE,vars.to.regress = "donor")
    m.overall <- FindVariableFeatures(object =m.overall, selection.method = "vst", top.genes = 2000)
    m.overall <- RunPCA(object =m.overall, npcs = 30, verbose = TRUE)
    m.overall <- RunTSNE(object =m.overall, reduction = "pca", dims = 1:15)
    m.overall <- RunUMAP(object =m.overall, reduction = "pca", dims = 1:15)
    plot_objects[[nam]]=m.overall
    print(nam)
  }
  #save(plot_objects,file="macaque_plotobjects.rda")
  #save(plot_objects,file="macaque_plotobjects_20200221.rda")

  combvec=paste0(allcl$donor,"_",allcl$roi)
  batchcol=rep("grey",length(combvec))
  batchcol[intersect(grep("Q17",combvec),grep("Magno",combvec))]="#0000FF"
  batchcol[intersect(grep("Q18",combvec),grep("Magno",combvec))]="#00AAFF"
  batchcol[intersect(grep("Q17",combvec),grep("Parvo",combvec))]="red"
  batchcol[intersect(grep("Q18",combvec),grep("Parvo",combvec))]="pink"

  pdf("macaque_umap_plot.pdf",useDingbats=F)
  temp=plot_objects[["all"]]
  temp=RunUMAP(temp,reduction="pca",dims=1:10)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         allcl$cluster_color[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         batchcol[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  
  temp=plot_objects[["mp"]]
  temp=RunUMAP(temp,reduction="pca",dims=1:15)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         allcl$cluster_color[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         batchcol[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  dev.off()
  
  ####proportion barplot###
  propvals=table(allcl$cluster_label,allcl$roi)
  propvals=sweep(propvals,1,rowSums(propvals),"/")*100
  pdf("macaque_proportions.pdf",useDingbats=F,width=10,height=6)
  barplot(t(propvals),col=c("blue","red","grey"),ylab="Fraction of nuclei")
  dev.off()
  
  load("macaque_plotobjects.rda")
  pdf("macaque_dendrograms.pdf")
  temp=plot_objects[["all"]]
  for (dimval in 3:3) {
    meanvals=matrix(0,ncol=dimval,nrow=length(unique(allcl$cluster_label)))
    rownames(meanvals)=unique(allcl$cluster_label)[order(unique(allcl$cluster_id))]
    for (ii in rownames(meanvals)) {
      for (jj in 1:ncol(meanvals)) {
        meanvals[ii,jj]=mean(temp@reductions$pca@cell.embeddings[intersect(allcl$sample_name[allcl$cluster_label==ii],names(temp@active.ident)),jj])
      }
    }
    distmat=dist((meanvals))
    dend=as.dendrogram(hclust(distmat,method="ward"))
    dend=reorder(dend,(1:11)^2)
    plot(dend,hang=-1,main=dimval)
    distmat=dist((meanvals[grep("^O",rownames(meanvals),invert=T),]))
    dend=as.dendrogram(hclust(distmat,method="ward"))
    dend=reorder(dend,(1:11)^2)
    plot(dend,hang=-1,main=dimval)
  }
  rdev.off()
  
}


spec="human"
if (spec=="human") {
  load("collected_data_20191220.rda")
  keepgen=grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|[0-9]P$|-AS|-PS|DUX4L",rownames(datlist[[spec]]),val=T,invert=T)
  keepdons=which(metalist$human$external_donor_name_label %in% c("H200.1023","H200.1025","H200.1030"))
  datlist[[spec]]=datlist[[spec]][keepgen,keepdons]
  datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
  load("human_clusters_subclusters_20200126.rda")
  allcl=data.frame(sample_name=colnames(m.int@assays$RNA),cluster_id=1,cluster_label=0,cluster_color="black",stringsAsFactors = F)
  rownames(allcl)=allcl$sample_name
  allcl[colnames(m.int7@assays$RNA),"cluster_label"]=paste0("GABA1_",m.int7@active.ident)
  allcl[colnames(m.int6@assays$RNA),"cluster_label"]=paste0("GABA2_",m.int6@active.ident)
  allcl[colnames(m.int5@assays$RNA),"cluster_label"]=paste0("K",m.int5@active.ident)
  allcl[colnames(m.int4@assays$RNA),"cluster_label"]=paste0("MP",m.int4@active.ident)
  allcl[colnames(m.int3@assays$RNA)[m.int3@active.ident %in% c(4)],"cluster_label"]=paste0("O",m.int3@active.ident[m.int3@active.ident %in% c(4)])
  table(allcl$cluster_label)
  ##remove outlier clusters - glia and low quality cells###
  allcl=allcl[!(allcl$cluster_label %in% c(0,"O4")),]
  ###merge clusters with sex-specific genes
  allcl$cluster_label=gsub("K0","K2",allcl$cluster_label)
  allcl$cluster_label=gsub("MP0","MP1",allcl$cluster_label)
  allcl$cluster_label=gsub("GABA1_0","GABA1",allcl$cluster_label)
  allcl$cluster_label=gsub("GABA1_1|GABA1_2","GABA2",allcl$cluster_label)
  allcl$cluster_label=gsub("GABA2_0|GABA2_1|GABA2_2","GABA3",allcl$cluster_label)
  allcl$cluster_label=gsub("GABA2_3","GABA4",allcl$cluster_label)
  table(allcl$cluster_label)
  
  

  
  ###assign order & colors###
  allcl$cluster_id=match(allcl$cluster_label,c("GABA1","GABA2","GABA3","GABA4","K1","K2","MP1","O1"))#,"O1","O2")) 
  colvec=c("#6A3D9A","#F781BF","#E31A1C","#FF7F00","#33A02C","#00DD00","#377EB8","#00CCDD","#DDDDDD")#,"#888888")
  allcl$cluster_color=colvec[allcl$cluster_id]
  
  #allcl$cluster_id=as.numeric(as.factor(allcl$cluster_label))
  #allcl$cluster_color=rainbow(length(unique(allcl$cluster_label)),start=0,end=4/6)[as.numeric(as.factor(allcl$cluster_label))]
  
  allcl$donor=metalist[[spec]]$external_donor_name_label[match(allcl$sample_name,rownames(metalist[[spec]]))]
  allcl$roi=metalist[[spec]]$roi_label[match(allcl$sample_name,rownames(metalist[[spec]]))]
  save(allcl,file="human_allcl_20200126.rda")
  
  sampdat=data.frame(sample_name=allcl$sample_name,t(as.matrix(datlist[[spec]][,allcl$sample_name])))
  pdf("human_violin_plot_20200126.pdf",height=10,width=3)
  group_violin_plot(sampdat,allcl,grouping="cluster",
                    genes=c("SNAP25","GAD2","SLC17A6","LAMP5","NXPH2","TRPC4","ALCAM","KRT80","PVALB","HGF",
                            "PRKCG","CALB1","CASQ2","GRB14","PENK","BCHE",
                            "BTNL9"),
                    log_scale=F)
  dev.off()
  
  
  table(paste(allcl$donor,allcl$roi),allcl$cluster_label)
  ####UMAP plots###
  ###new overall umap plot###
  plot_objects=list()
  subcells=list()
  subcells[["all"]]=1:nrow(allcl)
  subcells[["mp"]]=grep("^MP",allcl$cluster_label)
  subdat=datlist[[spec]][,allcl$sample_name]
  subdat=sweep(subdat,2,Matrix::colSums(subdat),"/")*10^6
  for (nam in names(subcells)) {
    keepcells=subcells[[nam]]
    m.overall <- CreateSeuratObject(counts = subdat[,allcl$sample_name[keepcells]], meta.data = allcl[keepcells,])
    m.overall <- NormalizeData(object =m.overall)
    m.overall <- ScaleData(object =m.overall, verbose = FALSE)
    m.overall <- FindVariableFeatures(object =m.overall, selection.method = "vst", top.genes = 2000)
    m.overall <- RunPCA(object =m.overall, npcs = 30, verbose = FALSE)
    m.overall <- RunTSNE(object =m.overall, reduction = "pca", dims = 1:15)
    m.overall <- RunUMAP(object =m.overall, reduction = "pca", dims = 1:15)
    plot_objects[[nam]]=m.overall
    print(nam)
  }
  save(plot_objects,file="human_plotobjects_20200126.rda")
  
  combvec=paste0(allcl$donor,"_",allcl$roi)
  batchcol=rep("grey",length(combvec))
  batchcol[intersect(grep("H17",combvec),grep("MC",combvec))]="#0000FF"
  batchcol[intersect(grep("1023",combvec),grep("MC",combvec))]="#0044FF"
  batchcol[intersect(grep("1024",combvec),grep("MC",combvec))]="#0088FF"
  batchcol[intersect(grep("1025",combvec),grep("MC",combvec))]="#00AAFF"
  batchcol[intersect(grep("1030",combvec),grep("MC",combvec))]="#00FFFF"
  
  batchcol[intersect(grep("H17",combvec),grep("PC",combvec))]="#FF0000"
  batchcol[intersect(grep("1023",combvec),grep("PC",combvec))]="#AA0000"
  batchcol[intersect(grep("1024",combvec),grep("PC",combvec))]="#880000"
  batchcol[intersect(grep("1025",combvec),grep("PC",combvec))]="#440000"
  batchcol[intersect(grep("1030",combvec),grep("PC",combvec))]="#440044"
  
  batchcol[intersect(grep("H17",combvec),grep("KC",combvec))]="#00FF00"
  batchcol[intersect(grep("1023",combvec),grep("KC",combvec))]="#00AA00"
  batchcol[intersect(grep("1024",combvec),grep("KC",combvec))]="#008800"
  batchcol[intersect(grep("1025",combvec),grep("KC",combvec))]="#004400"
  batchcol[intersect(grep("1030",combvec),grep("KC",combvec))]="#444400"
  
  
  pdf("human_umap_plot_20200126.pdf",useDingbats=F)
  temp=plot_objects[["all"]]
  dimval=26
  temp=RunUMAP(temp,reduction="pca",dims=1:dimval)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         allcl$cluster_color[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2",main=dimval)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         batchcol[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  temp=plot_objects[["mp"]]
  temp=RunUMAP(temp,reduction="pca",dims=1:15)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         allcl$cluster_color[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         batchcol[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  dev.off()
  
  ####proportion barplot###
  propvals=table(allcl$cluster_label,gsub("-a|-p","",allcl$roi))
  propvals=sweep(propvals,1,rowSums(propvals),"/")*100
  pdf("human_proportions_20200126.pdf",useDingbats=F,width=10,height=6)
  barplot(t(propvals[,c(2,3,1)]),col=c("blue","red","grey"),ylab="Fraction of nuclei")
  dev.off()
  
  ###hierarchical plot###
  pdf("human_dendrogram_20200126.pdf")
  temp=plot_objects[["all"]]
  meanvals=matrix(0,ncol=5,nrow=length(unique(allcl$cluster_label)))
  rownames(meanvals)=unique(allcl$cluster_label)[order(unique(allcl$cluster_id))]
  for (ii in rownames(meanvals)) {
    for (jj in 1:ncol(meanvals)) {
        meanvals[ii,jj]=mean(temp@reductions$pca@cell.embeddings[allcl$sample_name[allcl$cluster_label==ii],jj])
    }
  }
  distmat=dist((meanvals))
  dend=as.dendrogram(hclust(distmat,method="ward"))
  dend=reorder(dend,(1:7)^2)
  plot(dend,hang=-1)
  dev.off()
}


spec="mouse"
if (spec=="mouse") {
  load("collected_data_20191220.rda")
  keepgen=grep("^LOC|^mt-|^Rpl|^Rps",rownames(datlist[[spec]]),val=T,invert=T)
  datlist[[spec]]=datlist[[spec]][keepgen,]
  datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6
  load("mouse_clusters_subclusters.rda")
  allcl=data.frame(sample_name=colnames(m.int@assays$RNA),cluster_id=1,cluster_label=0,cluster_color="black",stringsAsFactors = F)
  rownames(allcl)=allcl$sample_name
  allcl[colnames(m.int6@assays$RNA)[m.int6@active.ident %in% c(0:2)],"cluster_label"]=paste0("Chrna",m.int6@active.ident[m.int6@active.ident %in% c(0:2)])
  allcl[colnames(m.int7@assays$RNA)[m.int7@active.ident %in% c(1,3:7)],"cluster_label"]=paste0("GABA",m.int7@active.ident[m.int7@active.ident %in% c(1,3:7)])
  allcl[colnames(m.int5@assays$RNA),"cluster_label"]="LP"
  allcl[colnames(m.int4@assays$RNA),"cluster_label"]="LGD"
  allcl[colnames(m.int2@assays$RNA)[m.int2@active.ident %in% c(13)],"cluster_label"]="LGV"
  table(allcl$cluster_label)
  ###remove non-neuronal cells#
  allcl=allcl[allcl$cluster_label!="0",]
  ###merge clusters
  allcl$cluster_label=gsub("Chrna0","Chrna3",allcl$cluster_label)
  #allcl$cluster_label=gsub("GABA0","GABA2",allcl$cluster_label)
  
  
  ###assign order & colors###
  allcl$cluster_id=match(allcl$cluster_label,c("Chrna1","Chrna2","Chrna3",paste0("GABA",c(1,4,5,3,6,7)),"LGV","LP","LGD"))#,"O1","O2")) 
  colvec=c("#DD6091","#FF88AD","#FDDCEE","#E6BB3D","#E6E93D","#FFFF00","#FFB307","#FE8A3E","#C15313","#A6E6E3","#00DDC5","#007B9D")
  #colvec=c("#6A3D9A","#F781BF","#E31A1C","#FF7F00","#33A02C","#00DD00","#377EB8","#00CCDD","#DDDDDD")#,"#888888")
  allcl$cluster_color=colvec[allcl$cluster_id]
  

  #allcl$donor=metalist[[spec]]$external_donor_name_label[match(allcl$sample_name,rownames(metalist[[spec]]))]
  allcl$roi=metalist[[spec]]$roi[match(allcl$sample_name,rownames(metalist[[spec]]))]
  save(allcl,file="mouse_allcl.rda")
  
  sampdat=data.frame(sample_name=allcl$sample_name,t(as.matrix(datlist[[spec]][,allcl$sample_name])))
  pdf("mouse_violin_plot.pdf",height=10,width=3)
  group_violin_plot(sampdat,allcl,grouping="cluster",
                    genes=c("Snap25","Gad2","Slc17a6","Dmbx1","Chrna6","Dscaml1","Alyref","Dapk3","Zmynd19","Pax7",
                            "Pvalb","Spata20","Rspo2","Htr1a","Eltd1","Rspo3",
                            "Neurod6","Slc17a7","Glra3","Scn4b"),
                    log_scale=F)
  dev.off()
  
  
  table(paste(allcl$donor,allcl$roi),allcl$cluster_label)
  ####UMAP plots###
  ###new overall umap plot###
  plot_objects=list()
  subcells=list()
  subcells[["all"]]=1:nrow(allcl)
  subcells[["LGD"]]=grep("LGD",allcl$cluster_label)
  subdat=datlist[[spec]][,allcl$sample_name]
  subdat=sweep(subdat,2,Matrix::colSums(subdat),"/")*10^6
  for (nam in names(subcells)) {
    keepcells=subcells[[nam]]
    m.overall <- CreateSeuratObject(counts = subdat[,allcl$sample_name[keepcells]], meta.data = allcl[keepcells,])
    m.overall <- NormalizeData(object =m.overall)
    m.overall <- ScaleData(object =m.overall, verbose = FALSE)
    m.overall <- FindVariableFeatures(object =m.overall, selection.method = "vst", top.genes = 2000)
    m.overall <- RunPCA(object =m.overall, npcs = 30, verbose = FALSE)
    m.overall <- RunTSNE(object =m.overall, reduction = "pca", dims = 1:15)
    m.overall <- RunUMAP(object =m.overall, reduction = "pca", dims = 1:15)
    plot_objects[[nam]]=m.overall
    print(nam)
  }
  save(plot_objects,file="mouse_plotobjects.rda")
  
  combvec=paste0(allcl$roi)
  batchcol=rep("grey",length(combvec))
  batchcol[grep("Core",combvec)]="#580077"
  batchcol[grep("Shell",combvec)]="#D066E0"
  batchcol[grep("LGv",combvec)]="#8FBE3F"
  batchcol[grep("LP",combvec)]="#1C75BC"
  
  
  
  pdf("mouse_umap_plot.pdf",useDingbats=F)
  temp=plot_objects[["all"]]
  dimval=15
  temp=RunUMAP(temp,reduction="pca",dims=1:dimval)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         allcl$cluster_color[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2",main=dimval)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         batchcol[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  temp=plot_objects[["LGD"]]
  temp=RunUMAP(temp,reduction="pca",dims=1:15)
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         allcl$cluster_color[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  plot(temp@reductions$umap@cell.embeddings,pch=19,col=
         batchcol[match(rownames(temp@reductions$tsne@cell.embeddings),allcl$sample_name)],cex=0.5,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
  dev.off()
  
  ####proportion barplot###
  propvals=table(allcl$cluster_label,gsub("Mouse ","",allcl$roi))
  propvals=sweep(propvals,1,rowSums(propvals),"/")*100
  pdf("mouse_proportions.pdf",useDingbats=F,width=10,height=6)
  barplot(t(propvals[c("Chrna1","Chrna2","Chrna3",paste0("GABA",c(1,4,5,3,6,7)),"LGV","LP","LGD")
                     ,c(1,2,4,3)]),col=c("#580077","#D066E0","#1C75BC","#8FBE3F"),ylab="Fraction of nuclei")
  dev.off()

  ###hierarchical plot###
  pdf("mouse_dendrogram.pdf")
  temp=plot_objects[["all"]]
  meanvals=matrix(0,ncol=30,nrow=length(unique(allcl$cluster_label)))
  rownames(meanvals)=unique(allcl$cluster_label)[order(unique(allcl$cluster_id))]
  for (ii in rownames(meanvals)) {
    for (jj in 1:ncol(meanvals)) {
      meanvals[ii,jj]=mean(temp@reductions$pca@cell.embeddings[allcl$sample_name[allcl$cluster_label==ii],jj])
    }
  }
  dimval=10
  distmat=dist((meanvals[c("LP","LGV","LGD"),1:dimval]))
  dend=as.dendrogram(hclust(distmat,method="ward"))
  plot(dend,hang=-1,main=dimval)
  distmat=dist((meanvals[grep("GABA|Chrna",rownames(meanvals)),1:dimval]))
  dend=as.dendrogram(hclust(distmat,method="ward"))
  plot(dend,hang=-1,main=dimval)
  #distmat=dist((meanvals[,1:dimval]))
  #dend=as.dendrogram(hclust(distmat,method="ward"))
  #plot(dend,hang=-1,main=dimval)
  dev.off()
}




