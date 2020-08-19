
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

load("collected_data_20191220.rda")

st=format(Sys.time(), "%Y%m%d_%H%M_")



########################
## Figure 4: Mouse plots
########################

spec="mouse"

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
  
  

load("mouse_allcl.rda")
cluster_palette <- read.csv("cluster_palette.csv")

###assign order & colors###
allcl$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Mouse"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Mouse"])]

sampdat=data.frame(sample_name=allcl$sample_name,t(as.matrix(datlist[[spec]][,allcl$sample_name])))


###violin plot###
plot <- group_violin_plot(sampdat,allcl,
          grouping="cluster",
          genes=c("Snap25","Gad2","Slc17a6","Dmbx1","Chrna6",
                "Dscaml1","Alyref","Dapk3","Zmynd19","Pax7",
                "Pvalb","Spata20","Rspo2","Htr1a","Eltd1",
                "Rspo3","Neurod6","Slc17a7","Glra3","Scn4b"),
                    log_scale=F,
                    label_height=15,
                    max_width=25)
ggsave(plot,file=paste0(st,"mouse_violin_plot.pdf"),height=10,width=4, useDingbats=FALSE)



####UMAP plots###
load("mouse_plotobjects.rda")

  temp=plot_objects[["all"]]
  dimval=15
#nrs <-c(41,44,45,46)
  temp=RunUMAP(temp,reduction="pca",dims=1:dimval, seed.use=44)
umap_all <- temp@reductions$umap@cell.embeddings
plot(umap_all)

  temp=plot_objects[["LGD"]]
  temp=RunUMAP(temp,reduction="pca",dims=1:15,, seed.use=46)
umap_lgd <- temp@reductions$umap@cell.embeddings
plot(umap_lgd)
save(umap_all, umap_lgd, file="mouse_umap.rda")  


g.list <- list()
for(i in 1:50) {
  print(i)
  temp=plot_objects[["LGD"]]
  dimval=15
  #nrs <-c(41,44,45,46)
  temp=RunUMAP(temp,reduction="pca",dims=1:dimval, seed.use=i)
  umap_lgd <- as.data.frame(temp@reductions$umap@cell.embeddings)
  p1 <- ggplot(umap_lgd, aes(x=UMAP_1, y=UMAP_2)) + geom_point() + ggtitle(label=i)
  g.list[[i]] <- p1
}
ggsave(file=paste0(st,"mouse_umap.test.pdf"),gridExtra::marrangeGrob(g.list,nrow = 2, ncol=2), height=8,width=8)



[1] "Mouse LGN - Core"  "Mouse LGN - Shell"
[3] "Mouse LGv"         "Mouse LP"


load("mouse_umap.rda")
colnames(umap_all) = c("Lim1","Lim2")
clcol <- droplevels(setNames(cluster_palette$cluster_color[cluster_palette$species_label=="Mouse"], cluster_palette$cluster_id[cluster_palette$species_label=="Mouse"]))
g1= plot_RD_meta(umap_all, factor(allcl[row.names(umap_all),] %>% pull(cluster_id),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g1 = g1+  ggtitle("cluster") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

roi_palette <- read.csv("roi_palette.csv",stringsAsFactors = FALSE)
rcol <- setNames(roi_palette$roi_color[roi_palette$species=="Mouse"], roi_palette$dissection_roi[roi_palette$species=="Mouse"])
#rcol <- rcol[!duplicated(rcol)]
allcl$dissection_roi <- roi_palette$dissection_roi[match(allcl$roi, roi_palette$roi)]

g2= plot_RD_meta(umap_all, factor(allcl[row.names(umap_all),] %>% pull(dissection_roi),levels=names(rcol)), meta.col=rcol,alpha=0.5, cex=2)
g2 = g2 + ggtitle("roi") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

g3 <- g2 + facet_wrap(~meta, nrow=1)

#ggsave(file=paste0(st,"mouse_umap.pdf"),gridExtra::marrangeGrob(list(g1,g2),nrow = 1, ncol=1), height=5,width=5)
#ggsave(file=paste0(st,"mouse_umap.split.pdf"),g3, height=5,width=20)

colnames(umap_lgd) = c("Lim1","Lim2")
g4= plot_RD_meta(umap_lgd, factor(allcl[row.names(umap_lgd),] %>% pull(cluster_id),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g4 = g4+  ggtitle("cluster") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

g5= plot_RD_meta(umap_lgd, factor(allcl[row.names(umap_lgd),] %>% pull(dissection_roi),levels=names(rcol)), meta.col=rcol,alpha=0.5, cex=2)
g5 = g5 + ggtitle("roi") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

g6 <- g5 + facet_wrap(~meta, nrow=1)


ggsave(file=paste0(st,"mouse_umap.pdf"),gridExtra::marrangeGrob(list(g1,g2,g4,g5),nrow = 1, ncol=1), height=5,width=5)
ggsave(file=paste0(st,"mouse_umap.split.pdf"),gridExtra::marrangeGrob(list(g3,g6),nrow = 1, ncol=1), height=5,width=20)



  
  ####proportion barplot###
  propvals=table(allcl$cluster_label,gsub("Mouse ","",allcl$roi))
  propvals=sweep(propvals,1,rowSums(propvals),"/")*100
  pdf(paste0(st,"mouse_proportions.pdf"),useDingbats=F,width=10,height=6)
  barplot(t(propvals[c("Chrna1","Chrna2","Chrna3",paste0("GABA",c(1,4,5,3,6,7)),"LGV","LP","LGD")
                     ,c(1,2,4,3)]),col=c("#580077","#D066E0","#7285A5","#8FBE3F"),ylab="Fraction of nuclei")
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





  ####de genes heatmap###

setwd("C:/Users/menonv/Dropbox (HHMI)/AIBS/Transcriptomics/Manuscripts/LGd/continuum_analysis/")
require(feather)
require(dplyr)
require(WGCNA)
require(gplots)
require(edgeR)

setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")



#mdat=read_feather("../Data/mouse_LGN_STAR/data.feather")
mdat=read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_LGN_20180105_STAR/data.feather")
manno=read_feather("../old/Data/mouse_LGN_STAR/anno.feather")
class1=which(manno$cluster_label=="LGd_Slc17a6_Scn4b" & manno$roi_label=="LGd_Core")
class2=which(manno$cluster_label=="LGd_Slc17a6_Scn4b" & manno$roi_label=="LGd_Shell")
classvals=as.factor(rep(c(1,2),times=c(length(class1),length(class2))))
startmat=t(sweep(mdat[c(class1,class2),-ncol(mdat)],2,rowSums(mdat[c(class1,class2),-ncol(mdat)]),"/")*10^6)
e_design=model.matrix(~classvals)
y2 = DGEList(counts=t(mdat[c(class1,class2),-ncol(mdat)]))
y2 = estimateDisp(y2, e_design)
fit = glmQLFit(y2, e_design)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
outval2=topTags(qlf.2vs1,n=nrow(startmat),p.value=1)
outval2$table=outval2$table[intersect(rownames(outval2$table),rownames(startmat)),]
mean1=apply(startmat[rownames(outval2$table),classvals==1],1,mean)
mean2=apply(startmat[rownames(outval2$table),classvals==2],1,mean)
frac1=rowSums(startmat[rownames(outval2$table),classvals==1]>0)/length(class1)
frac2=rowSums(startmat[rownames(outval2$table),classvals==2]>0)/length(class2)
outmat=cbind(outval2$table,mean1,mean2,frac1,frac2)

suboutmat=outmat[outmat$FDR<0.001 & apply(outmat[,c("frac1","frac2")],1,max)>0.25,]
suboutmat=suboutmat[-grep("Gm|LOC",rownames(suboutmat)),]


keepgen=c("Stxbp6","Phldb2","Tnc","Cdc42ep4","Il1rapl2","Drd1","Necab1","Rgs6","Mgat4c","Mbp","Adcy2","Angpt1","Fam19a1","Dcx","Tmem236","Scnn1a","Zfp385c","Proser2","G630055G22Rik","Wwtr1","Iyd","Bace2","Slc35f4","Acp5","Etl4","Fam65b","Lama2","Sec14l4","Enpep","Cped1")
keepgen<-append(keepgen, c("Pvalb", "Calb1"))

submat <- suboutmat[keepgen,]
submat <- submat[order(-submat$logFC),]

keepgen1=rownames(submat)[order(-submat$logFC)[1:16]]
keepgen2=rownames(submat)[order(submat$logFC)[1:16]]

plotmat=startmat[c(keepgen1,keepgen2),]
colvec=rep(c("#580077","#B000CC"),times=c(length(class1),length(class2)))
breakvals=seq(from=-2,to=2,length.out=201)
colnames(plotmat)=rep("",ncol(plotmat))

pdf(paste0(st, "core_shell_de_genes.pdf"),useDingbats=F, height=5, width=8)
heatmap.2(log10(plotmat+1),breaks=breakvals,col=bluered(200),scale="row",trace='none',Rowv=F,Colv=F,ColSideColors=colvec)
dev.off()