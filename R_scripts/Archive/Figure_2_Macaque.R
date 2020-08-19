
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
require(Matrix)
require(scrattch.vis)
library(gridExtra)
source("doublet.finder.r")
source("hicat_modules.r")


st=format(Sys.time(), "%Y%m%d_%H%M_")


########################
## Figure 2: Macaque plots
########################


load("collected_data_20191220.rda")
load("macaque_plotobjects_20200221.rda")

load("macaque_allcl.rda")
cluster_palette <- read.csv("cluster_palette.csv")

allcl$class <- cluster_palette$class[cluster_palette$species_label == "Macaque"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Macaque"])]

allcl$cluster_label <- cluster_palette$new_cluster_label[cluster_palette$species_label == "Macaque"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Macaque"])]

allcl <- allcl[allcl$class != "Low Quality",]

### dendrograms ###

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
  dev.off()
  
}




### violin plot ###
spec="macaque"
  orthos=read.csv("ortholog_table_20191122.csv",as.is=T)
  nonlocgenes=orthos$rhesus_symbol[intersect(grep("^LOC",orthos$human_symbol,invert=T),grep("^LOC",orthos$rhesus_symbol))]
  keepgen=setdiff(rownames(datlist[[spec]]),setdiff(grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|^TRNA",rownames(datlist[[spec]]),val=T),nonlocgenes))
  datlist[[spec]]=datlist[[spec]][keepgen,]
  datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6


###assign order & colors###
allcl$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Macaque"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Macaque"])]

sampdat=data.frame(sample_name=allcl$sample_name,t(as.matrix(datlist[[spec]][,allcl$sample_name])))

plot <- group_violin_plot(sampdat,allcl,group_order=c(1,2,3,4,5,6,9,10,7,8),
					grouping="cluster",
                    genes=c("SNAP25","GAD2","SLC17A6","DLX1","NPY","TNFAIP8L3",
                    	"LAMP5","KRT80","DSEL","SLITRK2","COL2A1","NXPH2",
                    	"PRKCG","CAMK2A" ,
                    	#"G0S2",
                    	"PENK","SGCD","DAB2","LVRN",
                    	"CASQ2","IL15RA","FOXP2","EBF1","FGF11","ROBO3","GRIK3",
                    	"PLPP4","RUFY4","PTPN20"),
                    log_scale=F,
                    label_height=15,
                    max_width=25)
ggsave(plot,file=paste0(st,"macaque_violin_plot.pdf"),height=10,width=3, useDingbats=FALSE)



### UMAP plots ###

temp=plot_objects[["all"]]
temp <- subset(temp, subset= cluster_label != "O3")
dimval=26
temp=RunUMAP(temp,reduction="pca",dims=1:dimval)
umap_all <- temp@reductions$umap@cell.embeddings

temp=plot_objects[["mp"]]
temp=RunUMAP(temp,reduction="pca",dims=1:15)
umap_mp <- temp@reductions$umap@cell.embeddings
save(umap_all, umap_mp, file="macaque_umap.rda")  

load("macaque_umap.rda")

colnames(umap_all) = c("Lim1","Lim2")
clcol <- droplevels(setNames(cluster_palette$cluster_color[cluster_palette$species_label=="Macaque"], cluster_palette$cluster_id[cluster_palette$species_label=="Macaque"]))
g1= plot_RD_meta(umap_all, factor(allcl[row.names(umap_all),] %>% pull(cluster_id),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=1.5)
g1 = g1+	ggtitle("cluster") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
ggsave(g1, file=paste0(st,"macaque_umap_cluster.pdf"), height=5,width=5, useDingbats=FALSE)


donor_palette <- read.csv("donor_palette.csv")
allcl$donor_roi <- paste0(allcl$donor,"_",allcl$roi)
dcol <- setNames(donor_palette$donor_roi_color[donor_palette$species=="Macaque"], donor_palette$donor_roi[donor_palette$species=="Macaque"])
g2= plot_RD_meta(umap_all, factor(allcl[row.names(umap_all),] %>% pull(donor_roi),levels=names(dcol)), meta.col=dcol,alpha=0.5, cex=1.5)
g2 = g2 + ggtitle("donor/roi") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
ggsave(g2, file=paste0(st,"macaque_umap_donor.roi.pdf"), height=5,width=5, useDingbats=FALSE)

g3 <- g2 + facet_wrap(~meta)

library(stringr)
#allcl$dissection_roi <- str_extract(allcl$roi, "[^-]+")
g2[["data"]]$roi <- allcl$roi

g4 <- g2 + facet_wrap(~roi)


colnames(umap_mp) = c("Lim1","Lim2")
g5= plot_RD_meta(umap_mp, factor(allcl[row.names(umap_mp),] %>% pull(cluster_id),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g5 = g5+	ggtitle("cluster") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
ggsave(g5, file=paste0(st,"macaque_umap_mp_cluster.pdf"), height=5,width=5, useDingbats=FALSE)

g6= plot_RD_meta(umap_mp, factor(allcl[row.names(umap_mp),] %>% pull(donor_roi),levels=names(dcol)), meta.col=dcol,alpha=0.5, cex=2)
g6 = g6 + ggtitle("donor/roi") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
ggsave(g6, file=paste0(st,"macaque_umap_mp_donor.roi.pdf"), height=5,width=5, useDingbats=FALSE)

g7 <- g6 + facet_wrap(~meta)


ggsave(file=paste0(st,"macaque_umap.pdf"),marrangeGrob(list(g1,g2,g3,g4,g5,g6,g7),nrow = 1, ncol=1), height=5,width=6)



####proportion barplot###
allcl$cluster_label <- factor(allcl$cluster_label, levels=labels(dend))

  propvals=table(allcl$cluster_label,allcl$roi)
  propvals=sweep(propvals,1,rowSums(propvals),"/")*100
  pdf(paste0(st,"macaque_proportions.pdf"),useDingbats=F,width=10,height=6)
  barplot(t(propvals),col=c("#174596","#FF7F00","#BEBEBE"),ylab="Fraction of nuclei")
  dev.off()
  




