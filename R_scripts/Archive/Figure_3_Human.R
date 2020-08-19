
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
## Figure 3: Human plots
########################


spec="human"
load("collected_data_20191220.rda")
keepgen=grep("^LOC|^MT-|^RPL|^RPS|^SNAR|^CT45|^MIR|^MTND|[0-9]P$|-AS|-PS|DUX4L",rownames(datlist[[spec]]),val=T,invert=T)
keepdons=which(metalist$human$external_donor_name_label %in% c("H200.1023","H200.1025","H200.1030"))
datlist[[spec]]=datlist[[spec]][keepgen,keepdons]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6

load("human_allcl_20200126.rda")
cluster_palette <- read.csv("cluster_palette.csv")

###assign order & colors###
allcl$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Human"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Human"])]

sampdat=data.frame(sample_name=allcl$sample_name,t(as.matrix(datlist[[spec]][,allcl$sample_name])))

genes<-c("BRD4", "CAV1", "EEF1A2", "INA", "KCNA1", "NEFH", "NEFL", "PPP2R2C", "SFRP2", "TCF7L2"  )

plot <- group_violin_plot(sampdat,allcl,
					grouping="cluster",
                    genes=c("SNAP25","GAD2","SLC17A6","LAMP5","NXPH2",
                    		"TRPC4","ALCAM","KRT80","PVALB","HGF",
                            "PRKCG","CALB1","CASQ2","GRB14","PENK",
                            "BCHE","BTNL9"),
                    log_scale=F,
                    label_height=15,
                    max_width=25)
ggsave(plot,file=paste0(st,"human_violin_plot.pdf"),height=10,width=3, useDingbats=FALSE)
  
  
table(paste(allcl$donor,allcl$roi),allcl$cluster_label)

####UMAP plots###
load("human_plotobjects_20200126.rda")

temp=plot_objects[["all"]]
dimval=26
temp=RunUMAP(temp,reduction="pca",dims=1:dimval)
umap_all <- temp@reductions$umap@cell.embeddings

temp=plot_objects[["mp"]]
temp=RunUMAP(temp,reduction="pca",dims=1:15)
umap_mp <- temp@reductions$umap@cell.embeddings
save(umap_all, umap_mp, file="human_umap.rda")  

load("human_umap.rda")
colnames(umap_all) = c("Lim1","Lim2")
clcol <- droplevels(setNames(cluster_palette$cluster_color[cluster_palette$species_label=="Human"], cluster_palette$cluster_id[cluster_palette$species_label=="Human"]))
g1= plot_RD_meta(umap_all, factor(allcl[row.names(umap_all),] %>% pull(cluster_id),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g1 = g1+	ggtitle("cluster") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

donor_palette <- read.csv("donor_palette.csv")
dcol <- droplevels(setNames(donor_palette$donor_roi_color[donor_palette$species=="Human"], donor_palette$donor[donor_palette$species=="Human"]))
g2= plot_RD_meta(umap_all, factor(allcl[row.names(umap_all),] %>% pull(donor),levels=names(dcol)), meta.col=dcol,alpha=0.5, cex=2)
g2 = g2 + ggtitle("donor") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

g3 <- g2 + facet_wrap(~meta)

library(stringr)
allcl$dissection_roi <- str_extract(allcl$roi, "[^-]+")
g2[["data"]]$roi <- allcl$dissection_roi

g4 <- g2 + facet_wrap(~roi)

colnames(umap_mp) = c("Lim1","Lim2")
g5= plot_RD_meta(umap_mp, factor(allcl[row.names(umap_mp),] %>% pull(cluster_id),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g5 = g5+	ggtitle("cluster") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

g6= plot_RD_meta(umap_mp, factor(allcl[row.names(umap_mp),] %>% pull(donor),levels=names(dcol)), meta.col=dcol,alpha=0.5, cex=2)
g6 = g6 + ggtitle("donor") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

g7 <- g6 + facet_wrap(~meta)
mp.cl <- allcl[grep("MP",allcl$cluster_label),]
g6[["data"]]$roi <- mp.cl$dissection_roi
g8 <- g6 + facet_wrap(~roi)


#ggsave(file=paste0(st,"human_umap.pdf"),gridExtra::marrangeGrob(list(g1,g2,g5,g6),nrow = 1, ncol=1), height=5,width=5)
#ggsave(file=paste0(st,"human_umap.split.pdf"),gridExtra::marrangeGrob(list(g3,g4,g7,g8),nrow = 1, ncol=1), height=5,width=15)


##extra
roi_palette <- read.csv("roi_palette.csv",stringsAsFactors = FALSE)
rcol <- setNames(roi_palette$roi_color[roi_palette$species=="Human"], roi_palette$dissection_roi[roi_palette$species=="Human"])
rcol <- rcol[!duplicated(rcol)]
allcl$dissection_roi <- str_extract(allcl$roi, "[^-]+")

g9= plot_RD_meta(umap_all, factor(allcl[row.names(umap_all),] %>% pull(dissection_roi),levels=names(rcol)), meta.col=rcol,alpha=0.5, cex=2)
g9 = g9 + ggtitle("roi") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
g10 <- g9 + facet_wrap(~meta)


g9[["data"]]$donor <- allcl$donor
g11 <- g9 + facet_wrap(~donor)

colnames(umap_mp) = c("Lim1","Lim2")
g12= plot_RD_meta(umap_mp, factor(allcl[row.names(umap_mp),] %>% pull(dissection_roi),levels=names(rcol)), meta.col=rcol,alpha=0.5, cex=2)
g12 = g12 + ggtitle("roi") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
g13 <- g12 + facet_wrap(~meta)

mp.cl <- allcl[grep("MP",allcl$cluster_label),]
g12[["data"]]$donor <- mp.cl$donor
g14 <- g12 + facet_wrap(~donor)

ggsave(file=paste0(st,"human_umap.pdf"),gridExtra::marrangeGrob(list(g1,g2,g5,g6,g9,g12),nrow = 1, ncol=1), height=5,width=5)
ggsave(file=paste0(st,"human_umap.split.pdf"),gridExtra::marrangeGrob(list(g3,g4,g7,g8,g10,g11,g13,g14),nrow = 1, ncol=1), height=5,width=15)


temp=plot_objects[["mp"]]
temp=RunUMAP(temp,reduction="pca",dims=1:15, n.components=1)
umap_mp1 <- temp@reductions$umap@cell.embeddings

# plot 1d umap as grouped scatter plot per donor colored by dissection
umap_mp1 <- tibble::rownames_to_column(as.data.frame(umap_mp1), var="sample_name")
umap_mp1 <- left_join(umap_mp1, allcl)

ggplot(umap_mp1, aes(y=donor, x=UMAP_1)) + 
	geom_jitter(aes(color = dissection_roi),
    				stat = "identity", 
    				position = position_jitter(0.2)) +
	scale_color_manual(values=rcol)

ggplot(umap_mp1, aes(y=donor, x=UMAP_1)) + 
	geom_tile(aes(color = dissection_roi)) +
	scale_color_manual(values=rcol)

g15 <- ggplot(umap_mp1, aes(y=donor, x=UMAP_1)) + 
	geom_point(aes(color = dissection_roi), shape=124, size=8) +
	scale_color_manual(values=rcol)+
	theme_void() + 
 	theme(legend.position="none") +
 	theme(axis.text.y=element_text())
ggsave(g15, file=paste0(st,"human_mp_1dumap_donor_roi.pdf"), width=6, height=2.4)

g16 <- ggplot(umap_mp1, aes(y=dissection_roi, x=UMAP_1,color = dissection_roi)) + 
  #geom_jitter(stat = "identity", position = position_jitter(0.2)) +
  geom_quasirandom(groupOnX=FALSE)+
  scale_color_manual(values=rcol) +
  xlab('UMAP 1')+
  ylab('Dissection roi')+
  facet_wrap(~donor, ncol=1)+
  theme_bw()
ggsave(g16, file=paste0(st,"human_mp_1dumap_donor_roi_split.pdf"), width=5, height=7, useDingbats=FALSE)

#### testing various options for plotting the 1d umap result

library(cowplot)
library(dplyr)
library(readr)
cloud.dir <- "//allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/RainCloudPlots-master/tutorial_R"
source(file.path(cloud.dir,"R_rainclouds.R"))
#source(file.path(cloud.dir,"summarySE.R"))

p3 <- ggplot(umap_mp1, aes(x=donor, y=UMAP_1, fill=dissection_roi, color=dissection_roi))+
    geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust = 2, alpha=0.5)+
    geom_point(position = position_jitter(width = .15), size = .75)+
    #geom_quasirandom(dodge.width=0.1)+
    scale_color_manual(values=rcol)+
    scale_fill_manual(values=rcol)+
    ylab('UMAP 1')+
    xlab('Donor')+
    coord_flip()+
    theme_cowplot()+
    guides(fill = FALSE)


ggplot(umap_mp1, aes(x=donor, y=UMAP_1, fill=dissection_roi, color=dissection_roi))+
    geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust = 2, alpha=0.5)+
    geom_quasirandom()+
    scale_color_manual(values=rcol)+
    scale_fill_manual(values=rcol)+
    ylab('UMAP 1')+
    xlab('Donor')+
    coord_flip()+
    theme_cowplot()+
    guides(fill = FALSE)



p10 <- ggplot(umap_mp1, aes(x=donor, y=UMAP_1, fill=dissection_roi, color=dissection_roi)) +
    geom_flat_violin(aes(fill = dissection_roi),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
    geom_point(aes(x = donor, y = UMAP_1, colour = dissection_roi),position = position_jitter(width = .05), size = 1, shape = 20)+
    geom_boxplot(aes(x=donor, y=UMAP_1, fill=dissection_roi),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")
################





####proportion barplot###
  propvals=table(allcl$cluster_label,gsub("-a|-p","",allcl$roi))
  propvals=sweep(propvals,1,rowSums(propvals),"/")*100
  pdf(paste0(st,"human_proportions.pdf"),useDingbats=F,width=10,height=6)
  barplot(t(propvals[,c(2,3,1)]),col=c("#123466","#FF7F00","#235602"),ylab="Fraction of nuclei")
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


########################
## Figure 3: Human-Macaque integrated plots
########################
load("human_macaque_magpar_integrated.rda")

load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl_20200126.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")

####plot clusters###
coords=obj.integrated@reductions$umap@cell.embeddings
int_keepmac=intersect(rownames(coords)[coords[,1]< -5],allcl_mac$sample_name[allcl_mac$cluster_label %in% c("MP1_M","MP2_M")])
int_keephum=intersect(rownames(coords)[coords[,1]< -5],allcl_hum$sample_name[allcl_hum$cluster_label %in% c("MP1_H")])

umap_mp_mh<-tibble::rownames_to_column(as.data.frame(coords), var="sample_name")
allcl <- rbind(allcl_hum, allcl_mac[,c(1:4,6,7)] )

umap_mp_mh <- left_join(umap_mp_mh, allcl)
umap_mp_mh <- umap_mp_mh[umap_mp_mh$UMAP_1 < 0,]
umap_mp_mh <- umap_mp_mh[,c(2,3,1,4:8)] # for plot_RD_meta xy coords need to be in cols 1,2
umap_mp_mh <- umap_mp_mh[umap_mp_mh$cluster_label %in% c("MP1_M","MP2_M", "MP1_H"),]

clcol <- setNames(c("#2196F3","#0d47a1","#90caf9"),c("MP1_H","MP1_M","MP2_M")) 
g1= plot_RD_meta(umap_mp_mh, factor(umap_mp_mh$cluster_label,levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g1 = g1+  ggtitle("cluster") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

umap_mp_m <- umap_mp_mh[umap_mp_mh$cluster_label %in% c("MP1_M","MP2_M"),]
g2= plot_RD_meta(umap_mp_m, factor(umap_mp_m$cluster_label,levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g2 = g2+  ggtitle("cluster") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

umap_mp_h <- umap_mp_mh[umap_mp_mh$cluster_label %in% c("MP1_H"),]
g3= plot_RD_meta(umap_mp_h, factor(umap_mp_h$cluster_label,levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=2)
g3 = g3+  ggtitle("cluster") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

ggsave(file=paste0(st,"macHu_int_umap.pdf"),gridExtra::marrangeGrob(list(g1,g2,g3),nrow = 1, ncol=1), height=5,width=5)


########
## plotting greyed out populations
########

#### plot clusters with other species in grey
umap_mp_m <- umap_mp_mh[umap_mp_mh$cluster_label %in% c("MP1_M","MP2_M"),]
umap_mp_h <- umap_mp_mh[umap_mp_mh$cluster_label %in% c("MP1_H"),]

#clcol <- setNames(c("#2196F3","#0d47a1","#90caf9"),c("MP1_H","MP1_M","MP2_M")) 
clcol <- setNames(c("grey85","#0d47a1","#90caf9"),c("MP1_H","MP1_M","MP2_M")) 

plot1 <- ggplot() +  
          geom_point(data=umap_mp_h, aes(UMAP_1, y=UMAP_2, color=cluster_label), alpha=0.5,size=2, shape=19)  +  
          geom_point(data=umap_mp_m, aes(x=UMAP_1, y=UMAP_2, color=cluster_label),alpha=0.75,size=2, shape=19)+ 
         
          scale_colour_manual(values=clcol)+ theme_void()+
          #guides(colour = guide_legend(override.aes = list(size=5)))+     
          theme(legend.position="none") 
    ggsave("MP.mac.cluster.png",plot1,height=5,width=5)
    ggsave("MP.mac.cluster.pdf",plot1,height=5,width=5)

clcol <- setNames(c("#2196F3","grey85","grey85"),c("MP1_H","MP1_M","MP2_M")) 
plot2 <- ggplot() + 
          geom_point(data=umap_mp_m, aes(x=UMAP_1, y=UMAP_2, color=cluster_label),alpha=0.5,size=2)+ 
          geom_point(data=umap_mp_h, aes(UMAP_1, y=UMAP_2, color=cluster_label), alpha=0.75,size=2)  +  
          scale_colour_manual(values=clcol)+ theme_void()+
          #guides(colour = guide_legend(override.aes = list(size=5)))+     
          theme(legend.position="none") 
    ggsave("MP.hum.cluster.png",plot2,height=5,width=5)
    ggsave("MP.hum.cluster.pdf",plot2,height=5,width=5)


#### gene plot with other species greyed out

select.genes = c("FOXP2", "ROBO2")
norm.dat.m <- as.data.frame(t(as.matrix(datlist[["macaque"]][select.genes,int_keepmac])))
norm.dat.m <- tibble::rownames_to_column(norm.dat.m, var="sample_name")
norm.dat.h <- as.data.frame(t(as.matrix(datlist[["human"]][select.genes,int_keephum])))
norm.dat.h <- tibble::rownames_to_column(norm.dat.h, var="sample_name")

plots.m <- list()
for (gene in select.genes){ 
  print(gene)
  plots.m[[gene]] <- local({ 
    pc.dat <- left_join(umap_mp_m, select(norm.dat.m,"sample_name", gene))

  #  p <- ggplot() + 
  #        geom_point(data=umap_mp_h, aes(x=UMAP_1, y=UMAP_2),color="grey85",size=2) +
  #        geom_point(data=pc.dat, aes(x=UMAP_1, y=UMAP_2,color=pc.dat[,gene]),size=2) +
  #        scale_color_gradient(low="#482677FF",high="yellow") 
   
    p <- ggplot() + 
          geom_point(data=umap_mp_h, aes(x=UMAP_1, y=UMAP_2),color="grey85",size=2, alpha=0.5) + 
          geom_point(data=pc.dat, aes(x=UMAP_1, y=UMAP_2,color=pc.dat[,gene]),size=2, alpha=0.75) +  
          scale_color_gradient(low="gray50",high="red")
    p <- p+ theme_void()+ theme(legend.position="none") 
    
    ggsave(paste0(gene,".MP.mac.pdf"),p,height=5,width=5)
     })
  }

plots.h <- list()
for (gene in select.genes){ 
  print(gene)
  plots.h[[gene]] <- local({ 
    pc.dat <- left_join(umap_mp_h, select(norm.dat.h,"sample_name", gene))

  #  p <- ggplot() + 
  #       geom_point(data=umap_mp_m, aes(x=UMAP_1, y=UMAP_2),color="grey85",size=2) +
  #       geom_point(data=pc.dat, aes(x=UMAP_1, y=UMAP_2,color=pc.dat[,gene]),size=2) +
  #        scale_color_gradient(low="#482677FF",high="yellow") 
   
    p <- ggplot() + 
          geom_point(data=umap_mp_m, aes(x=UMAP_1, y=UMAP_2),color="grey85",size=2, alpha=0.5) + 
          geom_point(data=pc.dat, aes(x=UMAP_1, y=UMAP_2,color=pc.dat[,gene]),size=2,alpha=0.75) +  
          scale_color_gradient(low="gray50",high="red")
    p <- p+ theme_void()+ theme(legend.position="none") 
    
    ggsave(paste0(gene,".MP.hum.pdf"),p,height=5,width=5)
     })
  }



 

####plot FOXP2 and other genes expression###
rownames(umap_mp_m) <- umap_mp_m$sample_name
rownames(umap_mp_h) <- umap_mp_h$sample_name

select.genes = c("FOXP2","CDH8", "EPHA7","SUSD2","EPHA1","GRID2","EBF1","CRH", "ROBO2","SPON1")
select.genes = c("FOXP2", "ROBO3")
## max val=AD0404
g=plot_RD_gene(umap_mp_m, norm.dat=as.matrix(datlist[["macaque"]][,int_keepmac]), genes=select.genes, cex=2)
ggsave(paste0(st,".genes_mac_mpcluster.select.pdf"),gridExtra::marrangeGrob(g,nrow = 1, ncol=1),height=5, width=5)  


g=plot_RD_gene(umap_mp_h, norm.dat=as.matrix(datlist[["human"]][,int_keephum]), genes=select.genes, cex=2)
ggsave(paste0(st,".genes_hum_mpcluster.select.pdf"),gridExtra::marrangeGrob(g,nrow = 1, ncol=1),height=5, width=5)  



pdf("human_macaque_MP_plots.pdf",useDingbats=F,width=24,height=12)
par(mfrow=c(1,2))
plot(coords[int_keepmac,1:2],
     col=allcl_mac$cluster_color[match(int_keepmac,allcl_mac$sample_name)],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
plot(coords[int_keephum,1:2],
     col=allcl_hum$cluster_color[match(int_keephum,allcl_hum$sample_name)],pch=19,xlab="UMAP Axis 1",ylab="UMAP Axis 2")
for (gene in c("FOXP2","CDH8", "EPHA7","SUSD2","EPHA1","GRID2","EBF1","CRH", "ROBO2","SPON1")) {
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