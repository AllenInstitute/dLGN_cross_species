
########################
## Setup
########################

setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")

library(dplyr)
library(ggplot2)
library(cowplot)
#library(scrattch.hicat)
library(feather)
source("doublet.finder.r")
require(Matrix)
source("hicat_modules.r")
#require(scrattch.vis)
library(Seurat)
library(grid)
library(gridExtra)
library(scales)
library(reshape2)

load("collected_data_20191220.rda")

st=format(Sys.time(), "%Y%m%d_%H%M_")



########################
## Figure 5: Cross species plots
########################


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
#### convert to human ortholog
#rownames(datlist[[spec]] <- orthos$human_symbol[match(rownames(datlist[[spec]]),orthos$mouse_symbol)]
#keepgen <- intersect(orthos$human_symbol, rownames(datlist[[spec]])

datlist[[spec]]=datlist[[spec]][keepgen,]
datlist[[spec]]@x=datlist[[spec]]@x/rep.int(Matrix::colSums(datlist[[spec]]),diff(datlist[[spec]]@p))*10^6




#### setup metadata
cluster_palette <- read.csv("cluster_palette.csv")
roi_palette <- read.csv("roi_palette.csv",stringsAsFactors = FALSE)

load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Macaque"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Macaque"])]
allcl_mac$species="Macaque"
allcl_mac$species_color="#939598"
allcl_mac$dissection_roi <- roi_palette$dissection_roi[match(allcl_mac$roi, roi_palette$roi)]
allcl_mac$roi_color <- roi_palette$roi_color[match(allcl_mac$roi, roi_palette$roi)]

load("human_allcl_20200126.rda")
allcl_hum=allcl
allcl_hum$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Human"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Human"])]
allcl_hum$species="Human"
allcl_hum$species_color="#000000"
allcl_hum$dissection_roi <- roi_palette$dissection_roi[match(allcl_hum$roi, roi_palette$roi)]
allcl_hum$roi_color <- roi_palette$roi_color[match(allcl_hum$roi, roi_palette$roi)]

load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_color <- cluster_palette$cluster_color[cluster_palette$species_label == "Mouse"][match(allcl$cluster_id, cluster_palette$cluster_id[cluster_palette$species_label == "Mouse"])]
allcl_mus$species="Mouse"
allcl_mus$species_color="#FF73D0"
allcl_mus$dissection_roi <- roi_palette$dissection_roi[match(allcl_mus$roi, roi_palette$roi)]
allcl_mus$roi_color <- roi_palette$roi_color[match(allcl_mus$roi, roi_palette$roi)]

####integrate all neuronal cell types: note that the macaque donors need to be integrated separately###
#classkeep=c("GA","GL")
#keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$class %in% classkeep & cluster_palette$species_label=="Macaque"])
keepclasses=as.character(cluster_palette$cluster_label[cluster_palette$species_label=="Macaque"])

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

dimval <- 5
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
#p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T,split.by="species",cols=colvec)
#print(p1)
obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8)
obj.integrated <- FindClusters(object = obj.integrated, resolution = 0.2) #0.4
DimPlot(obj.integrated, reduction = "umap")

save(obj.integrated, file=paste0(st, "integrated_all.rda"))



##### umap plots original vals ####

umap_int_all <- obj.integrated@reductions$umap@cell.embeddings
seurat.meta<-obj.integrated@meta.data
seurat.meta$cluster.id_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_id)
seurat.meta$cluster.label_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_label)

clcols <- unique(seurat.meta[,c("cluster.id_all", "cluster.label_all", "cluster_color")])
clcol <- setNames(clcols$cluster_color, clcols$cluster.label_all)

## colored by original cluster
colnames(umap_int_all) = c("Lim1","Lim2")
g1= plot_RD_meta(umap_int_all, factor(seurat.meta[row.names(umap_int_all),] %>% pull(cluster.label_all),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=1)
g1 = g1+	ggtitle("cluster") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

g1[["data"]]$species <- seurat.meta$species
g2 <- g1 + facet_wrap(~species, ncol=1)


## colored by species
species.cols <- setNames(c("#939598","#000000","#FF73D0"), c("Macaque","Human","Mouse"))
g3 <- plot_RD_meta(umap_int_all, factor(seurat.meta[row.names(umap_int_all),] %>% pull(species),levels=names(species.cols)), meta.col=species.cols,alpha=0.5, cex=1)
g3 <- g3+	ggtitle("Species") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

g4 <- g3 + facet_wrap(~meta, ncol=1)


## colored by dissection roi
#roi_palette <- read.csv("roi_palette.csv",stringsAsFactors = FALSE)

rcol <- unique(seurat.meta[,c("dissection_roi", "roi_color")])
rcol <- setNames(rcol$roi_color, rcol$dissection_roi)

g5 <- plot_RD_meta(umap_int_all, factor(seurat.meta[row.names(umap_int_all),] %>% pull(dissection_roi),levels=names(rcol)), meta.col=rcol,alpha=0.5, cex=1)
g5 <- g5 +	ggtitle("Species") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
g5[["data"]]$species <- seurat.meta$species
g6 <- g5 + facet_wrap(~species, ncol=1)


ggsave(file=paste0(st,"species_all_umap.pdf"),gridExtra::marrangeGrob(list(g1,g3,g5),nrow = 1, ncol=1), height=5,width=4)
ggsave(file=paste0(st,"species_all_umap.split.pdf"),gridExtra::marrangeGrob(list(g2,g4,g6),nrow = 1, ncol=1), height=15,width=4)


###### Fig 5C/D ######

load("20200720_1316_integrated_all.rda")
umap_int_all <- as.data.frame(obj.integrated@reductions$umap@cell.embeddings)
colnames(umap_int_all) = c("Lim1","Lim2")

seurat.meta<-obj.integrated@meta.data

cluster_palette <- read.csv("cluster_palette.csv")
roi_palette <- read.csv("roi_palette.csv",stringsAsFactors = FALSE)

seurat.meta <- left_join(seurat.meta, select(cluster_palette, cluster_color, class,new_cluster_label))
rownames(seurat.meta) <- seurat.meta$sample_name

seurat.meta$cluster.id_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_id)
seurat.meta$cluster.label_all <- paste0(seurat.meta$species,"_", seurat.meta$new_cluster_label)
clcols <- unique(seurat.meta[,c("cluster.id_all", "cluster.label_all", "cluster_color")])
clcol <- setNames(clcols$cluster_color, clcols$cluster.label_all)

rcol <- unique(seurat.meta[,c("dissection_roi", "roi_color")])
rcol <- setNames(rcol$roi_color, rcol$dissection_roi)

#meta <- seurat.meta[seurat.meta$class == "GL",]
#umap <- umap_int_all[seurat.meta$sample_name,]

umap <- umap_int_all[umap_int_all$Lim1 > -4,] 
meta <- seurat.meta[seurat.meta$sample_name %in% rownames(umap),]

g7 <- plot_RD_meta(umap, factor(meta[row.names(umap),] %>% pull(dissection_roi),levels=names(rcol)), meta.col=rcol,alpha=0.5, cex=1)
g7 <- g7 +  ggtitle("Species") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme_void() + 
    theme(legend.position="none") 

ggsave(file=paste0(st,"all_species_GL_roi_umap.pdf"),g7, height=5,width=3, useDingbats=FALSE)



umap <- tibble::rownames_to_column(umap, var="sample_name")
umap <- umap %>% left_join(meta)

rcol <- unique(umap[,c("dissection_roi", "roi_color")])
rcol <- setNames(rcol$roi_color, rcol$dissection_roi)

library(ggridges)
x_axis <- ggplot(umap, aes(x = Lim1, y = dissection_roi, fill=roi_color)) +   
      geom_density_ridges(scale=10, aes(alpha=0.5), rel_min_height=0.01) + 
      scale_fill_identity() +
      #theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
      theme_void() +
      theme(legend.position="none") 

y_axis <- ggplot(umap, aes(x = Lim2, y = dissection_roi, fill=roi_color)) +   
      geom_density_ridges(scale=10, aes(alpha=0.5), rel_min_height=0.01) + 
      scale_fill_identity() +
      #theme_ridges(grid = FALSE, center_axis_labels = FALSE) +
      theme_void() +
      theme(legend.position="none") 
ggsave(file=paste0(st,"umap_ridges_5c.pdf"),gridExtra::marrangeGrob(list(x_axis, y_axis),nrow = 1, ncol=1), height=2,width=5)



###############
## identify genes correlated with umap-axes per species
###############

keepcells_mac=meta$sample_name[which(meta$species =="Macaque")]
keepcells_hum=meta$sample_name[which(meta$species =="Human")]
keepcells_mus=meta$sample_name[which(meta$species =="Mouse")]
coords=obj.integrated@reductions$umap@cell.embeddings[,1:2]
coords <- coords[meta$sample_name,]

require(Matrix)
axis_genelists=list()

### Macaque axis 1
mac_cells_axis1=coords[keepcells_mac,1]
genvals=t(cor(mac_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells_mac])),method="spearman"))
axis.genes <- as.data.frame(genvals)
colnames(axis.genes) <- "Mac.axis1"
axis_genelists[["macaque_axis1"]]=genvals


### Macaque axis 2
mac_cells_axis2=coords[keepcells_mac,2]
genvals=t(cor(mac_cells_axis2,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells_mac])),method="spearman"))
axis.genes$"Mac.axis2" <- genvals
axis_genelists[["macaque_axis2"]]=genvals


### Human axis 1
hum_cells_axis1=coords[keepcells_hum,1]
genvals=t(cor(hum_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells_hum])),method="spearman"))
axis.genes$"Hu.axis1" <- genvals
axis_genelists[["human_axis1"]]=genvals


### Human axis 2
hum_cells_axis2=coords[keepcells_hum,2]
genvals=t(cor(hum_cells_axis2,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells_hum])),method="spearman"))
axis.genes$"Hu.axis2" <- genvals
axis_genelists[["human_axis2"]]=genvals


### Mouse axis 1
mus_cells_axis1=coords[keepcells_mus,1]
genvals=t(cor(mus_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells_mus])),method="spearman"))
axis.genes$"Mo.axis1" <- genvals
axis_genelists[["mouse_axis1"]]=genvals

### Mouse axis 2
mus_cells_axis2=coords[keepcells_mus,2]
genvals=t(cor(mus_cells_axis2,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells_mus])),method="spearman"))
axis.genes$"Mo.axis2" <- genvals
axis_genelists[["mouse_axis2"]]=genvals

keepgenesplot <- axis.genes %>%
        tibble::rownames_to_column('gene') %>%
        filter_if(is.numeric, any_vars(. > 0.25 | .< -0.25))

genesplot <- melt(keepgenesplot)
genesplot$variable <- factor(genesplot$variable, levels=c("Mac.axis1","Hu.axis1","Mo.axis1", "Mac.axis2",  "Hu.axis2", "Mo.axis2"))
p.gene.cor <- ggplot(genesplot) + 
      geom_tile(aes(y=reorder(gene, value,median, order=TRUE), x=variable, fill =value)) +
      #scale_fill_gradient2(low = "midnightblue", mid = "white", high ="red", oob=squish)
      scale_fill_gradient(low = "#132B43", high ="#56B1F7", oob=squish) #"#8accff"

x <- ggplot_build(p.gene.cor)
gene.order.plot <- x[["layout"]][["panel_params"]][[1]][["y.sec"]][["limits"]]

genes.to.plot <- head(gene.order.plot, 30) 
genes.to.plot2 <- tail(gene.order.plot, 30)
extra <- c("foxp2","robo2", "camk2a")
genes.to.plot <- append(genes.to.plot,append(extra ,genes.to.plot2))
genesplot.s <- genesplot[genesplot$gene %in% genes.to.plot,]
p.gene.cor <- ggplot(genesplot.s) + 
      geom_tile(aes(y=reorder(gene, value,median, order=TRUE), x=variable, fill =value)) +
      scale_fill_gradient(low = "#132B43", high ="#56B1F7", oob=squish) #"#8accff"
ggsave(file=paste0(st, "kpml_integrated_heatmap_corr.axis.pdf"), p.gene.cor, width=5, height=8, useDingbats=FALSE)















###############
##### integrated comparison  Fig S4##### 
###############

seurat.meta <- obj.integrated@meta.data

species.bar <- seurat.meta %>% 
			select(species, seurat_clusters) %>% 
			group_by(seurat_clusters, species) %>%
			mutate(n = n()) %>%
			ungroup() %>%
			group_by(seurat_clusters) %>%
  			mutate(cluster_n = n(),
         			frac = n/cluster_n) %>%
			unique()

species.bar.plot <- ggplot(species.bar, aes(x=seurat_clusters,y=frac,fill=species)) + geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=c("Macaque"="#939598","Human"="#000000","Mouse"="#FF73D0")) 
ggsave(file=paste0(st,"integrated_bar.pdf"),species.bar.plot, height=4,width=6)

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
seurat.cols <- unique(seurat.meta[,c("seurat_clusters")])
g7 <- plot_RD_meta(umap_int_all, seurat.meta$seurat_clusters, alpha=0.5, cex=1)
g7 <- g7 +	ggtitle("seurat_clusters") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() #+ 
 		#theme(legend.position="none") 
ggsave(file=paste0(st,"integrated_umap.pdf"),g7, height=5,width=5)

g7[["data"]]$species <- seurat.meta$species
g8 <- g7 + facet_wrap(~species, ncol=1)
ggsave(file=paste0(st,"integrated_split_umap.pdf"),g8, height=15,width=4)


roi.bar <- seurat.meta %>% 
      select(roi, seurat_clusters) %>% 
      group_by(seurat_clusters, roi) %>%
      mutate(n = n()) %>%
      ungroup() %>%
      group_by(seurat_clusters) %>%
        mutate(cluster_n = n(),
              frac = n/cluster_n) %>%
      unique()
roi.cols <- setNames(roi_palette$roi_color, roi_palette$roi)
roi.bar.plot <- ggplot(roi.bar, aes(x=seurat_clusters,y=frac,fill=species)) + geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=roi.cols) 

### river plot ####

anno <- seurat.meta %>%                              
         select(sample_name, cluster.label_all,seurat_clusters, cluster_color) %>%
     arrange(cluster.label_all,seurat_clusters)
 

plot_order <- anno %>%
     group_by( cluster.label_all,seurat_clusters, cluster_color) %>%
     summarise(n_cells = n()) %>%
     ungroup() %>%
     group_by(cluster.label_all) %>%
     arrange(-n_cells) %>%
     filter(row_number() == 1) %>%
     ungroup() %>%
     arrange(cluster.label_all) %>%
     mutate(cluster_id = 1:n()) %>%
     select(cluster_id, cluster.label_all)
anno <- anno %>%
    #select(-cluster.label_all) %>%
    left_join(plot_order)

x=brewer.pal(n=9, name="Paired")
x <- setNames(x, c(0:8))
anno$seurat_color <- x[match(anno$seurat_clusters, names(x))]
anno$seurat_id <- anno$seurat_clusters
anno$seurat_label <- anno$seurat_clusters
colnames(anno)[2]<- "cluster_label"

nodes <- make_plot_nodes(make_group_nodes(anno, c("cluster","seurat")), pad = 0.2)
group_links <- make_group_links(anno, c("cluster","seurat"), nodes)
  
  
group1.cl.size =table(anno$cluster_id)
group2.cl.size =table(anno$seurat_id)

min.cells=4 
min.frac=0.1 
  group_links = group_links %>%
    mutate(group1_frac = n /group1.cl.size[as.character(group1_id)]) %>%
    mutate(group2_frac = n /group2.cl.size[as.character(group2_id)]) %>%
    filter(n >= min.cells) %>%
    filter(pmax(group1_frac, group2_frac)> min.frac)
  
  
  links <- make_plot_links(group_links, fill = "cluster")
  
  
  # Get nodes for labeling
  left_nodes <- nodes %>%
    filter(group == "cluster")
  right_nodes <- nodes %>%
    filter(group == "seurat")
  
  grob.1.x <- as.numeric(left_nodes$xpos[1])
  grob.2.x <- right_nodes$xpos[1]

library(grid)  
  grob1 <- grobTree(textGrob("Species cluster", x=0.29,  y=0.98, hjust=1,  gp=gpar(col="black", fontsize=8)))
  grob2 <- grobTree(textGrob("Integrated cluster", x=0.7,  y=0.98, hjust=0,  gp=gpar(col="black", fontsize=8)))
  
  
  river <- ggplot() +
    geom_rect(data = nodes,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = color),
              color = "#808080") +
    geom_ribbon(data = links,
                aes(x = x, ymax = y,
                    ymin = ymin,
                    group = link_id,
                    fill = fill),
                color = "#808080",
                alpha = 0.4) +
    geom_text(data = left_nodes,
              aes(x = xmin,
                  y = (ymin + ymax)/2,
                  label = name,
                  color = color),
              hjust = 1,
              nudge_x= -0.005,
              size = 8/6) +
    geom_text(data = right_nodes,
              aes(x = xmax,
                  y = (ymin + ymax)/2,
                  label = name,
                  color = color),
              hjust = 0,
              nudge_x= 0.005,
              size = 8/6) +
    annotation_custom(grob1) + 
    annotation_custom(grob2) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_y_reverse() +
    theme_void() +
    xlim(0.5, 2.5)

ggsave(file=paste0(st,"river_Species_integrated.pdf"),river, height=7,width=5)



### compare_annotate ####

ref.cl <- setNames(seurat.meta$cluster.label_all,rownames(seurat.meta))
cl <- setNames(seurat.meta$seurat_clusters,rownames(seurat.meta))
cldf <- unique(seurat.meta[,c("species", "cluster.id_all","cluster.label_all")])
rownames(cldf) <- cldf$cluster.label_all
cldf$cluster_label <- cldf$cluster.label_all

species <- unique(cldf$species)
for(i in species) {
	#ref.cl.spec <- droplevels(ref.cl[ref.cl %in% row.names(cldf)[cldf$species == i]])
	ref.cl.spec <- ref.cl[ref.cl %in% row.names(cldf)[cldf$species == i]]

	cl.df.spec <- cldf[cldf$species == i,]

	x<- compare_annotate(cl, ref.cl.spec, cl.df.spec, reorder=FALSE)
ggsave(file=paste0(st,i,"_annot_comp.pdf"),x$g, height=5,width=5)

}

x<- compare_annotate(cl, ref.cl, cldf, reorder=FALSE)
annot_plot <- x$g

library(stringr)
annot_plot[["data"]]$species <- str_extract(annot_plot[["data"]]$ref.cl.label , "[^_]+")

p <- annot_plot + coord_flip()
#p <- p + facet_wrap(~species)

levels <- c("Macaque_GABA1","Macaque_GABA2","Macaque_GABA3","Macaque_GABA4", "Macaque_K1",     "Macaque_K2", "Macaque_MP1", "Macaque_MP2","Macaque_O1","Macaque_O2","Macaque_O3","Human_GABA1" ,"Human_GABA2" , "Human_GABA3", "Human_GABA4","Human_K1" , "Human_K2" , "Human_MP1" ,  "Mouse_Chrna1",  "Mouse_Chrna2",  "Mouse_Chrna3", "Mouse_GABA1" ,"Mouse_GABA4","Mouse_GABA5", "Mouse_GABA3",   "Mouse_GABA6","Mouse_GABA7", "Mouse_LGD", "Mouse_LP" ,"Mouse_LGV" )
annot_plot[["data"]]$ref.cl.label <- ordered(annot_plot[["data"]]$ref.cl.label, levels=levels)
p <- annot_plot + coord_flip()
ggsave(file=paste0(st,"all_annot_comp.pdf"),p, height=4,width=5)



##Cross correlation between integrated clusters and Mouse only clusters
load("mouse_clusters_subclusters.rda")
rownames(m.int@assays[["RNA"]]@data)=tolower(rownames(m.int@assays[["RNA"]]@data))

genes<-intersect(rownames(obj.integrated@assays$integrated@data),rownames(m.int@assays[["RNA"]]@data))

mouse_exp <- m.int@assays[["RNA"]]@data[genes,]
integrated_exp<-obj.integrated@assays$RNA@data[genes,]

mouse_exp<-t(scale(t(as.matrix(mouse_exp))))
integrated_exp<-t(scale(t(as.matrix(integrated_exp))))

integrated.cl <- unique(seurat.meta$seurat_clusters)
mo.cl <- unique(seurat.meta$cluster.label_all[seurat.meta$species == "Mouse"])

corr<- data.frame(matrix(NA, ncol=length(mo.cl),nrow=length(integrated.cl)))
rownames(corr)<-factor(integrated.cl,levels=integrated.cl)
colnames(corr)<-factor(mo.cl,levels=mo.cl)

for (i in mo.cl){
    for (j in integrated.cl){
    	select.cells.i <- seurat.meta$sample_name[seurat.meta$cluster.label_all == i]
    	select.cells.j <- seurat.meta$sample_name[seurat.meta$seurat_clusters == j]
    corr[j,i]<-cor(rowMeans(as.matrix(mouse_exp[,select.cells.i])),rowMeans(as.matrix(integrated_exp[,select.cells.j])))
}}


p_value<- data.frame(matrix(NA, ncol=length(mo.cl),nrow=length(integrated.cl)))
rownames(p_value)<-factor(integrated.cl,levels=integrated.cl)
colnames(p_value)<-factor(mo.cl,levels=mo.cl)

for (i in mo.cl){
    for (j in integrated.cl){
    	select.cells.i <- seurat.meta$sample_name[seurat.meta$cluster.label_all == i]
    	select.cells.j <- seurat.meta$sample_name[seurat.meta$seurat_clusters == j]
   p_value[j,i]<-cor.test(rowMeans(as.matrix(mouse_exp[,select.cells.i])),rowMeans(as.matrix(integrated_exp[,select.cells.j])))$p.value
}}

test.m <- melt(as.matrix(corr))
test.p <- melt(as.matrix(p_value))


for(i in 1:dim(test.m)[1]){
    if (test.p$value[i]>0.05){
        test.m$value[i]<-0
    }else{
        test.m$value[i]<-test.m$value[i]
    }
}

test.m$X1 <- as.character(test.m$X1)

library("scales")
p.mo <- ggplot(test.m, aes(X2, X1)) + 
			geom_tile(aes(fill =value),colour = "white")+ 
			scale_fill_continuous(limits=c(0, 0.7), breaks=seq(0,0.75,by=0.25),low = "white",high = "midnightblue", oob=squish)
#scale_x_discrete(breaks=c(0,5994,7712,9820,11428,11674,11836,11878))
p.mo<-p.mo+
		ylab("integrated")+
		xlab("Mouse")+
		theme(axis.text.x=element_text(size=10, angle=90, hjust=1),
				axis.title=element_text(size=10,face="bold"),
				plot.title = element_text(size=10),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black", size = 0.1),
				axis.ticks.y = element_blank())+
		ggtitle("correlation between integrated and mouse clusters")
print(p.mo)
ggsave(file=paste0(st,"mouse_cl_comp_cor.pdf"),p.mo, height=5,width=5.5)


###############
### Magno/Parvo/Konio, Pulv and LGd/LP
###############

#keepclasses=unique(as.character(grep("K|MP|LP|LGD",cluster_palette$cluster_label,val=T)))
#keepclasses=unique(as.character(grep("K|MP|LP|LGD|O",cluster_palette$cluster_label,val=T)))
keepclasses=unique(as.character(grep("K|MP|LGD",cluster_palette$cluster_label,val=T)))

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


dimval <- 5
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
#p1=DimPlot(obj.integrated, reduction = "umap", group.by = "cluster_label",label=T,split.by="species",cols=colvec)
#print(p1)
obj.integrated <- FindNeighbors(object = obj.integrated, k.param = 8)
obj.integrated <- FindClusters(object = obj.integrated, resolution = 0.2) #0.4
DimPlot(obj.integrated, reduction = "umap")
save(obj.integrated,file=paste0(st,"integrated_kmpl.rda"))

### test to determine if core/shell can be separated
uniqr=unique(obj.integrated@meta.data$roi)
uniqr=uniqr[order(uniqr)]
colvec=obj.integrated@meta.data$roi_color[match(uniqr,obj.integrated@meta.data$roi)]

colvec <- c("")

pdf(paste0(st,"hu_ma_mo_kpml_integrated_umap_dimensions_roi.pdf"),width=15)
for (dimval in c(5,10,15,20,25,30)) {
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval)
  p1=DimPlot(obj.integrated, reduction = "umap", group.by = "dissection_roi",label=T,split.by="species")
  print(p1)
}
dev.off()

obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:10)

##### umap plots original vals ####
umap_int_kmpl <- obj.integrated@reductions$umap@cell.embeddings
seurat.meta<-obj.integrated@meta.data
seurat.meta$cluster.id_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_id)
seurat.meta$cluster.label_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_label)

clcols <- unique(seurat.meta[,c("cluster.id_all", "cluster.label_all", "cluster_color")])
clcol <- setNames(clcols$cluster_color, clcols$cluster.label_all)

## colored by original cluster
colnames(umap_int_kmpl) = c("Lim1","Lim2")
g1= plot_RD_meta(umap_int_kmpl, factor(seurat.meta[row.names(umap_int_kmpl),] %>% pull(cluster.label_all),levels=names(clcol)), meta.col=clcol,alpha=0.5, cex=1)
g1 = g1+	ggtitle("cluster") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

g1[["data"]]$species <- seurat.meta$species
g2 <- g1 + facet_wrap(~species, ncol=1)

#library(ggplus)
#facet_multiple(plot = p, facets = 'id', ncol = 2, nrow = 2)

## colored by species
species.cols <- setNames(c("#939598","#000000","#FF73D0"), c("Macaque","Human","Mouse"))
g3 <- plot_RD_meta(umap_int_kmpl, factor(seurat.meta[row.names(umap_int_kmpl),] %>% pull(species),levels=names(species.cols)), meta.col=species.cols,alpha=0.5, cex=1)
g3 <- g3+	ggtitle("Species") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 

g4 <- g3 + facet_wrap(~meta, ncol=1)


## colored by dissection roi
#roi_palette <- read.csv("roi_palette.csv",stringsAsFactors = FALSE)

rcol <- unique(seurat.meta[,c("dissection_roi", "roi_color")])
rcol <- setNames(rcol$roi_color, rcol$dissection_roi)

g5 <- plot_RD_meta(umap_int_kmpl, factor(seurat.meta[row.names(umap_int_kmpl),] %>% pull(dissection_roi),levels=names(rcol)), meta.col=rcol,alpha=0.5, cex=1)
g5 <- g5 +	ggtitle("Species") +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
 		theme_void() + 
 		theme(legend.position="none") 
g5[["data"]]$species <- seurat.meta$species
g6 <- g5 + facet_wrap(~species, ncol=1)


p <- ggExtra::ggMarginal(g5, type = "density",groupColour = TRUE, groupFill = TRUE, size=5)
ggsave(file=paste0(st,"roi_marg_test.pdf"),p, height=6,width=6)

ggsave(file=paste0(st,"species_kmpl_umap.pdf"),gridExtra::marrangeGrob(list(g1,g3,g5),nrow = 1, ncol=1), height=5,width=5)
ggsave(file=paste0(st,"species_kmpl_umap.split.pdf"),gridExtra::marrangeGrob(list(g2,g4,g6),nrow = 1, ncol=1), height=15,width=5)

library(ggridges)
umap <- tibble::rownames_to_column(as.data.frame(umap_int_kmpl), var="sample_id")
meta <- tibble::rownames_to_column(seurat.meta, var="sample_id")
umap <- left_join(umap, meta)

exclud <- c("Pulv", "LP", "LGv")
umap.s <-  umap[!(umap$dissection_roi %in% exclud),]

x_axis <- ggplot(umap.s, aes(x = Lim1, y = dissection_roi, fill=roi_color)) +   
      geom_density_ridges(scale=10, aes(alpha=0.5), rel_min_height=0.01) + 
      scale_fill_identity() +
      #theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
      theme_void() +
      theme(legend.position="none") 

y_axis <- ggplot(umap.s, aes(x = Lim2, y = dissection_roi, fill=roi_color)) +   
      geom_density_ridges(scale=10, aes(alpha=0.5), rel_min_height=0.01) + 
      scale_fill_identity() +
      #theme_ridges(grid = FALSE, center_axis_labels = FALSE) +
      theme_void() +
      theme(legend.position="none") 
ggsave(file=paste0(st,"umap_ridges_s4.pdf"),gridExtra::marrangeGrob(list(x_axis, y_axis),nrow = 1, ncol=1), height=2,width=5)



### marker gene expression

genes <- c("foxp2","cntn5","kcnb2","grm7","ebf1","necab1", "pvalb", "greb1l","btnl9","dgkb" )



select.genes = c("foxp2","cntn5","kcnb2","grm7","ebf1","necab1", "pvalb", "greb1l","btnl9","dgkb" )
## max val=AD0404
g=plot_RD_gene(umap_int_kmpl, norm.dat=as.matrix(obj.integrated@assays$RNA@data), genes=select.genes, cex=2)
ggsave(paste0(st,".genes_integrated.kmpl.select.pdf"),gridExtra::marrangeGrob(g,nrow = 1, ncol=1),height=5, width=5)  


g=plot_RD_gene(umap_mp_h, norm.dat=as.matrix(datlist[["human"]][,int_keephum]), genes=select.genes, cex=2)
ggsave(paste0(st,".genes_hum_mpcluster.select.pdf"),gridExtra::marrangeGrob(g,nrow = 1, ncol=1),height=5, width=5)  

pacman::p_unlock(lib.loc = pacman::p_path())
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')

##### Find markers that are conserved between species
#load("20200413_1635_integrated_kmpl.rda")
#load("20200720_1316_integrated_kmpl.rda")
seurat.meta<-obj.integrated@meta.data
seurat.meta$cluster.id_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_id)
seurat.meta$cluster.label_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_label)

#trace(FindConservedMarkers, edit=TRUE) #to add min.cell.group=1 to prevent error
levs <- levels(seurat.meta$seurat_clusters)
markers <- list()
for(i in levs){
markers[i] <- FindConservedMarkers(obj.integrated, ident.1 = i,grouping.var = "species", verbose = TRUE)}

FindMarkers(obj.integrated, ident.1 = 1, ident.2 = 0, verbose = TRUE)


Cell_type<-factor(c("Macaque_K1",     "Macaque_K2", "Macaque_MP1", "Macaque_MP2","Human_K1" , "Human_K2" , "Human_MP1" , "Mouse_LGD", "Mouse_LP" ), levels=c("Macaque_K1",     "Macaque_K2", "Macaque_MP1", "Macaque_MP2","Human_K1" , "Human_K2" , "Human_MP1" , "Mouse_LGD", "Mouse_LP" ))

gene_list<-factor(rev(c("foxp2","cntn5","kcnb2","grm7","ebf1","necab1", "pvalb", "greb1l","btnl9","dgkb" )),levels=rev(c("foxp2","cntn5","kcnb2","grm7","ebf1","necab1", "pvalb", "greb1l","btnl9","dgkb")))

rownames(obj.integrated@meta.data)# the matrix has info of pct, avg, gene, cluster

Cell_number<- data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)

for (i in 1:length(gene_list)){
	L<-length(Cell_type)

	Cell_number_t<- data.frame("cluster" =Cell_type, "gene"=(rep(gene_list[i],L))) # do not use c if the gene is factorizsed
#used normalized uncorrected data
	for (p in 1:length(Cell_type)){
		cl.lab <- Cell_number_t$cluster[p]
		select.cells.i <- seurat.meta$sample_name[seurat.meta$cluster.label_all == cl.lab]
		Cell_number_t$pct[p]<-100*sum(obj.integrated@assays$RNA@data[gene_list[i],select.cells.i]>0)/length(select.cells.i)
		Cell_number_t$avg[p]<-(mean(obj.integrated@assays$RNA@data[gene_list[i],select.cells.i])-mean(obj.integrated@assays$RNA@data[gene_list[i],]))/sd(obj.integrated@assays$RNA@data[gene_list[i],])
	#for avg only consider the expressed cell
	 #t<-LHb.integrated@assays$RNA@data[gene_list[i],eval(parse(text=paste(Cell_type[p],"_barcode",sep="")))]>0
	# Cell_number_t$avg[p]<-mean(LHb.integrated@assays$RNA@data[gene_list[i],t])/sd(LHb.integrated@assays$RNA@data[gene_list[i],eval(parse(text=paste(Cell_type[p],"_barcode",sep="")))])    
    	}
	Cell_number<-rbind(Cell_number_t,Cell_number)}

as.factor(Cell_number$cluster)

ggplot(Cell_number, aes(gene, cluster)) + geom_point(aes(size = pct, colour=avg)) +  scale_x_discrete(limits = rev(levels(Cell_number$gene)))+scale_y_discrete(limits =rev(levels(Cell_number$cluster)))+
scale_color_gradient(low = "white", high = "darkblue",limits = c(-1.5,1.5),oob=squish) + 
geom_point(aes(size = pct), pch=21,, lwd=0,stroke=0)+ scale_size_continuous(range = c(0,7),limits=c(0,100),breaks=seq(0,100,25))+
theme(axis.title.y=element_text(size=15),axis.text.y=element_text(size=10,colour = "black"),axis.text.x=element_text(size=8,angle = 50, hjust =1,colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("% cells")
ggsave(file=".pdf",height=4, width=16 , paper = "letter")





###############
## identify genes correlated with umap-axes per species
###############

st=format(Sys.time(), "%Y%m%d_%H%M_")

load("collected_data_20191220.rda")
#load("20200414_1315_integrated_all.rda")
#seurat.meta<-obj.integrated@meta.data
#seurat.meta$cluster.id_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_id)
#seurat.meta$cluster.label_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_label)

#load("20200413_1635_integrated_kmpl.rda")
load("macaque_allcl.rda")
allcl_mac=allcl
allcl_mac$cluster_label=paste0(allcl_mac$cluster_label,"_M")
load("human_allcl.rda")
allcl_hum=allcl
allcl_hum$cluster_label=paste0(allcl_hum$cluster_label,"_H")
load("mouse_allcl.rda")
allcl_mus=allcl
allcl_mus$cluster_label=paste0(allcl_mus$cluster_label,"_Ms")




grepchar=""
keepcells_mac1=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q17"))]
keepcells_mac2=allcl_mac$sample_name[intersect(grep(grepchar,allcl_mac$cluster_label),which(allcl_mac$donor=="Q18"))]
keepcells_hum=allcl_hum$sample_name[grep(grepchar,allcl_hum$cluster_label)]
keepcells_mus=allcl_mus$sample_name[grep(grepchar,allcl_mus$cluster_label)]
coords=obj.integrated@reductions$umap@cell.embeddings[,1:2]

require(Matrix)
axis_genelists=list()

### Macaque axis 1
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP|K",allcl_mac$cluster_label)])]
mac_cells_axis1=coords[keepcells,1]
#genvals=t(cor(mac_cells_axis1,as.matrix(t(datlist[["macaque"]][,keepcells])),method="spearman"))
genvals=t(cor(mac_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells])),method="spearman"))
axis.genes <- as.data.frame(genvals)
colnames(axis.genes) <- "Mac.axis1"
axis_genelists[["macaque_axis1"]]=genvals


### Macaque axis 2
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mac$sample_name[grep("MP",allcl_mac$cluster_label)] & coords[,2]>-2)]
mac_cells_axis2=coords[keepcells,2]
#genvals=t(cor(mac_cells_axis1,as.matrix(t(datlist[["macaque"]][,keepcells])),method="spearman"))
genvals=t(cor(mac_cells_axis2,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells])),method="spearman"))
axis.genes$"Mac.axis2" <- genvals
axis_genelists[["macaque_axis2"]]=genvals


### Human axis 1
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP|K",allcl_hum$cluster_label)])]
hum_cells_axis1=coords[keepcells,1]
genvals=t(cor(hum_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells])),method="spearman"))
axis.genes$"Hu.axis1" <- genvals
axis_genelists[["human_axis1"]]=genvals


### Human axis 2
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_hum$sample_name[grep("MP",allcl_hum$cluster_label)] & coords[,2]>-2)]
hum_cells_axis1=coords[keepcells,2]
genvals=t(cor(hum_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells])),method="spearman"))
axis.genes$"Hu.axis2" <- genvals
axis_genelists[["human_axis2"]]=genvals


### Mouse axis 1
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)] & coords[,2]>-2)]
mus_cells_axis1=coords[keepcells,1]
genvals=t(cor(mus_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells])),method="spearman"))
axis.genes$"Mo.axis1" <- genvals
axis_genelists[["mouse_axis1"]]=genvals

### Mouse axis 2
keepcells=rownames(coords)[which(rownames(coords) %in% allcl_mus$sample_name[grep("LGD",allcl_mus$cluster_label)])]
mus_cells_axis1=coords[keepcells,2]
genvals=t(cor(mus_cells_axis1,as.matrix(t(obj.integrated@assays[["integrated"]]@data[,keepcells])),method="spearman"))
axis.genes$"Mo.axis2" <- genvals
axis_genelists[["mouse_axis2"]]=genvals

keepgenesplot <- axis.genes %>%
				tibble::rownames_to_column('gene') %>%
				filter_if(is.numeric, any_vars(. > 0.25 | .< -0.25))

genesplot <- melt(keepgenesplot)
genesplot$variable <- factor(genesplot$variable, levels=c("Mac.axis1","Hu.axis1","Mo.axis1", "Mac.axis2",  "Hu.axis2", "Mo.axis2"))
p.gene.cor <- ggplot(genesplot) + 
			geom_tile(aes(y=reorder(gene, value,median, order=TRUE), x=variable, fill =value)) +
			#scale_fill_gradient2(low = "midnightblue", mid = "white", high ="red", oob=squish)
      scale_fill_gradient(low = "#132B43", high ="#56B1F7", oob=squish) #"#8accff"

x <- ggplot_build(p.gene.cor)
gene.order.plot <- x[["layout"]][["panel_params"]][[1]][["y.sec"]][["limits"]]

genes.to.plot <- head(gene.order.plot, 30) 
genes.to.plot2 <- tail(gene.order.plot, 30)
extra <- c("foxp2","robo2", "camk2a")
genes.to.plot <- append(genes.to.plot,append(extra ,genes.to.plot2))
genesplot.s <- genesplot[genesplot$gene %in% genes.to.plot,]
p.gene.cor <- ggplot(genesplot.s) + 
			geom_tile(aes(y=reorder(gene, value,median, order=TRUE), x=variable, fill =value)) +
			scale_fill_gradient(low = "#132B43", high ="#56B1F7", oob=squish) #"#8accff"
ggsave(file=paste0(st, "kpml_integrated_heatmap_corr.axis.pdf"), p.gene.cor, width=5, height=8, useDingbats=FALSE)


#scale_x_discrete(breaks=c(0,5994,7712,9820,11428,11674,11836,11878))
p.humo<-p.humo+
		ylab("Human")+
		xlab("Mouse")+
		theme(axis.text.x=element_text(size=10, angle=90, hjust=1),
				axis.title=element_text(size=10,face="bold"),
				plot.title = element_text(size=10),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black", size = 0.1),
				axis.ticks.y = element_blank())+
		ggtitle("correlation between human and mouse clusters")



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



########################
#### 1d umap of exc
########################

obj.integrated.1d <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dimval,n.components=1)
umap.1d<-obj.integrated.1d@reductions[["umap"]]@cell.embeddings
anno <- tibble::rownames_to_column(as.data.frame(umap.1d),var="sample_name")
anno <- left_join(anno, seurat.meta)


g13 <- ggplot(anno, aes(y=species, x=UMAP_1)) + 
	geom_jitter(aes(color = dissection_roi),
    				stat = "identity", 
    				position = position_jitter(0.2)) +
	scale_color_manual(values=rcol)

g14 <- ggplot(anno, aes(y=species, x=UMAP_1)) + 
	geom_tile(aes(color = dissection_roi)) +
	scale_color_manual(values=rcol)

g15 <- ggplot(anno, aes(y=species, x=UMAP_1)) + 
	geom_point(aes(color = dissection_roi), shape=124, size=8) +
	scale_color_manual(values=rcol)+
	theme_void() + 
 	theme(legend.position="none") +
 	theme(axis.text.y=element_text())

g16 <- ggplot(anno, aes(y=dissection_roi, x=UMAP_1,color = dissection_roi)) + 
  #geom_jitter(stat = "identity", position = position_jitter(0.2)) +
  geom_quasirandom(groupOnX=FALSE)+
  scale_color_manual(values=rcol) +
  xlab('UMAP 1')+
  ylab('Dissection roi')+
  facet_wrap(~species, ncol=1)+
  theme_bw()

library(ggridges)
g17 <-	ggplot(anno) +
    geom_density_ridges(aes(y = species, x = UMAP_1,
                            group = interaction(species, dissection_roi),
                            fill = dissection_roi), alpha=0.4) +
    scale_fill_manual(values=rcol)

ggsave(file=paste0(st,"species_kmpl_1dumap.pdf"), gridExtra::marrangeGrob(list(g13,g14,g15,g16,g17),nrow = 1, ncol=1), height=5,width=5, useDingbats=FALSE)


########################
## Gene specificty clusters between species
########################

setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")

library(dplyr)
library(ggplot2)
library(cowplot)
#library(scrattch.hicat)
library(feather)
source("doublet.finder.r")
require(Matrix)
source("hicat_modules.r")
#require(scrattch.vis)
library(Seurat)
library(grid)
library(gridExtra)
library(scales)
library(reshape2)

load("collected_data_20191220.rda")

st=format(Sys.time(), "%Y%m%d_%H%M_")

#load("20200414_1315_integrated_all.rda")
load("20200720_1316_integrated_all.rda")


seurat.meta<-obj.integrated@meta.data
seurat.meta$cluster.id_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_id)
seurat.meta$cluster.label_all <- paste0(seurat.meta$species,"_", seurat.meta$cluster_label)

########################
#### Cross correlation of cluster genes 
#### data prep
########################

load("human_clusters_subclusters_20200126.rda")
hum.int <- m.int
rownames(hum.int@assays[["RNA"]]@data)=tolower(rownames(hum.int@assays[["RNA"]]@data))

load("mouse_clusters_subclusters.rda")
mo.int <- m.int
rownames(mo.int@assays[["RNA"]]@data)=tolower(rownames(mo.int@assays[["RNA"]]@data))

load("macaque_clusters_subclusters.rda")
ma.int <- m.int
rownames(ma.int@assays[["RNA"]]@data)=tolower(rownames(ma.int@assays[["RNA"]]@data))

genes<-Reduce(intersect, list(rownames(mo.int@assays[["RNA"]]@data),rownames(ma.int@assays[["RNA"]]@data),rownames(hum.int@assays[["RNA"]]@data)))

mac_exp <- ma.int@assays[["RNA"]]@data[genes,]
mouse_exp <- mo.int@assays[["RNA"]]@data[genes,]
human_exp<-hum.int@assays[["RNA"]]@data[genes,]

mac_exp <- t(scale(t(as.matrix(mac_exp))))
mouse_exp<-t(scale(t(as.matrix(mouse_exp))))
human_exp<-t(scale(t(as.matrix(human_exp))))

#mac.cl <- unique(seurat.meta$cluster.label_all[seurat.meta$species == "Macaque"])
#mo.cl <- unique(seurat.meta$cluster.label_all[seurat.meta$species == "Mouse"])
#human.cl <- unique(seurat.meta$cluster.label_all[seurat.meta$species == "Human"])
mac.cl <- c("Macaque_GABA1","Macaque_GABA2","Macaque_GABA3","Macaque_GABA4", "Macaque_K1", "Macaque_K2", "Macaque_MP1", "Macaque_MP2", "Macaque_O1","Macaque_O2","Macaque_O3")
mo.cl <- c("Mouse_Chrna1",  "Mouse_Chrna2",  "Mouse_Chrna3", "Mouse_GABA1" ,"Mouse_GABA4","Mouse_GABA5", "Mouse_GABA3",   "Mouse_GABA6","Mouse_GABA7", "Mouse_LGD", "Mouse_LP" ,"Mouse_LGV" )
human.cl <- c("Human_GABA1" ,"Human_GABA2" , "Human_GABA3", "Human_GABA4","Human_K1" , "Human_K2" , "Human_MP1" )


##### human vs mouse #####
corr<- data.frame(matrix(NA, ncol=length(mo.cl),nrow=length(human.cl)))
rownames(corr)<-factor(human.cl,levels=human.cl)
colnames(corr)<-factor(mo.cl,levels=mo.cl)

p_value <- corr

for (i in mo.cl){
    for (j in human.cl){
    	select.cells.i <- seurat.meta$sample_name[seurat.meta$cluster.label_all == i]
    	select.cells.j <- seurat.meta$sample_name[seurat.meta$cluster.label_all == j]
    	x=rowMeans(as.matrix(mouse_exp[,select.cells.i]))
    	#x <- x[!is.na(x)]
    	y=rowMeans(as.matrix(human_exp[,select.cells.j]))
    	#y <- y[!is.na(y)]
    	to.select <- intersect(names(x[!is.na(x)]), names(y[!is.na(y)]))
    	x <- x[to.select]
    	y <- y[to.select]
    	
    	corr[j,i]<-cor(x,y) 

    	p_value[j,i]<-cor.test(x,y)$p.value
    }
}

test.m <- melt(as.matrix(corr))
test.p <- melt(as.matrix(p_value))

for(i in 1:dim(test.m)[1]){
    if (test.p$value[i]>0.05){
        test.m$value[i]<- 0
    }else{
        test.m$value[i]<-test.m$value[i]
    }
}

p.humo <- ggplot(test.m, aes(Var2, Var1)) + 
			geom_tile(aes(fill =value),colour = "white")+ 
			scale_fill_continuous(limits=c(0, 0.3), breaks=seq(0,0.3,by=0.01,low = "white",high = "midnightblue", oob=squish)
			#scale_fill_gradient2(low = "ivory4", mid = "white", high = "midnightblue", oob=squish)
#scale_x_discrete(breaks=c(0,5994,7712,9820,11428,11674,11836,11878))
p.humo<-p.humo+
		ylab("Human")+
		xlab("Mouse")+
		theme(axis.text.x=element_text(size=10, angle=90, hjust=1),
				axis.title=element_text(size=10,face="bold"),
				plot.title = element_text(size=10),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black", size = 0.1),
				axis.ticks.y = element_blank())+
		ggtitle("correlation between human and mouse clusters")
print(p.humo)


##### human vs macaque #####

corr<- data.frame(matrix(NA, ncol=length(mac.cl),nrow=length(human.cl)))
rownames(corr)<-factor(human.cl,levels=human.cl)
colnames(corr)<-factor(mac.cl,levels=mac.cl)

p_value <- corr

for (i in mac.cl){
    for (j in human.cl){
    	select.cells.i <- seurat.meta$sample_name[seurat.meta$cluster.label_all == i]
    	select.cells.j <- seurat.meta$sample_name[seurat.meta$cluster.label_all == j]
    	x=rowMeans(as.matrix(mac_exp[,select.cells.i]))
    	#x <- x[!is.na(x)]
    	y=rowMeans(as.matrix(human_exp[,select.cells.j]))
    	#y <- y[!is.na(y)]
    	to.select <- intersect(names(x[!is.na(x)]), names(y[!is.na(y)]))
    	x <- x[to.select]
    	y <- y[to.select]
    	
    	corr[j,i]<-cor(x,y) 

    	p_value[j,i]<-cor.test(x,y)$p.value
    }
}

test.m <- melt(as.matrix(corr))
test.p <- melt(as.matrix(p_value))

for(i in 1:dim(test.m)[1]){
    if (test.p$value[i]>0.05){
        test.m$value[i]<-0
    }else{
        test.m$value[i]<-test.m$value[i]
    }
}

p.huma <- ggplot(test.m, aes(Var2, Var1)) + 
			geom_tile(aes(fill =value),colour = "white")+ 
			scale_fill_continuous(limits=c(0, 0.6), breaks=seq(0,0.6,by=0.2),low = "white",high = "midnightblue", oob=squish)
			#scale_fill_gradient2(low = "ivory4", mid = "white", high = "midnightblue", oob=squish)
p.huma <-p.huma+
		ylab("Human")+
		xlab("Macaque")+
		theme(axis.text.x=element_text(size=10, angle=90, hjust=1),
				axis.title=element_text(size=10,face="bold"),
				plot.title = element_text(size=10),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black", size = 0.1),
				axis.ticks.y = element_blank())+
		ggtitle("correlation between human and macaque clusters")
print(p.huma)



##### macaque vs mouse #####

corr<- data.frame(matrix(NA, ncol=length(mac.cl),nrow=length(mo.cl)))
rownames(corr)<-factor(mo.cl,levels=mo.cl)
colnames(corr)<-factor(mac.cl,levels=mac.cl)

p_value <- corr

for (i in mac.cl){
    for (j in mo.cl){
    	select.cells.i <- seurat.meta$sample_name[seurat.meta$cluster.label_all == i]
    	select.cells.j <- seurat.meta$sample_name[seurat.meta$cluster.label_all == j]
    	x=rowMeans(as.matrix(mac_exp[,select.cells.i]))
    	#x <- x[!is.na(x)]
    	y=rowMeans(as.matrix(mouse_exp[,select.cells.j]))
    	#y <- y[!is.na(y)]
    	to.select <- intersect(names(x[!is.na(x)]), names(y[!is.na(y)]))
    	x <- x[to.select]
    	y <- y[to.select]
    	
    	corr[j,i]<-cor(x,y) 

    	p_value[j,i]<-cor.test(x,y)$p.value
    }
}

test.m <- melt(as.matrix(corr))
test.p <- melt(as.matrix(p_value))

for(i in 1:dim(test.m)[1]){
    if (test.p$value[i]>0.05){
        test.m$value[i]<-0
    }else{
        test.m$value[i]<-test.m$value[i]
    }
}

p.moma <- ggplot(test.m, aes(Var1, Var2)) + 
			geom_tile(aes(fill =value),colour = "white")+ 
			scale_fill_continuous(limits=c(0, 0.3), breaks=seq(0,0.3,by=0.1),low = "white",high = "midnightblue", oob=squish)
			#scale_fill_gradient2(low = "ivory4", mid = "white", high = "midnightblue", oob=squish)
p.moma <-p.moma+
		xlab("Mouse")+
		ylab("Macaque")+
		theme(axis.text.x=element_text(size=10, angle=90, hjust=1),
				axis.title=element_text(size=10,face="bold"),
				plot.title = element_text(size=10),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black", size = 0.1),
				axis.ticks.y = element_blank())+
		ggtitle("correlation between mouse and macaque clusters")
print(p.moma)

ggsave(file=paste0(st,"cluster_cor_species.pdf"),gridExtra::marrangeGrob(list(p.huma, p.humo, p.moma),nrow = 1, ncol=1), height=5,width=5)








