#https://cole-trapnell-lab.github.io/monocle3/docs/differential/

########################
## Setup for Monocle3
########################

setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")

library(dplyr)
library(Seurat)
library(VGAM)
library(monocle3)

load("collected_data_20191220.rda")

st=format(Sys.time(), "%Y%m%d_%H%M_")

load("mouse_allcl.rda")
load("mouse_plotobjects.rda")


########################
## ALL
########################

data <- datlist$mouse
data <- as(as.matrix(plot_objects[["all"]]@assays[["RNA"]]@data), 'sparseMatrix')

gene_metadata <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

pd <- plot_objects[["all"]]@meta.data



##create cds object
cds <- new_cell_data_set(data,
                         cell_metadata =pd,
                         gene_metadata = gene_metadata)

#View data
pData(cds)
fData(cds)

##reduce dimensions or import if already peformed
#cds <- preprocess_cds(cds, num_dim = 100)
#cds <- reduce_dimension(cds)
umap2d <- plot_objects[["all"]]@reductions[["umap"]]@cell.embeddings
reducedDims(cds)$UMAP = plot_objects[["all"]]@reductions[["umap"]]@cell.embeddings

##clustering
cds <- cluster_cells(cds)

# or use Seurat cluster result
cl <- setNames(plot_objects[["all"]]@meta.data$cluster_id, plot_objects[["all"]]@meta.data$sample_name )
cds@clusters$UMAP$clusters = cl

cds <- learn_graph(cds)
cds <- order_cells(cds,
					reduction_method = "UMAP",
                    root_pr_nodes=NULL,
                    root_cells=NULL,
                    verbose = TRUE)






########################
## dLGN
########################

pd <- plot_objects[["LGD"]]@meta.data
colnames(pd)[4] <- "exp_comp_name"
pd <- pd[grep("LGN", pd$roi),]


data <- as(as.matrix(plot_objects[["LGD"]]@assays[["RNA"]]@data), 'sparseMatrix')
data <- data[,rownames(pd)]

gene_metadata <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))


##create cds object
cds <- new_cell_data_set(data,
                         cell_metadata =pd,
                         gene_metadata = gene_metadata)

#View data
pData(cds)
fData(cds)

##reduce dimensions or import if already peformed
#cds <- preprocess_cds(cds, num_dim = 100)
#cds <- reduce_dimension(cds)
umap.2d <- plot_objects[["LGD"]]@reductions[["umap"]]@cell.embeddings
umap.2d <- umap.2d[rownames(pd),]
reducedDims(cds)$UMAP = umap.2d

##clustering
cds <- cluster_cells(cds)

# optional: add Seurat cluster result
cl <- setNames(plot_objects[["LGD"]]@meta.data$cluster_id, plot_objects[["LGD"]]@meta.data$sample_name )
cl <- cl[rownames(pd)]
cds@clusters$UMAP$clusters = cl

cds <- learn_graph(cds)
cds <- order_cells(cds,
					reduction_method = "UMAP",
                    root_pr_nodes=NULL,
                    root_cells=NULL,
                    verbose = TRUE)

plot_cells(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "roi")

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#Monocle works by fitting a regression model to each gene. You can specify this model to account for various factors in your experiment (time, treatment, and so on).
gene_fits <- fit_models(cds, model_formula_str = "~roi"#,
 						#expression_family = "negbinominal"
 						)
#see which of these genes have time-dependent expression
fit_coefs <- coefficient_table(gene_fits)

# just extract the time terms
roi_terms <- fit_coefs %>% filter(term != "(Intercept)")

g_table <- roi_terms %>% filter(q_value < 0.05) %>%
         select(gene_short_name, term, q_value, estimate)


x=evaluate_fits(gene_fits)


pr_graph_test_res <- graph_test(cds, neighbor_graph="knn")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

#Finding modules of co-regulated genes
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)

CS_genes <- c("Mgat4c", "Etl4", "Trhde", "Necab1","Gpc5")
CS_roi_cds <- cds[rowData(cds)$gene_short_name %in% CS_genes, ]

plot_genes_in_pseudotime(CS_roi_cds,
                         color_cells_by="roi",
                         min_expr=0.5)






########################
## Setup for Monocle
########################

setwd("C:/Users/cindy.vanvelthoven/Dropbox/AIBS/Transcriptomics/Manuscripts/LGd/code_and_data/")

library(dplyr)
library(Seurat)
library(VGAM)
library(monocle)

load("collected_data_20191220.rda")

st=format(Sys.time(), "%Y%m%d_%H%M_")

load("mouse_allcl.rda")
load("mouse_plotobjects.rda")


data <- datlist$mouse
data <- as(as.matrix(plot_objects[["LGD"]]@assays[["RNA"]]@data), 'sparseMatrix')


pd <- plot_objects[["LGD"]]@meta.data
colnames(pd)[4] <- "exp_comp_name"
pd <- new("AnnotatedDataFrame", data = pd)

gene_metadata <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = gene_metadata)

cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      #lowerDetectionLimit = 0.5,
                      expressionFamily = uninormal())

pData(cds)
fData(cds)


#Run ordering algorithm
var_genes <- plot_objects[["LGD"]][["RNA"]]@var.features
ordering_genes <- var_genes

cds <- setOrderingFilter(cds, ordering_genes)
print(dim(exprs(cds)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
cds <- reduceDimension(cds,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=4,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=0)

# First decide what you want to color your cells by
print(head(pData(cds)))

## order cells change colors and theta to match your plot
cds <- orderCells(cds)
plot_cell_trajectory(cds, 
                     color_by = "roi",
                     theta = -15,
                     show_branch_points = FALSE,
                     show_tree = TRUE,
                     cell_size = 4) + scale_color_manual(breaks = c("X", "Y", "Z"), 
                         values=c("red", "green", "blue")) + theme(legend.position = "right")

# now onwards to pseudotemporal plots
genex <- c("geneX", "geneY")

sig_gene_names <- (genex)
head(sig_gene_names)
pseudotemporalplot <- plot_pseudotime_heatmap(cds[sig_gene_names],
                        num_clusters = 9, 
                        cores = 4,
                        hmcols = NULL,
                        show_rownames = T)

pseudotemporalplot