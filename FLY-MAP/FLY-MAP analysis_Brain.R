library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(viridis)
library(nichenetr)
source("utils.R")
library(biomaRt)
library("Hmisc")
library(corrplot)

###1. Construct the FLY-MAP pathways
pathways <- gmtPathways("Drosophila.KEGG.gmt")
ensembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
##. Example list of CG numbers from pathways
cg_list <- unique(unlist(pathways))  # Extract unique CG IDs
##. Retrieve gene names using 'flybase_annot_id' as filter
gene_mapping <- getBM(
  attributes = c("flybase_annot_id", "external_gene_name"),
  filters = "flybase_annot_id",
  values = cg_list,
  mart = ensembl
)
##. Convert to a named vector for easy lookup
cg_to_gene <- setNames(gene_mapping$external_gene_name, gene_mapping$flybase_annot_id)
##. Function to convert CG IDs to gene names
convert_pathway <- function(cg_genes) {
  gene_names <- cg_to_gene[cg_genes]  # Lookup in named vector
  gene_names[is.na(gene_names)] <- cg_genes[is.na(gene_names)]  # Keep CG ID if no match
  return(gene_names)
}
##. Apply to all pathways
pathways_named <- lapply(pathways, convert_pathway)
##. Keep only metabolism related pathways
pathways_final <- pathways_named[-(1:49)]  # Remove first 49 items that are not related to metabolism

###2. FLY-MAP clustering analysis on Brain
##. Read-in the brain dataset
Brain_All <- readRDS(file = "Brain_All_Filtered.rds")
Brain <- subset(x = Brain_All, idents = c("Brain_880", "Brain_Ctrl"))

metabolics <- unique(as.vector(unname(unlist(pathways_final))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(Brain)),row.names = rownames(Brain))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE
Brain_Met_metabolics <- subset(row_data,subset=metabolic=="TRUE")
Brain_Met_metabolics <- rownames(Brain_Met_metabolics)

##. Seurat standard clustering using metabolic gene set
Brain_Met <- NormalizeData(Brain, normalization.method = "LogNormalize", scale.factor = 10000)
Brain_Met <- ScaleData(Brain_Met, features = Brain_Met_metabolics)
Brain_Met <- RunPCA(Brain_Met, features = Brain_Met_metabolics)
ElbowPlot(Brain_Met, ndims = 50)
Brain_Met <- FindNeighbors(Brain_Met, dims = 1:13)
Brain_Met <- FindClusters(Brain_Met, resolution = 0.2)
Brain_Met <- RunUMAP(Brain_Met, dims = 1:13)

DimPlot(Brain_Met, reduction = "umap",pt.size=0.1,cols = c( '0' ='#1f77b4', '1' ='#d62728', '2' ='#F2A104','3' = '#c5b0d5', '4' ='#2ca02c',
                                                            '5' = '#c7c7c7', '6' ='#c5b0d5','7' = '#f7b6d2', '8' ='#9467bd', '9' ='#bcbd22','10' ='#aec7e8'))
##. Calculate the cell proportion in each metabolic cluster
Met_cluster_Proportion <- prop.table(table(Brain_Met@meta.data$RNA_snn_res.0.2, Brain_Met@meta.data$Condition), margin = 2)
write.csv(Met_cluster_Proportion,file="Brain_Metablism cell proportion of each condition in metabolic each cluster.csv")

###3. Calculate the pathway activities
pathway_names <- names(pathways_final)
# replace "/" with " " in all the pathway names
pathway_names <- gsub("/", "", pathway_names)
pathway_names <- gsub("\"", "", pathway_names)
pathway_names <- gsub("-", "_", pathway_names) 
pathway_names <- gsub(" ", "_", pathway_names)
pathway_names <- gsub(",", "_", pathway_names)

for(p in 1:49){
  genes <- pathways_final[[p]]
  genes_comm <- list (intersect(genes, Brain_Met_metabolics))
  if(length(genes_comm) < 1) next
  
  Path_name <-pathway_names[p]
  Brain_Met <- AddModuleScore(object = Brain_Met, features = genes_comm, name = Path_name)
}
# Genes defining the 50th pathway (D Arginine and D ornithine metabolism) were not detected in the dataset, thus this pathway was not included.
for(p in 51:82){
  genes <- pathways_final[[p]]
  genes_comm <- list (intersect(genes, Brain_Met_metabolics))
  if(length(genes_comm) < 1) next
  
  Path_name <-pathway_names[p]
  Brain_Met <- AddModuleScore(object = Brain_Met, features = genes_comm, name = Path_name)
}

###4. Extract expression and run global correlation matrix analysis for 880 and Ctrl separately
Met_list_plot <- pathway_names[-50]
Met_list_plot <- paste0(Met_list_plot, "1")

Idents(Brain_Met) <- "Condition"
Brain_Met_880 <-subset(x = Brain_Met, idents ="B880")
Metexp_880 <- FetchData(object = Brain_Met_880, vars = Met_list_plot)
res <-rcorr(as.matrix(Metexp_880),type="pearson")
corrplot(res$r,  col=rev(COL2('RdBu', 100)), type = "upper", order = "hclust", 
         p.mat = res$P, tl.col = "black",tl.cex = 0.5,sig.level = 0.05, insig = "blank")

Brain_Met_Ctrl <-subset(x = Brain_Met, idents ="BCtrl")
Metexp_Ctrl <- FetchData(object = Brain_Met_Ctrl, vars = Met_list_plot)
res <-rcorr(as.matrix(Metexp_Ctrl),type="pearson")
corrplot(res$r,  col=rev(COL2('RdBu', 100)), type = "lower", order = "hclust", 
         p.mat = res$P, tl.col = "black",tl.cex = 0.5,sig.level = 0.05, insig = "blank")