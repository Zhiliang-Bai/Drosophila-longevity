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

###2. FLY-MAP clustering analysis on Gut
##. Read-in the gut dataset
Gut_All <- readRDS(file = "Gut_All_Clustered.rds")
Idents(Gut_All) <- "RNA_snn_res.0.2"
new_names <- paste0("G", 0:15)
old_names <- as.character(0:15)
rename_mapping <- setNames(new_names, old_names)
Gut_All <- RenameIdents(object = Gut_All, rename_mapping)
Gut_All@meta.data$Cluster <- Idents(Gut_All)

Idents(Gut_All) <- "Condition"
Gut <- subset(x = Gut_All, idents = c("G880", "GCtrl"))

metabolics <- unique(as.vector(unname(unlist(pathways_final))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(Gut)),row.names = rownames(Gut))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE
Gut_Met_metabolics <- subset(row_data,subset=metabolic=="TRUE")
Gut_Met_metabolics <- rownames(Gut_Met_metabolics)

##. Seurat standard clustering using metabolic gene set
Gut_Met <- NormalizeData(Gut, normalization.method = "LogNormalize", scale.factor = 10000)
Gut_Met <- ScaleData(Gut_Met, features = Gut_Met_metabolics)
Gut_Met <- RunPCA(Gut_Met, features = Gut_Met_metabolics)
ElbowPlot(Gut_Met, ndims = 50)
Gut_Met <- FindNeighbors(Gut_Met, dims = 1:30)
Gut_Met <- FindClusters(Gut_Met, resolution = 0.2)
Gut_Met <- RunUMAP(Gut_Met, dims = 1:30)

DimPlot(Gut_Met, reduction = "umap",pt.size=0.1,cols = c( '0' ='#dbdb8d', '1' ='#d62728', '2' ='#ffbb78','3' = '#2ca02c', '4' ='#1f77b4',
                                                          '5' = '#c7c7c7', '6' ='#c5b0d5','7' = '#f7b6d2', '8' ='#9467bd', '9' ='#bcbd22','10' ='#aec7e8'))

##. Calculate the cell proportion in each metabolic cluster
Met_cluster_Proportion <- prop.table(table(Gut_Met@meta.data$RNA_snn_res.0.2, Gut_Met@meta.data$Condition), margin = 2)
write.csv(Met_cluster_Proportion,file="Gut_Metablism cell proportion of each condition in each metabolic cluster.csv")

##. Extract the cell-type annotations from the prior whole-transcriptome clustering and calculate composition
Cell_type_Proportion <- prop.table(table(Gut_Met@meta.data$Cluster, Gut_Met@meta.data$RNA_snn_res.0.2), margin = 2)
write.csv(Cell_type_Proportion,file="Gut_Metablism cell Proportion of each cell type in each metabolic cluster.csv")

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
  genes_comm <- list (intersect(genes, Gut_Met_metabolics))
  if(length(genes_comm) < 1) next
  
  Path_name <-pathway_names[p]
  Gut_Met <- AddModuleScore(object = Gut_Met, features = genes_comm, name = Path_name)
}
# Genes defining the 50th pathway (D Arginine and D ornithine metabolism) were not detected in the dataset, thus this pathway was not included.
for(p in 51:82){
  genes <- pathways_final[[p]]
  genes_comm <- list (intersect(genes, Gut_Met_metabolics))
  if(length(genes_comm) < 1) next
  
  Path_name <-pathway_names[p]
  Gut_Met <- AddModuleScore(object = Gut_Met, features = genes_comm, name = Path_name)
}

FeaturePlot(object = Gut_Met, features = "Valine__leucine_and_isoleucine_biosynthesis1",order=TRUE, min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="magma")
FeaturePlot(object = Gut_Met, features = "Oxidative_phosphorylation1",order=TRUE, min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="viridis")

###4. Analyze and visualize the pathway activities
Met_list_plot <- pathway_names[-50]
Met_list_plot <- paste0(Met_list_plot, "1")

Gut_Met_Sub <- subset(Gut_Met,idents = '7', invert = TRUE) # The cell number in this cluster was extremely low, thus was not included for the analysis. 
Idents(Gut_Met_Sub) <- "RNA_snn_res.0.2"
levels (Gut_Met_Sub)
DotPlot(Gut_Met_Sub, features = Met_list_plot,dot.min=0.15,dot.scale=6)+ RotatedAxis()+scale_color_viridis(discrete = FALSE, option="viridis")

###5. Extract expression and run global correlation matrix analysis for 880 and Ctrl separately
Idents(Gut_Met) <- "Condition"
Gut_Met_Sub_880 <-subset(x = Gut_Met, idents ="G880")
Metexp_880 <- FetchData(object = Gut_Met_Sub_880, vars = Met_list_plot)
res <-rcorr(as.matrix(Metexp_880),type="pearson")
corrplot(res$r,  col=rev(COL2('RdBu', 100)), type = "upper", order = "hclust", 
         p.mat = res$P, tl.col = "black",tl.cex = 0.5,sig.level = 0.05, insig = "blank")


Gut_Met_Sub_Ctrl <-subset(x = Gut_Met_Sub, idents ="GCtrl")
Metexp_Ctrl <- FetchData(object = Gut_Met_Sub_Ctrl, vars = Met_list_plot)
res <-rcorr(as.matrix(Metexp_Ctrl),type="pearson")
corrplot(res$r,  col=rev(COL2('RdBu', 100)), type = "lower", order = "hclust", 
         p.mat = res$P, tl.col = "black",tl.cex = 0.5,sig.level = 0.05, insig = "blank")