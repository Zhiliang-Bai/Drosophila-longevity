library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(viridis)
library(Connectome)
library(cowplot)
library(tidyverse)
library(devtools)
library(ggpubr)

###1. Read-in data object of Brain with cluster information
Brain <- readRDS(file = "Brain_All_Clustered.rds")
##. Rename identity classes
#. Define a vector of new cluster names
new_names <- paste0("B", 0:21)
#. Create a named vector where old cluster names ("0" to "21") map to new names ("B0" to "B21")
old_names <- as.character(0:21)
rename_mapping <- setNames(new_names, old_names)
#. Rename the clusters
Brain <- RenameIdents(object = Brain, rename_mapping)
#. Add renamed clusters to metadata
Brain@meta.data$Cluster <- Idents(Brain)

###2. Read-in data object of Gut with cluster information
Gut <- readRDS(file = "Gut_All_Clustered.rds")
##. Rename identity classes
new_names <- paste0("G", 0:15)
old_names <- as.character(0:15)
rename_mapping <- setNames(new_names, old_names)
Gut <- RenameIdents(object = Gut, rename_mapping)
Gut@meta.data$Cluster <- Idents(Gut)

Gut <- subset(Gut,idents = c('G3','G13'), invert=TRUE) # Both 880 or Ctrl have few cells in the two clusters.
Gut@meta.data$Cluster <- Idents(Gut)

###3. Prepare integrated data object for L-R analysis
Fly <- merge(Brain, y = list(Gut))
Fly@meta.data[["RNA_snn_res.0.1"]] <- NULL
Fly@meta.data[["RNA_snn_res.0.2"]] <- NULL
Fly@meta.data[["seurat_clusters"]] <- NULL

##. Read in L-R pair file
Customlist <- read.csv(file = "Ligand_receptor_pair_high_confident.txt", sep = "\t")

##. Select only L-R genes
gene_list <- unique(c(Customlist$Gene_secreted, Customlist$Gene_receptor))
all.genes <- rownames(Fly)
common_genes <- intersect(gene_list, all.genes)

##. Normalization and Scale data
Idents(Fly) <- "Cluster"
Fly <- NormalizeData(Fly)
Fly <- ScaleData(Fly,features = common_genes)

##. Extract different conditions
Idents(Fly) <- "Condition"
levels(Fly)

Fly_880 <- subset(Fly,idents = c('B880','G880'))
Idents(Fly_880) <- "Cluster"

Fly_868 <- subset(Fly,idents = c('B868','G868'))
Idents(Fly_868) <- "Cluster"

Fly_Ctrl <- subset(Fly,idents = c('BCtrl','GCtrl'))
Idents(Fly_Ctrl) <- "Cluster"

###4. Calculate L-R connectome for each condition
Fly_880.con <- CreateConnectome(Fly_880,LR.database='custom',assay = "RNA",custom.list= Customlist,p.values = T,calculate.DOR = F)
Fly_868.con <- CreateConnectome(Fly_868,LR.database='custom',assay = "RNA",custom.list= Customlist,p.values = T,calculate.DOR = F)
Fly_Ctrl.con <- CreateConnectome(Fly_Ctrl,LR.database='custom',assay = "RNA",custom.list= Customlist,p.values = T,calculate.DOR = F)

##. Filter only significant pairs for 880
Fly_880.con_filtered <- FilterConnectome(Fly_880.con,min.pct = 0.1,min.z = 0.1,max.p = 0.05, remove.na = T)

Source.cluster <- unique(Fly_880.con_filtered$source)
#. Get elements that start with "G"
Source.cluster <- Source.cluster[grep("^G", Source.cluster)]

Target.cluster <- unique(Fly_880.con_filtered$target)
#. Get elements that start with "B"
Target.cluster <- Target.cluster[grep("^B", Target.cluster)] 

CircosPlot(Fly_880.con_filtered, sources.include = Source.cluster, targets.include = Target.cluster, min.z = 0.3,lab.cex = 0.6,cols.use=c('G0' = '#8497B0', 'G1' = '#878787', 'G2' = '#F0CE58', 'G4' = '#0FFFFF','G5' = '#DCEAF7', 
                                                                                                                                        'G6' = '#C3FF00','G7' = '#00475F','G8' = '#F297A7','G9' = '#ABDDDE','G10' = '#D7EF9B',
                                                                                                                                        'G11' = '#B487B7','G12' = '#1D65A6','G14' = '#FC6FCF','G15' = '#FFEE00',
                                                                                                                                        'B0' ='#72A2C0', 'B1' ='#ff7f0e', 'B2' ='#ffbb78','B3' = '#2ca02c', 'B4' ='#1f77b4',
                                                                                                                                        'B5' ='#9edae5', 'B6' ='#d62728', 'B7' ='#17becf', 'B8' ='#dbdb8d', 'B9' ='#ff9896', 
                                                                                                                                        'B10' ='#8c564b', 'B11' ='#7f7f7f', 'B12' = '#98df8a', 'B13' ='#c49c94', 
                                                                                                                                        'B14' = '#c7c7c7', 'B15' ='#c5b0d5','B16' = '#f7b6d2', 'B17' ='#9467bd', 
                                                                                                                                        'B18' ='#bcbd22','B19' ='#aec7e8', 'B20' ='#F3E96B','B21' ='#e377c2'))



##. Filter only significant pairs for Ctrl
Fly_Ctrl.con_filtered <- FilterConnectome(Fly_Ctrl.con,min.pct = 0.1,min.z = 0.1,max.p = 0.05, remove.na = T)

Source.cluster2 <- unique(Fly_Ctrl.con_filtered$source)
#. Get elements that start with "G"
Source.cluster2 <- Source.cluster2[grep("^G", Source.cluster2)]

Target.cluster2 <- unique(Fly_Ctrl.con_filtered$target)
#. Get elements that start with "B"
Target.cluster2 <- Target.cluster2[grep("^B", Target.cluster2)] 

CircosPlot(Fly_Ctrl.con_filtered, sources.include = Source.cluster2, targets.include = Target.cluster2, min.z = 0.5,lab.cex = 0.6,cols.use=c('G0' = '#8497B0', 'G1' = '#878787', 'G2' = '#F0CE58', 'G4' = '#0FFFFF','G5' = '#DCEAF7', 
                                                                                                                                           'G6' = '#C3FF00','G7' = '#00475F','G8' = '#F297A7','G9' = '#ABDDDE','G10' = '#D7EF9B',
                                                                                                                                           'G11' = '#B487B7','G12' = '#1D65A6','G14' = '#FC6FCF','G15' = '#FFEE00',
                                                                                                                                           'B0' ='#72A2C0', 'B1' ='#ff7f0e', 'B2' ='#ffbb78','B3' = '#2ca02c', 'B4' ='#1f77b4',
                                                                                                                                           'B5' ='#9edae5', 'B6' ='#d62728', 'B7' ='#17becf', 'B8' ='#dbdb8d', 'B9' ='#ff9896', 
                                                                                                                                           'B10' ='#8c564b', 'B11' ='#7f7f7f', 'B12' = '#98df8a', 'B13' ='#c49c94', 
                                                                                                                                           'B14' = '#c7c7c7', 'B15' ='#c5b0d5','B16' = '#f7b6d2', 'B17' ='#9467bd', 
                                                                                                                                           'B18' ='#bcbd22','B19' ='#aec7e8', 'B20' ='#F3E96B','B21' ='#e377c2'))

###5. Calculate differential L-R pairs and plot correlation between 880 and Ctrl

##. Check if two lists contain the same contents
setequal(Fly_880.con$edge, Fly_Ctrl.con$edge)

##. Create data.frames with edge and weight_sc columns
df1 <- data.frame(edge = Fly_880.con$edge, weight_880 = Fly_880.con$weight_sc)
df2 <- data.frame(edge = Fly_Ctrl.con$edge, weight_Ctrl = Fly_Ctrl.con$weight_sc)

##. Merge by "edge" to align values
merged_df <- merge(df1, df2, by = "edge", sort = FALSE)

##. Compute the difference (Fly_880 - Fly_Ctrl)
merged_df$difference <- merged_df$weight_880 - merged_df$weight_Ctrl

ggscatter(merged_df, x = "weight_Ctrl", y = "weight_880",
          color = "difference", size=0.1,
          add = "reg.line", add.params = list(color = "#424242",size= 0.5), conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson") +
  scale_color_gradient2(low = "#215F9A", mid = "white", high = "#C00000", midpoint = 0)

###6. Extract top-upregulated and downregulated L-R pairs in 880 and Plot
##. Select the top upregulated edges based on difference
top_edges <- merged_df %>%
  arrange(desc(difference)) %>%  # Sort by difference in descending order
  slice_head(n = 10000) %>%              # Select the top upregulated rows
  pull(edge)                           # Extract edge names as a vector

##. Subset Fly_880.con based on the top edges
Fly_880_top <- Fly_880.con %>% filter(edge %in% top_edges)

Source.cluster_top <- unique(Fly_880_top$source)
Source.cluster_top <- Source.cluster_top[grep("^G", Source.cluster_top)]

Target.cluster_top <- unique(Fly_880_top$target)
Target.cluster_top <- Target.cluster_top[grep("^B", Target.cluster_top)] 

CircosPlot(Fly_880_top, sources.include = Source.cluster_top, targets.include = Target.cluster_top, max.p = 0.05,lab.cex = 1,cols.use=c('G0' = '#8497B0', 'G1' = '#878787', 'G2' = '#F0CE58', 'G4' = '#0FFFFF','G5' = '#DCEAF7', 
                                                                                                                                             'G6' = '#C3FF00','G7' = '#00475F','G8' = '#F297A7','G9' = '#ABDDDE','G10' = '#D7EF9B',
                                                                                                                                             'G11' = '#B487B7','G12' = '#1D65A6','G14' = '#FC6FCF','G15' = '#FFEE00',
                                                                                                                                             'B0' ='#72A2C0', 'B1' ='#ff7f0e', 'B2' ='#ffbb78','B3' = '#2ca02c', 'B4' ='#1f77b4',
                                                                                                                                             'B5' ='#9edae5', 'B6' ='#d62728', 'B7' ='#17becf', 'B8' ='#dbdb8d', 'B9' ='#ff9896', 
                                                                                                                                             'B10' ='#8c564b', 'B11' ='#7f7f7f', 'B12' = '#98df8a', 'B13' ='#c49c94', 
                                                                                                                                             'B14' = '#c7c7c7', 'B15' ='#c5b0d5','B16' = '#f7b6d2', 'B17' ='#9467bd', 
                                                                                                                                             'B18' ='#bcbd22','B19' ='#aec7e8', 'B20' ='#F3E96B','B21' ='#e377c2'))

##. Select the top most downregulated edges (smallest difference)
down_edges <- merged_df %>%
  filter(difference < 0) %>%  # Keep only negative differences
  arrange(difference) %>%     # Sort in ascending order (smallest first)
  slice_head(n = 5000) %>%      # Select the top most negative values
  pull(edge)                  # Extract edge names as a vector

##. Subset Fly_880.con based on the most downregulated edges
Fly_880_down <- Fly_880.con %>% filter(edge %in% down_edges)

Source.cluster_down <- unique(Fly_880_down$source)
Source.cluster_down <- Source.cluster_down[grep("^G", Source.cluster_down)]

Target.cluster_down <- unique(Fly_880_down$target)
Target.cluster_down <- Target.cluster_down[grep("^B", Target.cluster_down)] 

CircosPlot(Fly_880_down, sources.include = Source.cluster_down, targets.include = Target.cluster_down, max.p = 0.05,lab.cex = 1,cols.use=c('G0' = '#8497B0', 'G1' = '#878787', 'G2' = '#F0CE58', 'G4' = '#0FFFFF','G5' = '#DCEAF7', 
                                                                                                                                       'G6' = '#C3FF00','G7' = '#00475F','G8' = '#F297A7','G9' = '#ABDDDE','G10' = '#D7EF9B',
                                                                                                                                       'G11' = '#B487B7','G12' = '#1D65A6','G14' = '#FC6FCF','G15' = '#FFEE00',
                                                                                                                                       'B0' ='#72A2C0', 'B1' ='#ff7f0e', 'B2' ='#ffbb78','B3' = '#2ca02c', 'B4' ='#1f77b4',
                                                                                                                                       'B5' ='#9edae5', 'B6' ='#d62728', 'B7' ='#17becf', 'B8' ='#dbdb8d', 'B9' ='#ff9896', 
                                                                                                                                       'B10' ='#8c564b', 'B11' ='#7f7f7f', 'B12' = '#98df8a', 'B13' ='#c49c94', 
                                                                                                                                       'B14' = '#c7c7c7', 'B15' ='#c5b0d5','B16' = '#f7b6d2', 'B17' ='#9467bd', 
                                                                                                                                       'B18' ='#bcbd22','B19' ='#aec7e8', 'B20' ='#F3E96B','B21' ='#e377c2'))

###7. Get statistical L-R pair numbers in clusters
##. Filter edges with difference > 0
upregulated_edges <- merged_df %>%
  filter(difference > 0) %>%  
  pull(edge)  # Extract edge names as a vector

# Subset Fly_880.con to keep only edges with difference > 0
Fly_880.con_Up <- Fly_880.con %>% filter(edge %in% upregulated_edges)

# Keep only items with adjusted p-values < 0.05
Fly_880.con_Up <- Fly_880.con_Up %>%
  filter(p_val_adj.lig < 0.05 & p_val_adj.rec < 0.05)

# Keep only items where percent.source and percent.target > 0.1
Fly_880.con_Up <- Fly_880.con_Up %>%
  filter(percent.source > 0.1 & percent.target > 0.1)

# Count number of edges for each unique source
source_counts <- Fly_880.con_Up %>%
  count(source, name = "num_edges_source")
source_counts_Up <- source_counts %>%
  filter(grepl("^G", source))
# Count number of edges for each unique target
target_counts <- Fly_880.con_Up %>%
  count(target, name = "num_edges_target")
target_counts_Up <- target_counts %>%
  filter(grepl("^B", target))

write.csv(source_counts_Up, file = "Source cluster_counts_Up.csv")
write.csv(target_counts_Up, file = "Target cluster_counts_Up.csv")

##. Filter edges with difference < 0
downregulated_edges <- merged_df %>%
  filter(difference < 0) %>%  
  pull(edge)  # Extract edge names as a vector

Fly_880.con_Down <- Fly_880.con %>% filter(edge %in% downregulated_edges)

Fly_880.con_Down <- Fly_880.con_Down %>%
  filter(p_val_adj.lig < 0.05 & p_val_adj.rec < 0.05)

Fly_880.con_Down <- Fly_880.con_Down %>%
  filter(percent.source > 0.1 & percent.target > 0.1)

source_counts2 <- Fly_880.con_Down %>%
  count(source, name = "num_edges_source")
source_counts_Down <- source_counts2 %>%
  filter(grepl("^G", source))

target_counts2 <- Fly_880.con_Down %>%
  count(target, name = "num_edges_target")
target_counts_Down <- target_counts2 %>%
  filter(grepl("^B", target))

write.csv(source_counts_Down, file = "Source cluster_counts_Down.csv")
write.csv(target_counts_Down, file = "Target cluster_counts_Down.csv")

###8. Calculate differential L-R pairs and Plot correlation between 880 and 868
##. Check if two lists contain the same contents
setequal(Fly_880.con$edge, Fly_868.con$edge)

##. Create data.frames with edge and weight_sc columns
df1 <- data.frame(edge = Fly_880.con$edge, weight_880 = Fly_880.con$weight_sc)
df2 <- data.frame(edge = Fly_868.con$edge, weight_868 = Fly_868.con$weight_sc)

##. Merge by "edge" to align values
merged_df <- merge(df1, df2, by = "edge", sort = FALSE)

##. Compute the difference (Fly_880 - Fly_868)
merged_df$difference <- merged_df$weight_880 - merged_df$weight_868

ggscatter(merged_df, x = "weight_868", y = "weight_880",
          color = "difference", size=0.1,
          add = "reg.line", add.params = list(color = "#424242",size= 0.5), conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson") +
  scale_color_gradient2(low = "#215F9A", mid = "white", high = "#C00000", midpoint = 0)