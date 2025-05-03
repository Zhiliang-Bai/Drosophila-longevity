library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(viridis)
library(nichenetr)

###1. Read-in the integrated dataset downloaded from our GEO deposition
Gut_All <- readRDS(file = "Gut_All.rds")

###2. Clustering following Seurat standard workflow
Gut_All <- NormalizeData(Gut_All, normalization.method = "LogNormalize", scale.factor = 10000)
Gut_All <- FindVariableFeatures(Gut_All, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Gut_All)
Gut_All <- ScaleData(Gut_All, features = all.genes)
Gut_All <- RunPCA(Gut_All, features = VariableFeatures(object = Gut_All))
ElbowPlot(Gut_All, ndims = 50)
Gut_All <- FindNeighbors(Gut_All, dims = 1:30)
Gut_All <- FindClusters(Gut_All, resolution = 0.2)
Gut_All <- RunUMAP(Gut_All, dims = 1:30)
saveRDS(Gut_All, file = "Gut_All_Clustered.rds")

DimPlot(Gut_All, reduction = "umap",pt.size=0.1,cols = c('0' = '#8497B0', '1' = '#878787', '2' = '#F0CE58', '3' = '#EB545C','4' = '#0FFFFF', '5' = '#DCEAF7', '6' = '#C3FF00','7' = '#00475F',
                                                         '8' = '#F297A7','9' = '#ABDDDE','10' = '#D7EF9B','11' = '#B487B7','12' = '#1D65A6','13' = '#F2A104','14' = '#FC6FCF','15' = '#FFEE00'))

###3. Find markers for every cluster compared to all remaining cells, report only the positive ones
Gut.markers <- FindAllMarkers(Gut_All, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
Gut.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> top3
Gene_plot <- unique(c("betaTry","alphaTry","epsilonTry","whd","LpR1","GstT4","Pde1c","lectin-37Db","sens-2","pros","para","Sh","ZnT35C","CG3168","ct",
                      "wat","ELOVL","mgl","N","Dl","Sox21a","lab","CAH1","Ptx1","Yp1","apolpp","IM4","CG34324","CG34220","Muc68D","AstA","Tk","CCHa2",
                      "Scp2","Dgk","hth","Mhc","MsR1","sls","Obp44a","pdm3","zfh2","Octalpha2R","Lkr","tsh","5-HT2A","5-HT7","Sox21b"))

Marker.averages <- AverageExpression(Gut_All,features=Gene_plot,return.seurat = TRUE)
DoHeatmap(Marker.averages, features = Gene_plot, size = 3,
          draw.lines = FALSE)+ scale_fill_gradientn(colors=c("#192E5B","#1D65A6","#72A2C0","#F3E96B","#F2A104"))

###4. Get feature Score of a group of genes
Enterocytes <- list(c("Pde1c", "lectin-37Db", "sens-2"))
Gut_All <- AddModuleScore(object = Gut_All, features = Enterocytes, name = "Enterocytes")
FeaturePlot(object = Gut_All, features = "Enterocytes1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="mako")

Enteroblast <- list(c("N", "Dl", "Sox21a"))
Gut_All <- AddModuleScore(object = Gut_All, features = Enteroblast, name = "Enteroblast")
FeaturePlot(object = Gut_All, features = "Enteroblast1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="magma")

###5. Find DEGs between 880 and Ctrl in specific cluster
Gut_Cluster2 <- subset(Gut_All,idents = '2')
Idents(Gut_Cluster2)<-"Condition"
Cluster.markers <- FindMarkers(Gut_Cluster2, ident.1 ='G880', only.pos = FALSE,ident.2 ='GCtrl')
write.csv(Cluster.markers, file = "DEGs_G880vs.GCtrl_Cluster15.csv", row.names = TRUE, quote = TRUE)

###6. Get expression of selected genes in cluster 2 from each condition
features <- c("Atg1","Atg8a","Atg17","Atg13")
Idents(Gut_Cluster2)<-"Condition"

Cluster2_880 <-subset(x = Gut_Cluster2, idents ="G880")
Autophage_880 <- FetchData(object = Cluster2_880, vars = features)
write.csv(Autophage_880,file="Cluster2_880_Autophage genes in each cell.csv")

Cluster2_Ctrl <-subset(x = Gut_Cluster2, idents ="GCtrl")
Autophage_Ctrl <- FetchData(object = Cluster2_Ctrl, vars = features)
write.csv(Autophage_Ctrl,file="Cluster2_Ctrl_Autophage genes in each cell.csv")

Cluster2_868 <-subset(x = Gut_Cluster2, idents ="G868")
Autophage_868 <- FetchData(object = Cluster2_868, vars = features)
write.csv(Autophage_868,file="Cluster2_868_Autophage genes in each cell.csv")

###7. Get expression of EGFR-related ligand genes across conditions
Idents(Gut_All)<-"Condition"
Gut_880 <- subset(Gut_All,idents = 'G880')
EGFR_880 <- FetchData(object = Gut_880, vars = c("bnl","pyr","vn","ths"))
write.csv(EGFR_880,file="EGFR_880 ligand genes in each cell.csv")

Gut_Ctrl <- subset(Gut_All,idents = 'GCtrl')
EGFR_Ctrl <- FetchData(object = Gut_Ctrl, vars = c("bnl","pyr","vn","ths"))
write.csv(EGFR_Ctrl,file="EGFR_Ctrl ligand genes in each cell.csv")

Gut_868 <-subset(x = Gut_All, idents ="G868")
EGFR_868 <- FetchData(object = Gut_868, vars = c("bnl","pyr","vn","ths"))
write.csv(EGFR_868,file="EGFR_868 ligand genes in each cell.csv")

###8. Get expression of Circadian rhythm-related genes across conditions
Idents(Gut_All)<-"Condition"
Gut_880 <- subset(Gut_All,idents = 'G880')
Circadian_880 <- FetchData(object = Gut_880, vars = c("Clk","Pdp1","sgg","tim","vri"))
write.csv(Circadian_880,file="Circadian_880 genes in each cell.csv")

Gut_868 <-subset(x = Gut_All, idents ="G868")
Circadian_868 <- FetchData(object = Gut_868, vars = c("Clk","Pdp1","sgg","tim","vri"))
write.csv(Circadian_868,file="Circadian_868 genes in each cell.csv")