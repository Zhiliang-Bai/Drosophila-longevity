library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(viridis)
library(nichenetr)

###1. Read-in the integrated dataset downloaded from our GEO deposition
Brain_All <- readRDS(file = "Brain_All.rds")

###2. Clustering following Seurat standard workflow
Brain_All <- NormalizeData(Brain_All, normalization.method = "LogNormalize", scale.factor = 10000)
Brain_All <- FindVariablefeatures(Brain_All, selection.method = "vst", nFeatures_C3 = 2000)
all.genes <- rownames(Brain_All)
Brain_All <- ScaleData(Brain_All, Features_C3 = all.genes)
Brain_All <- RunPCA(Brain_All, Features_C3 = VariableFeatures_C3(object = Brain_All))
ElbowPlot(Brain_All, ndims = 50)
Brain_All <- FindNeighbors(Brain_All, dims = 1:30)
Brain_All <- FindClusters(Brain_All, resolution = 0.1)
Brain_All <- RunUMAP(Brain_All, dims = 1:30)
saveRDS(Brain_All, file = "Brain_All_Clustered.rds")

DimPlot(Brain_All, reduction = "umap",pt.size=0.1,cols = c('0' ='#72A2C0', '1' ='#ff7f0e', '2' ='#ffbb78','3' = '#2ca02c', '4' ='#1f77b4','5' ='#9edae5', '6' ='#d62728', '7' ='#17becf', '8' ='#dbdb8d', '9' ='#ff9896', '10' ='#8c564b', '11' ='#7f7f7f', '12' = '#98df8a', '13' ='#c49c94', 
                                                           '14' = '#c7c7c7', '15' ='#c5b0d5','16' = '#f7b6d2', '17' ='#9467bd', '18' ='#bcbd22','19' ='#aec7e8', '20' ='#F3E96B','21' ='#e377c2'))

###3. Find markers for every cluster compared to all remaining cells, report only the positive ones
Brain.markers <- FindAllMarkers(Brain_All, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
Gene_plot <- unique(c("pros","dati","Dop2R","VGlut","DIP-eta","CG46448","CG14989","Nos","Gad1","axo","CG42613","wtrw",
                      "Rgk1","Imp","sNPF","CG4797","dve","pdm3","SoxN","Drgx","acj6","shakB","grn","AstA-R1","sNPF","ey","Mef2",
                      "side-IV","Tk","kek3","Eaat1","CG2016","orb","nAChRalpha6","acj6","beat-Ic","Octbeta1R","bi","CG44837",
                      "Eaat1","Gat","alrm","Cyp4g15","e","CG43795","Dop1R2","jdp","mamo","sosie","Rh7","Mmp2","Rh2","Arr2","trpl",
                      "beat-IIb","Tk","beat-IIa","lncRNA:CR45223","acj6","kn","Vmat","DAT","DIP-delta","Dh31","5-HT7","Poxn"))
Marker.averages <- AverageExpression(Brain_All,Features_C3=Gene_plot,return.seurat = TRUE)
DoHeatmap(Marker.averages, Features_C3 = Gene_plot, size = 3,
          draw.lines = FALSE)+ scale_fill_gradientn(colors=c("#192E5B","#1D65A6","#72A2C0","#F3E96B","#F2A104"))

###4. Get feature Score of a group of genes
B3list <- list(c("axo", "CG42613", "wtrw"))
Brain_All <- AddModuleScore(object = Brain_All, Features_C3 = B3list, name = "B3list")
FeaturePlot(object = Brain_All, Features_C3 = "B3list1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="viridis")

B4list <- list(c("Rgk1", "Imp", "sNPF"))
Brain_All <- AddModuleScore(object = Brain_All, Features_C3 = B4list, name = "B4list")
FeaturePlot(object = Brain_All, Features_C3 = "B4list1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="mako")

B5list <- list(c("CG4797", "dve", "pdm3"))
Brain_All <- AddModuleScore(object = Brain_All, Features_C3 = B5list, name = "B5list")
FeaturePlot(object = Brain_All, Features_C3 = "B5list1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="turbo")

B8list <- list(c("crb", "ey", "Mef2"))
Brain_All <- AddModuleScore(object = Brain_All, Features_C3 = B8list, name = "B8list")
FeaturePlot(object = Brain_All, Features_C3 = "B8list1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="cividis")

###5. Find DEGs between 880 and Ctrl in a specific cluster
Brain_Cluster <- subset(Brain_All,idents = '3')
Idents(Brain_Cluster)<-"Condition"
Cluster.markers <- FindMarkers(Brain_Cluster, ident.1 ='B880', only.pos = FALSE,ident.2 ='BCtrl' , min.pct = 0.1,logfc.threshold=0.1)
write.csv(Cluster.markers, file = "DEGs_B880vs.BCtrl_Cluster3.csv", row.names = TRUE, quote = TRUE)

###6. Get expression of selected genes in a specific cluster
##. Get expression of selected genes in cluster 3
Brain_Cluster3 <- subset(Brain_All,idents = '3')
Features_C3 <- c("GABA-B-R1","Dop2R","Rdl","Desat1","Baldspot","Acsl")
Idents(Brain_Cluster3)<-"Condition"

Cluster3_880 <-subset(x = Brain_Cluster3, idents ="B880")
Select_Gene_C3_880 <- FetchData(object = Cluster3_880, vars = Features_C3)
write.csv(Select_Gene_C3_880,file="Cluster3_880_Select_Gene_C3 genes in each cell.csv")

Cluster3_Ctrl <-subset(x = Brain_Cluster3, idents ="BCtrl")
Select_Gene_C3_Ctrl <- FetchData(object = Cluster3_Ctrl, vars = Features_C3)
write.csv(Select_Gene_C3_Ctrl,file="Cluster3_Ctrl_Select_Gene_C3 genes in each cell.csv")

Cluster3_868 <-subset(x = Brain_Cluster3, idents ="B868")
Select_Gene_C3_868 <- FetchData(object = Cluster3_868, vars = Features_C3)
write.csv(Select_Gene_C3_868,file="Cluster3_868_Select_Gene_C3 genes in each cell.csv")

##. Get expression of selected genes in cluster 4
Brain_Cluster4 <- subset(Brain_All,idents = '4')
Features_C4 <- c("phu","cin","Sptr","Hex-A","Eno","Tpi")
Idents(Brain_Cluster4)<-"Condition"

Cluster4_880 <-subset(x = Brain_Cluster4, idents ="B880")
Select_Gene_C4_880 <- FetchData(object = Cluster4_880, vars = Features_C4)
write.csv(Select_Gene_C4_880,file="Cluster4_880_Select_Gene_C4 genes in each cell.csv")

Cluster4_Ctrl <-subset(x = Brain_Cluster4, idents ="BCtrl")
Select_Gene_C4_Ctrl <- FetchData(object = Cluster4_Ctrl, vars = Features_C4)
write.csv(Select_Gene_C4_Ctrl,file="Cluster4_Ctrl_Select_Gene_C4 genes in each cell.csv")

Cluster4_868 <-subset(x = Brain_Cluster4, idents ="B868")
Select_Gene_C4_868 <- FetchData(object = Cluster4_868, vars = Features_C4)
write.csv(Select_Gene_C4_868,file="Cluster4_868_Select_Gene_C4 genes in each cell.csv")

##. Get expression of selected genes in cluster 5
Brain_Cluster5 <- subset(Brain_All,idents = '5')
Features_C5 <- c("GABA-B-R1","Dop2R","Rdl","Hex-A","Gapdh1","Pgk")
Idents(Brain_Cluster5)<-"Condition"

Cluster5_880 <-subset(x = Brain_Cluster5, idents ="B880")
Select_Gene_C5_880 <- FetchData(object = Cluster5_880, vars = Features_C5)
write.csv(Select_Gene_C5_880,file="Cluster5_880_Select_Gene_C5 genes in each cell.csv")

Cluster5_Ctrl <-subset(x = Brain_Cluster5, idents ="BCtrl")
Select_Gene_C5_Ctrl <- FetchData(object = Cluster5_Ctrl, vars = Features_C5)
write.csv(Select_Gene_C5_Ctrl,file="Cluster5_Ctrl_Select_Gene_C5 genes in each cell.csv")

Cluster5_868 <-subset(x = Brain_Cluster5, idents ="B868")
Select_Gene_C5_868 <- FetchData(object = Cluster5_868, vars = Features_C5)
write.csv(Select_Gene_C5_868,file="Cluster5_868_Select_Gene_C5 genes in each cell.csv")

##. Get expression of selected genes in cluster 8
Brain_Cluster8 <- subset(Brain_All,idents = '8')
Features_C8 <- c("ND-ACP","COX7A","ATPsynG","RpL8","RpL35A","RpS26")
Idents(Brain_Cluster8)<-"Condition"

Cluster8_880 <-subset(x = Brain_Cluster8, idents ="B880")
Select_Gene_C8_880 <- FetchData(object = Cluster8_880, vars = Features_C8)
write.csv(Select_Gene_C8_880,file="Cluster8_880_Select_Gene_C8 genes in each cell.csv")

Cluster8_Ctrl <-subset(x = Brain_Cluster8, idents ="BCtrl")
Select_Gene_C8_Ctrl <- FetchData(object = Cluster8_Ctrl, vars = Features_C8)
write.csv(Select_Gene_C8_Ctrl,file="Cluster8_Ctrl_Select_Gene_C8 genes in each cell.csv")

Cluster8_868 <-subset(x = Brain_Cluster8, idents ="B868")
Select_Gene_C8_868 <- FetchData(object = Cluster8_868, vars = Features_C8)
write.csv(Select_Gene_C8_868,file="Cluster8_868_Select_Gene_C8 genes in each cell.csv")

###7. Get expression of EGFR receptor genes across clusters
Idents(Brain_All)<-"RNA_snn_res.0.1"
levels(x=Brain_All) <- c("21","20","19","18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1","0")
DotPlot(Brain_All, features = c("Egfr","sev","btl","htl"), cols = c("#215F99","#C01900"), dot.scale = 6,scale=FALSE) + RotatedAxis()

Idents(Brain_All)<-"Condition"
Brain_880 <- subset(Brain_All,idents = 'B880')
EGFR_880 <- FetchData(object = Brain_880, vars = c("Egfr","sev","btl","htl"))
write.csv(EGFR_880,file="EGFR_880 genes in each cell.csv")

Brain_Ctrl <- subset(Brain_All,idents = 'BCtrl')
EGFR_Ctrl <- FetchData(object = Brain_Ctrl, vars = c("Egfr","sev","btl","htl"))
write.csv(EGFR_Ctrl,file="EGFR_Ctrl genes in each cell.csv")

Brain_868 <-subset(x = Brain_All, idents ="B868")
EGFR_868 <- FetchData(object = Brain_868, vars = c("Egfr","sev","btl","htl"))
write.csv(EGFR_868,file="EGFR_868 genes in each cell.csv")
