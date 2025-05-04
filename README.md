# Long-lived yeast-derived lipids promote longevity in Drosophila

<img width="1857" alt="Schematic" src="https://github.com/user-attachments/assets/debb5ca8-ab86-4c2e-9a35-36b8ca2f2f75" />

To investigate how lipid products from engineered yeasts influence gut physiology and lifespan, we performed single-nucleus RNA sequencing (snRNA-seq) on gut tissues from 45-day-old flies fed with lipid extracts from different yeast strains, or a control diet without lipid supplementation, strictly following the [Fly Cell Atlas](https://flycellatlas.org/) protocol. Paired brain tissues were also collected, given the critical role of the gut-brain axis in Drosophila aging.

## 1. snRNA-seq data processing and analysis

Raw sequencing reads were aligned to the Drosophila melanogaster reference genome (FlyBase r6.31) using a pre-mRNA GTF annotation established by the Fly Cell Atlas. Alignment, barcode processing, and UMI counting were performed with Cell Ranger V6.1.2 (10x Genomics) to generate a digital gene expression matrix. Quality metrics—including mean reads per nucleus, fraction of valid barcodes, and sequencing saturation—were extracted from the Cell Ranger web summary reports. Subsequent data analysis was conducted using the [Seurat V4 pipeline](https://satijalab.org/seurat/articles/get_started.html). Low-quality nuclei and potential doublets were filtered based on gene count distributions. Specifically, nuclei with 200–1,600 detected genes were retained for brain tissues, whereas the upper threshold was extended to 3,000 genes for gut tissues due to their inherently higher gene content.

## 2. Unsupervised clustering of gut and brain dataset

Standard Seurat workflow functions were applied in the following order: NormalizeData, FindVariableFeatures, ScaleData, RunPCA, and ElbowPlot to determine the dimensionality of the dataset, followed by FindNeighbors, FindClusters (at a resolution of 0.2 for the gut dataset and 0.1 for the brain dataset), and RunUMAP for visualization. Differentially expressed genes were identified using the FindMarkers function through pairwise comparisons of clusters or yeast strain groups. Module scores for predefined gene sets were calculated using the AddModuleScore function to quantify expression signatures across nuclei.

<img width="1386" alt="Screenshot 2025-05-03 at 8 03 07 PM" src="https://github.com/user-attachments/assets/f921868f-8422-4756-9eb4-0793789668f4" />

Scripts are included in the "Unsupervised clustering" folder.

## 3. Gut-Brain ligand–receptor (L–R) interaction analysis


Scripts are included in the "Gut-Brain L–R interaction analysis"





