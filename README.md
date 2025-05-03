# Long-lived yeast-derived lipids promote longevity in Drosophila

<img width="1857" alt="Schematic" src="https://github.com/user-attachments/assets/debb5ca8-ab86-4c2e-9a35-36b8ca2f2f75" />

To investigate how lipid products from engineered yeasts influence gut physiology and lifespan, we performed single-nucleus RNA sequencing (snRNA-seq) on gut tissues from 45-day-old flies fed with lipid extracts from different yeast strains, or a control diet without lipid supplementation, strictly following the Fly Cell Atlas protocol. Paired brain tissues were also collected, given the critical role of the gut-brain axis in Drosophila aging.

## 1. snRNA-seq data processing and analysis

Raw sequencing reads were aligned to the Drosophila melanogaster reference genome (FlyBase r6.31) using a pre-mRNA GTF annotation established by the Fly Cell Atlas. Alignment, barcode processing, and UMI counting were performed with Cell Ranger V6.1.2 (10x Genomics) to generate a digital gene expression matrix. Quality metrics—including mean reads per nucleus, fraction of valid barcodes, and sequencing saturation—were extracted from the Cell Ranger web summary reports. Subsequent data analysis was conducted using the [Seurat V4 pipeline](https://satijalab.org/seurat/articles/get_started.html). Low-quality nuclei and potential doublets were filtered based on gene count distributions. Specifically, nuclei with 200–1,600 detected genes were retained for brain tissues, whereas the upper threshold was extended to 3,000 genes for gut tissues due to their inherently higher gene content.

## 2. Unsupervised clustering of gut and brain dataset
Scripts are included in the "Unsupervised clustering" folder.



