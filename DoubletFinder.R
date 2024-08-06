###### Perform DoubletFinder
# ZY
# 3 Nov 2023
# refer to: https://github.com/chris-mcginnis-ucsf/DoubletFinder

library(DoubletFinder)
library(Seurat)

## pk Identification(no ground-truth)-----------------------------------------------------
sweep.res.list_ohsv1 <- paramSweep_v3(ohsv1, PCs = 1:15, sct = TRUE)
sweep.stats_ohsv1 <- summarizeSweep(sweep.res.list_ohsv1, GT = FALSE)
bcmvn_ohsv1 <- find.pK(sweep.stats_ohsv1)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- ohsv1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(ohsv1@meta.data))          ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
ohsv1 <- doubletFinder_v3(ohsv1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

ohsv1 <- doubletFinder_v3(ohsv1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_2146", sct = TRUE)

## Plot results ---------------------------------------------------------------------------

ohsv1@meta.data[,"DF_hi.lo"] <- ohsv1@meta.data$DF.classifications_0.25_0.09_2146
ohsv1@meta.data$DF_hi.lo[which(ohsv1@meta.data$DF_hi.lo == "Doublet" & ohsv1@meta.data$DF.classifications_0.25_0.09_1967 == "Singlet")] <- "Doublet-Low Confidience"
ohsv1@meta.data$DF_hi.lo[which(ohsv1@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(ohsv1@meta.data$DF_hi.lo)
# Doublet-High Confidience  Doublet-Low Confidience                  Singlet 
# 198                       24                     5041 

## 结果展示，分类结果在pbmc@meta.data中
png("./output/2_cell_annotation/2_doubletFinder.png",2500,1800,res=300)
DimPlot(ohsv1, reduction = "tsne", group.by ="DF_hi.lo",cols =c("black","red","gold"))
dev.off()
