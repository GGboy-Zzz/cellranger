###### Perform standard Seurat pipeline
# ZY
# 3 Nov 2023
# refer to: https://satijalab.org/seurat/articles/pbmc3k_tutorial

library(Seurat)

setwd("/home/zhuyong/Wangshaowei_10X_6samples/results")
ohsv1 <- Read10X("/data/share_for_user/1207_zhangjunjLab/OHSV1/outs/counts_soupX")
ohsv1 <- CreateSeuratObject(ohsv1,project = "ohsv1",min.cells = 50,min.features = 500)
ohsv1[["percent.mt"]] <- PercentageFeatureSet(ohsv1,pattern = "^mt-")

png("./seurat_results/before_QC.png")
VlnPlot(object = ohsv1,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol = 3)
dev.off()

plot1 <- FeatureScatter(ohsv1, feature1 = "nCount_RNA", feature2 = "percent.mt")+NoLegend()
plot2 <- FeatureScatter(ohsv1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# subset , 阈值设置?
ohsv1 <- subset(ohsv1,nFeature_RNA > 500 & percent.mt < 20 & nCount_RNA>1000)
ohsv1 <- SCTransform(ohsv1,variable.features.n = 2000)

top10 <- head(VariableFeatures(ohsv1), 10)
plot1 <- VariableFeaturePlot(ohsv1)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

ohsv1 <- RunPCA(ohsv1)
ElbowPlot(ohsv1,ndims = 20)
ohsv1 <- FindNeighbors(ohsv1,dims = 1:15)
ohsv1 <- FindClusters(ohsv1)
ohsv1 <- RunTSNE(ohsv1)
DimPlot(ohsv1,reduction="tsne")

