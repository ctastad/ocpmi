# Dependencies, Initialize

library(dplyr)
library(Seurat)

setwd("~/starr_lab/ocpmi/code/")
pt59.data <- Read10X(data.dir = "./pt59/")
pt59 <- CreateSeuratObject(counts = pt59.data, project = "pt59", min.cells = 3, min.features = 200)
pt59


# Survey data

pt59.data[c("CD3D","TCL1A","MS4A1"), 1:30]


# Pre-processing workflow

pt59[["percent.mt"]] <- PercentageFeatureSet(object = pt59, pattern = "^MT-")


# Visualize QC metrics

head(x = pt59@meta.data, 5)
VlnPlot(object = pt59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = pt59, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = pt59, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))

pt59 <- subset(x = pt59, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing the data

pt59 <- NormalizeData(object = pt59, normalization.method = "LogNormalize", scale.factor = 1e4)

# Identification of highly variable features

pt59 <- FindVariableFeatures(object = pt59, selection.method = 'vst', nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = pt59), 10)
plot1 <- VariableFeaturePlot(object = pt59)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
CombinePlots(plots = list(plot1, plot2))

# Scaling the data

all.genes <- rownames(x = pt59)
pt59 <- ScaleData(object = pt59, features = all.genes)

# Linear dimensional reduction

pt59 <- RunPCA(object = pt59, features = VariableFeatures(object = pt59))
print(x = pt59[['pca']], dims = 1:5, nfeatures = 5)
