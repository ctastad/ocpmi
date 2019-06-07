# Dependencies, Initialize

library(dplyr)
library(Seurat)

pt59.data <- Read10X(data.dir = "~/pt59/")
pt59 <- CreateSeuratObject(counts = pt59.data, project = "pt59", min.cells = 3, min.features = 200)
pt59

# Survey data

pt59.data[c("CD3D","TCL1A","MS4A1"), 1:30]
dense.size <- object.size(x = as.matrix(x = pt59.data))
dense.size

# Pre-processing workflow

pt59[["percent.mt"]] <- PercentageFeatureSet(object = pt59, pattern = "^MT-")

## Visualize QC metrics

VlnPlot(object = pt59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


## Visualize QC metrics

VlnPlot(object = pt59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
