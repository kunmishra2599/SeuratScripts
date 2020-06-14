library(Seurat)
library(dplyr)

# Creates a subset of cells based on a certain cluster, can be modified for other purposes
obj.subcluster <- subset(obj, idents = X)

# Runs UMAP on subset
obj.subcluster <- FindNeighbors(obj.subcluster, ...)
obj.subcluster <- FindClusters(obj.subcluster, ...)
obj.subcluster <- RunUMAP(obj.subcluster, ...)
DimPlot(obj.subcluster, reduction = "umap", label = F, pt.size = 0.6, repel = T, label.size = 5)

# Adds back into the original cluster as a new identity "subcluster"
obj$subcluster <- as.character(Idents(obj))
obj$subcluster[Cells(obj.subcluster)] <- paste("subcluster", Idents(obj.subcluster))
DimPlot(obj, group.by = "subcluster", label = F)
