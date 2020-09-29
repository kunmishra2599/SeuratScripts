library(ClusterProfiler)
library(Seurat)
library(org.Hs.eg.db)

cell <-
  FindMarkers(
    subset(SeuratObj, subset = celltype == "Cell Name"),
    ident.1 = "Group 1",
    ident.2 = "Group 2",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

cell_adj <- subset(cell, subset = p_val_adj < 0.1)


cell_go <-
  enrichGO(
    rownames(cell_adj),
    keyType = "SYMBOL",
    ont = "BP",
    qvalueCutoff = 0.1,
    minGSSize = 5,
    readable = F,
    OrgDb = org.Hs.eg.db,
    universe = rownames(subset(SeuratObj, subset = celltype == "Cell Name"))[which(rowSums(
      subset(SeuratObj, subset = celltype == "Cell Name")@assays$RNA@data > 0
    ) > 0.01 * (ncol(
      subset(SeuratObj, subset = celltype == "Cell Name")
    )))]
  )
