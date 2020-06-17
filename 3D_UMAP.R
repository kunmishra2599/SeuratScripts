library(plotly)
library(Seurat)

# Read in raw data files ################################################################

obj_files <- Read10X(data.dir = "...")
obj <- CreateSeuratObject(counts = obj_files, project = "Project_Name", min.cells = X, min.features = Y)

# Running dimension reduction ###########################################################

obj <- NormalizeData(object = obj, normalization.method = "LogNormalize")
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
ElbowPlot(obj)  # Pick PC cut-off


# Running UMAP on data ##################################################################

obj <- RunUMAP(obj, dims = 1:20, n.components = 3L, seed.use = 2020)
DimPlot(obj, group.by = "PlateletsSubcluster",label = T, repel = T)  # Getting 2D plot with 3D


# Extracting out embeddings in each dimension  ##########################################

UMAP_1 <- obj@reductions$umap@cell.embeddings[,1]
UMAP_2 <- obj@reductions$umap@cell.embeddings[,2]
UMAP_3 <- obj@reductions$umap@cell.embeddings[,3]



# Create plot with data  ################################################################
plot.data <- FetchData(obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Cell_cluster_label"))

plot.data$label <- paste(plot.data$Cell_cluster_label)

# Cool function that replicates ggplot's color scheme for any number of points  #########

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colors = gg_color_hue(N)  # Where N is the number of clusters you have 

# Plotting the data ####################################################################

p <- plot_ly(data = plot.data, 
             x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
             colors = colors,    
             color = ~Cell_cluster_label, 
             type = "scatter3d", 
             mode = "markers", 
             marker = list(size = 2, width=3), # controls size of points
             text=~label, #This is that extra column we made earlier for which we will use for cell ID
             hoverinfo="text")

l <- list(
  font = list(
    size = 14),
  itemsizing = "constant")
p <- p %>% layout(legend = l)

# Saving plot as HTML file ############################################################
htmlwidgets::saveWidget(as_widget(p), "3dplot.html")
