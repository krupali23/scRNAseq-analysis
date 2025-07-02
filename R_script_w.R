install.packages("Seurat")
install.packages('Seurat')
library(Seurat)
library(Matrix)
library(hdf5r)
library(ggplot2)
#install.packages("hdf5r")
data <- read.csv("C:\\Users\\krupa\\Desktop\\Wafa_thyroid\\Data\\GSE183963_Single_cell_mouse_ndata.csv", header = TRUE, row.names = 1, sep = ";")
# Replace "path_to/GSE183963_Cellseq2barcodes_single_cell.csv" with the actual path to your cell barcode data
#cell_barcodes <- read.csv("C:\\Users\\krupa\\Desktop\\Wafa_thyroid\\Data\\GSE183963_Cellseq2barcodes_single_cell.csv", header = FALSE)
seurat_obj <- CreateSeuratObject(counts = data, project = "project_name")
my_matrix <- as.matrix(data)
data_sparse <- as(my_matrix, "dgCMatrix")
# Create a Seurat object for the first dataset
seurat_obj <- CreateSeuratObject(counts = data, project = "GSE183963")
# Check the number of barcodes and columns in your Seurat object
n_barcodes <- nrow(cell_barcodes)
n_columns <- ncol(seurat_obj)
# Split the V1 column based on the common prefix
cell_barcodes_split <- data.frame(do.call(rbind, strsplit(as.character(cell_barcodes$V1), "\t")))

# Ensure that the result has two columns
if (ncol(cell_barcodes_split) == 2) {
  colnames(cell_barcodes_split) <- c("Index", "Identifier")
} else {
  stop("Unexpected structure after splitting.")
}

# Check the structure of cell_barcodes_split
str(cell_barcodes_split)

# Check if common identifiers are now present
common_identifiers <- intersect(cell_barcodes_split$Identifier, substr(colnames(seurat_obj), nchar("mThyroidtissue1_") + 1, nchar(colnames(seurat_obj))))

# Check if common barcodes are now present
common_barcodes <- paste("mThyroidtissue1_", common_identifiers, sep = "")


# Create a metadata data frame with only common barcodes
metadata_df_common <- metadata_df[metadata_df$barcodes %in% common_barcodes, ]

# Assign metadata to the Seurat object
seurat_obj$meta.data <- metadata_df_common
  
# Continue with any additional processing steps as needed
# Add percent.mt information
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Normalize data
seurat_obj <- NormalizeData(seurat_obj)
# Identify variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Visualize QC metrics as a violin plot
#VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(all, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.2)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
#Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
dev.off()
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
VlnPlot(pbmc, features = c("Trpm5", "Pou2f3", "Chat", "Gnat3", "Dclk1"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("Trpm5", "Pou2f3", "Chat", "Gnat3", "Dclk1"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("Trpm5", "Pou2f3", "Chat", "Gnat3", "Dclk1"))
VlnPlot(pbmc, features = c("Lrmp", "Bsnd", "Foxi1", "Ascl2", "Ascl3"), slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c("Lrmp", "Bsnd", "Foxi1", "Ascl2", "Ascl3"))
VlnPlot(pbmc, features = c("Dclk1", "Bsnd", "Foxi1", "Ascl2", "Ascl3", "Ascl1", "Trpm5", "Uchl1", "Gnat3"))
dev.off()
table(pbmc$seurat_clusters)
# Subsetting cells from a specific cluster (replace 'cluster_id' with the actual cluster you want to subcluster)
sub_cluster_id <- 12  # Replace with the cluster you want to subcluster
sub_cluster_cells <- WhichCells(pbmc, ident = sub_cluster_id)

# Create a new Seurat object with only the cells from the selected cluster
seurat_subcluster <- subset(pbmc, cells = sub_cluster_cells)

# Preprocess the subclustered Seurat object
seurat_subcluster <- NormalizeData(seurat_subcluster)
seurat_subcluster <- FindVariableFeatures(seurat_subcluster, selection.method = "vst", nfeatures = 2000)
seurat_subcluster <- ScaleData(seurat_subcluster)

# Run PCA and cluster the cells within the subcluster
seurat_subcluster <- RunPCA(seurat_subcluster, verbose = FALSE)
seurat_subcluster <- FindNeighbors(seurat_subcluster, dims = 1:30)
seurat_subcluster <- FindClusters(seurat_subcluster, resolution = 0.2)
seurat_subcluster <- RunUMAP(seurat_subcluster, dims = 1:10)
# Visualize the subclustering results
DimPlot(seurat_subcluster, group.by = "seurat_clusters")  # Adjust group.by based on your cluster column name
VlnPlot(seurat_subcluster, features = c("Dclk1", "Trpm5",  "Chat", "Pou2f3", "Vil1", "Avil", "Svil", "Gnat3", "Plcb2"))


table(seurat_subcluster$seurat_clusters)
VlnPlot(seurat_subcluster, features = c("Dclk1", "Bsnd", "Foxi1", "Tas1r1", "Tas1r2", "Tas1r3", "Trpm5", "Uchl1", "Chat", "Pou2f3"))
table(seurat_subcluster$seurat_clusters)