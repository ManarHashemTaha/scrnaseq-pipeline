
# automated_integration_rpca.R

# ---------------------------- #
# Load Required Libraries
# ---------------------------- #
library(Seurat)
library(ggplot2)
library(patchwork)

# ---------------------------- #
# Load Input Object from Command Line
# ---------------------------- #
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("❌ Please provide path to a Seurat .rds file as input.", call. = FALSE)
}

input_path <- args[1]
cat("📥 Loading Seurat object from:", input_path, "\n")
MM_obj_subset <- readRDS(input_path)

# ---------------------------- #
# Split, Normalize, and PCA
# ---------------------------- #
MM_obj_subset.list <- SplitObject(MM_obj_subset, split.by = "orig.ident")

MM_obj_subset.list <- lapply(X = MM_obj_subset.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})

features <- SelectIntegrationFeatures(object.list = MM_obj_subset.list)

MM_obj_subset.list <- lapply(X = MM_obj_subset.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# ---------------------------- #
# Integration with RPCA
# ---------------------------- #
options(future.globals.maxSize = 2 * 1024^3)  # 2GB

MM.anchors <- FindIntegrationAnchors(object.list = MM_obj_subset.list, anchor.features = features, reduction = "rpca")
MM.combined <- IntegrateData(anchorset = MM.anchors)

DefaultAssay(MM.combined) <- "integrated"

# ---------------------------- #
# Downstream Processing
# ---------------------------- #
MM.combined <- ScaleData(MM.combined, verbose = FALSE)
MM.combined <- RunPCA(MM.combined, npcs = 30, verbose = FALSE)
MM.combined <- RunUMAP(MM.combined, reduction = "pca", dims = 1:15)
MM.combined <- FindNeighbors(MM.combined, reduction = "pca", dims = 1:15)

MM.combined <- FindClusters(MM.combined, resolution = 0.5)

# ---------------------------- #
# Save UMAP Visualizations
# ---------------------------- #
pdf("UMAP_rpca.pdf", width=10, height=10)
DimPlot(MM.combined, reduction = "umap", group.by = "orig.ident")
dev.off()

pdf("UMAP_rpca_2.pdf", width=10, height=10)
DimPlot(MM.combined, reduction = "umap")
dev.off()

pdf("UMAP_rpca-splited.IDs-2.pdf", width=10, height=10)
DimPlot(MM.combined, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident")
dev.off()

# ---------------------------- #
# Save Integrated Object
# ---------------------------- #
saveRDS(MM.combined, file = "integrated_MM_object.rds")
cat("✅ Integration complete. Output saved to: integrated_MM_object.rds\n")
