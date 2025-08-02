#!/usr/bin/env Rscript
# automated_seurat_analysis.R

# Usage:
# Rscript automated_seurat_analysis.R F:/Automation/Data-object.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("❌ Please provide path to a Seurat .rds file as input.", call. = FALSE)
}
input_path <- args[1]
cat("📥 Loading Seurat object from:", input_path, "\n")
Data_object <- readRDS(input_path)

#-----------------------------#
# EDITED: Create output dirs #
#-----------------------------#
library(tools)
analysis_name <- file_path_sans_ext(basename(input_path))
output_dir    <- file.path(getwd(), paste0(analysis_name, "_analysis"))
vis_dir       <- file.path(output_dir, "Visualization")
data_dir      <- file.path(output_dir, "Data")
dir.create(vis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------#
# EDITED: Parallel setup     #
#-----------------------------#
suppressPackageStartupMessages({
  library(future)
  library(RhpcBLASctl)
})
workers <- parallel::detectCores() - 1
options(future.globals.maxSize = 8 * 1024^3)
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = workers)
} else {
  plan(multicore,   workers = workers)
}
blas_set_num_threads(8)
omp_set_num_threads(8)
cat("🔧 Parallelization initialized with", workers, "workers\n")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
})

#-----------------------------#
# Initial Setup & QC
#-----------------------------#
Idents(Data_object) <- Data_object@meta.data$orig.ident
Data_object[["percent.mt"]] <- PercentageFeatureSet(Data_object, pattern = "^MT-")

# QC plots → JPGs in Visualization
ggsave(
  filename = file.path(vis_dir, "plot_features.jpg"),
  plot     = VlnPlot(Data_object, features = "nFeature_RNA", pt.size = 0)
)
ggsave(
  filename = file.path(vis_dir, "plot_nCount.jpg"),
  plot     = VlnPlot(Data_object, features = "nCount_RNA", pt.size = 0)
)
ggsave(
  filename = file.path(vis_dir, "plot_percent_mt.jpg"),
  plot     = VlnPlot(Data_object, features = "percent.mt", pt.size = 0)
)

Data_object <- subset(
  Data_object,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)

# After‐filtering QC plots → JPGs
ggsave(
  filename = file.path(vis_dir, "plot_features_after.jpg"),
  plot     = VlnPlot(Data_object, features = "nFeature_RNA", pt.size = 0)
)
ggsave(
  filename = file.path(vis_dir, "plot_nCount_after.jpg"),
  plot     = VlnPlot(Data_object, features = "nCount_RNA", pt.size = 0)
)
ggsave(
  filename = file.path(vis_dir, "plot_percent_mt_after.jpg"),
  plot     = VlnPlot(Data_object, features = "percent.mt", pt.size = 0)
)

#-----------------------------#
# Normalize & Feature Selection
#-----------------------------#
Data_object <- NormalizeData(Data_object)
Data_object <- FindVariableFeatures(Data_object, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Data_object), 10)

plot1 <- VariableFeaturePlot(Data_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

ggsave(
  filename = file.path(vis_dir, "variable_features.jpg"),
  plot     = plot1
)
ggsave(
  filename = file.path(vis_dir, "top10_variable_features.jpg"),
  plot     = plot2
)

#-----------------------------#
# Scaling and PCA
#-----------------------------#
all.genes <- rownames(Data_object)
Data_object <- ScaleData(Data_object)
Data_object <- RunPCA(Data_object, features = VariableFeatures(object = Data_object))

jpeg(file.path(vis_dir, "pca_summary.jpg"), width = 800, height = 600)
print(Data_object[["pca"]], dims = 1:5, nfeatures = 5)
dev.off()

jpeg(file.path(vis_dir, "pca_dim_loadings.jpg"), width = 800, height = 600)
VizDimLoadings(Data_object, dims = 1:2, reduction = "pca")
dev.off()

jpeg(file.path(vis_dir, "pca_dimplot.jpg"), width = 800, height = 600)
DimPlot(Data_object, reduction = "pca")
dev.off()

jpeg(file.path(vis_dir, "pca_heatmap1.jpg"), width = 800, height = 600)
DimHeatmap(Data_object, dims = 1, cells = 500, balanced = TRUE)
dev.off()

jpeg(file.path(vis_dir, "pca_heatmap_1to15.jpg"), width = 800, height = 600)
DimHeatmap(Data_object, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#-----------------------------#
# Dimensionality Estimation
#-----------------------------#
Data_object <- JackStraw(Data_object, num.replicate = 100)
Data_object <- ScoreJackStraw(Data_object, dims = 1:20)

jpeg(file.path(vis_dir, "jackstraw_plot.jpg"), width = 800, height = 600)
JackStrawPlot(Data_object, dims = 1:15)
dev.off()

jpeg(file.path(vis_dir, "elbow_plot.jpg"), width = 800, height = 600)
ElbowPlot(Data_object)
dev.off()

#-----------------------------#
# Clustering and UMAP
#-----------------------------#
Data_object <- FindNeighbors(Data_object, dims = 1:10)
Data_object <- FindClusters(Data_object, resolution = 1)
Data_object <- FindClusters(Data_object, resolution = 0.5)
Data_object <- RunUMAP(Data_object, dims = 1:10)

jpeg(file.path(vis_dir, "umap_clustered.jpg"), width = 800, height = 600)
DimPlot(Data_object, reduction = "umap")
dev.off()

jpeg(file.path(vis_dir, "umap_by_sample.jpg"), width = 800, height = 600)
DimPlot(Data_object, reduction = "umap", label = FALSE, group.by = "orig.ident")
dev.off()

jpeg(file.path(vis_dir, "umap_split_sample.jpg"), width = 800, height = 600)
DimPlot(Data_object, reduction = "umap", label = FALSE, group.by = "orig.ident", split.by = "orig.ident")
dev.off()

#-----------------------------#
# Marker Discovery
#-----------------------------#
Data_object <- JoinLayers(Data_object)
Data_object.markers <- FindAllMarkers(Data_object, only.pos = TRUE)

marker <- Data_object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(marker,
          file = file.path(data_dir, "top2_cluster_markers.csv"),
          row.names = TRUE)

cluster0.markers <- FindMarkers(Data_object, ident.1 = 0,
                                logfc.threshold = 0.25,
                                test.use = "roc", only.pos = TRUE)
write.csv(cluster0.markers,
          file = file.path(data_dir, "cluster0_markers.csv"),
          row.names = TRUE)

jpeg(file.path(vis_dir, "marker_FeaturePlot_set1.jpg"), width = 800, height = 600)
FeaturePlot(Data_object, features = c("BLOC1S5-TXNDC5", "IL7R", "NKG7", "IGLL5", "MT-ATP8", "H3-3B"))
dev.off()

jpeg(file.path(vis_dir, "marker_FeaturePlot_set2.jpg"), width = 800, height = 600)
FeaturePlot(Data_object, features = c("LYZ", "S100A8", "H4C3", "TCL1A", "HBB", "HSPB1"))
dev.off()

jpeg(file.path(vis_dir, "marker_FeaturePlot_set3.jpg"), width = 800, height = 600)
FeaturePlot(Data_object, features = c("CD79B", "MT-ND5", "HBA2", "STMN1", "G0S2", "MT-ND4L", "JCHAIN"))
dev.off()


#-----------------------------#
# Save Final Seurat Object
#-----------------------------#
saveRDS(Data_object,
     file = file.path(data_dir, "Data_object_Clustered.RDS"))

cat("✅ All outputs saved under:", output_dir, "\n")

