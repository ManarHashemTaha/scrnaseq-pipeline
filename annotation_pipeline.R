#!/usr/bin/env Rscript
# annotation_pipeline.R

# Usage:
#   In your shell, before running, increase stack and recursion:
#     ulimit -s unlimited
#   Then:
#     Rscript annotation_pipeline.R <ref_h5ad> <query_rds> <marker_csv>

# ---------------------------- #
# Increase recursion & C‐stack limits
# ---------------------------- #
options(expressions = 5e5)
Sys.setenv(R_CStackLimit = "10000000")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("❌ Usage: Rscript annotation_pipeline.R <ref_h5ad> <query_rds> <marker_csv>", call. = FALSE)
}
ref_h5ad    <- args[1]
query_rds   <- args[2]
marker_file <- args[3]

# ---------------------------- #
# Auto‐install missing Bioc packages
# ---------------------------- #
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (pkg in c("zellkonverter", "SingleCellExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
})

# ---------------------------- #
# Create output directories
# ---------------------------- #
suppressPackageStartupMessages(library(tools))
analysis_name <- file_path_sans_ext(basename(query_rds))
output_dir    <- file.path(getwd(), paste0(analysis_name, "_analysis"))
vis_dir       <- file.path(output_dir, "Visualization")
data_dir      <- file.path(output_dir, "Data")
dir.create(vis_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# Load libraries & parallelize
# ---------------------------- #
suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(Seurat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(dplyr)
})
workers <- parallel::detectCores() - 1
plan(multicore, workers = workers)
options(future.globals.maxSize = 4 * 1024^3)
message("🔧 Using ", workers, " cores; expressions=", getOption("expressions"), "; R_CStackLimit=", Sys.getenv("R_CStackLimit"))

# -------- 1. Read h5ad & convert --------
sce <- readH5AD(ref_h5ad)
symbols <- rowData(sce)$feature_name
names(symbols) <- rownames(sce)
ref_seurat <- as.Seurat(sce, counts = "X", data = NULL)
counts_matrix <- GetAssayData(ref_seurat, layer = "counts")
rownames(counts_matrix) <- symbols
ref <- CreateSeuratObject(counts = counts_matrix, meta.data = ref_seurat@meta.data)
DefaultAssay(ref) <- "RNA"

# -------- 2. Preprocess Reference --------
ref <- NormalizeData(ref, future.seed = TRUE)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
ref <- ScaleData(ref)
ref <- RunPCA(ref, features = VariableFeatures(ref), verbose = FALSE)
ref <- FindNeighbors(ref, dims = 1:10)
ref <- FindClusters(ref, resolution = 0.5)
ref <- RunUMAP(ref, dims = 1:10, return.model = TRUE)

if ("Cell_label" %in% colnames(ref@meta.data)) {
  jpeg(file.path(vis_dir, "ref_umap.jpg"), width = 10, height = 10, units = "in", res = 300)
  print(DimPlot(ref, reduction = "umap", group.by = "Cell_label", label = TRUE, repel = TRUE) + NoLegend())
  dev.off()
} else {
  warning("Column 'Cell_label' not found in reference metadata.")
}

# -------- 3. Label Transfer --------
Data_object <- readRDS(query_rds)
features <- intersect(rownames(ref), rownames(Data_object))
anchors_int <- FindTransferAnchors(
  reference           = ref,
  query               = Data_object,
  dims                = 1:50,
  features            = features,
  reference.reduction = "pca",
  normalization.method= "LogNormalize"
)
# pick safe k.weight
nAnch <- nrow(anchors_int@anchors)
kw    <- max(1, min(20, nAnch - 1))
message("🔧 Transfer k.weight = ", kw, " (", nAnch, " anchors)")
predictions <- TransferData(
  anchorset = anchors_int,
  refdata   = ref@meta.data$Cell_label,
  dims      = 1:20,
  k.weight  = kw
)
Data_object <- AddMetaData(Data_object, metadata = predictions)

jpeg(file.path(vis_dir, "query_predicted_labels.jpg"), width = 15, height = 10, units = "in", res = 300)
print(DimPlot(Data_object, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE))
dev.off()





saveRDS(
  Data_object,
  file = file.path(data_dir, "Data_object_int_annotated.rds")
)

cat("✅ Annotation pipeline complete. Outputs saved under:", output_dir, "\n")
