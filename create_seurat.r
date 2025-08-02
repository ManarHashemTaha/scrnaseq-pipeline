#!/usr/bin/env Rscript
# merge_seurat_with_metadata.R

# Usage:
# Rscript merge_seurat_with_metadata.R <input_dir> <metadata.csv> <output_rds>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("❌ Usage: Rscript merge_seurat_with_metadata.R <input_dir> <metadata.csv> <output_rds>", call. = FALSE)
}
input_path     <- args[1]
metadata_file  <- args[2]
output_file    <- args[3]

#-----------------------------#
# EDITED: Create output dirs  #
#-----------------------------#
library(tools)  
analysis_name <- basename(normalizePath(input_path))                            # ← EDITED
output_dir    <- file.path(getwd(), paste0(analysis_name, "_analysis"))         # ← EDITED
vis_dir       <- file.path(output_dir, "Visualization")                         # ← EDITED
data_dir      <- file.path(output_dir, "Data")                                  # ← EDITED
dir.create(vis_dir, recursive = TRUE, showWarnings = FALSE)                     # ← EDITED
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)                    # ← EDITED

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
})

# Read in each sample and build Seurat objects
files <- list.files(path = input_path)
seurat_objects <- list()
for (file in files) {
  countmatrix <- ReadMtx(
    mtx      = file.path(input_path, file, "matrix.mtx"),
    features = file.path(input_path, file, "genes.tsv"),
    cells    = file.path(input_path, file, "barcodes.tsv")
  )
  obj <- CreateSeuratObject(counts = countmatrix)
  obj$orig.ident <- substr(file, 1, 10)
  seurat_objects[[ substr(file, 1, 10) ]] <- obj
}

# Merge all samples
Data_object <- Reduce(function(x, y) merge(x, y), seurat_objects)

# Load and merge metadata
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
merged_metadata <- merge(
  Data_object@meta.data,
  metadata,
  by.x = "orig.ident",
  by.y = "sample_id",
  all.x = TRUE,
  sort  = FALSE
)
rownames(merged_metadata) <- rownames(Data_object@meta.data)
Data_object@meta.data <- merged_metadata

#-----------------------------#
# EDITED: Save outputs       #
#-----------------------------#
# Save the merged Seurat object into the Data folder
saveRDS(
  Data_object,
  file = file.path(data_dir, basename(output_file))                         # ← EDITED
)

cat("✅ Seurat object created with metadata and saved to:", file.path(data_dir, basename(output_file)), "\n")

