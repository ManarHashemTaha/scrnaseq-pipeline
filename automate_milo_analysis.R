#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# automate_milo_analysis.R
#
# Parallel Milo neighborhood analysis pipeline over saved Seurat RDS files
# from an earlier split-run script. Reads each *_*.rds, preprocesses,
# converts to SingleCellExperiment (RNA assay only), and runs Milo.
#
# Usage:
#   Rscript automate_milo_analysis.R <rds_directory> [group_column] [annotation_column]
#
# <rds_directory>: folder containing Seurat .rds files (e.g., MM_Healthy.rds,...)
# [group_column]: metadata column for group labels (default: "status")
# [annotation_column]: metadata column for cell annotations (default: "NewAnnotation")
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tools)
  library(parallel)
  library(Seurat)
  library(SingleCellExperiment)
  library(miloR)
  library(scater)
  library(scran)
  library(scuttle)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript automate_milo_analysis.R <rds_directory> [group_col] [annot_col]")
}
rds_dir    <- args[1]
group_col  <- ifelse(length(args)>=2, args[2], "stage")
annot_col  <- ifelse(length(args)>=3, args[3], "predicted.id")

# Output directories
out_dir  <- file.path(rds_dir, "milo_analysis")
vis_dir  <- file.path(out_dir, "Visualization")
data_dir <- file.path(out_dir, "Data")
dir.create(vis_dir,  recursive=TRUE, showWarnings=FALSE)
dir.create(data_dir, recursive=TRUE, showWarnings=FALSE)

# List all Seurat RDS files in the directory
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No .rds files found in ", rds_dir)

# Number of cores
total_cores <- max(1, detectCores() - 1)
message("Using ", total_cores, " cores")

# Function to process one Seurat RDS
process_rds <- function(rds_path) {
  prefix <- file_path_sans_ext(basename(rds_path))
  message("-- Processing: ", prefix)

  # Load Seurat object
  seu <- readRDS(rds_path)

  # Preprocessing on RNA assay
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=5000)
  hvgs <- VariableFeatures(seu)
  seu <- ScaleData(seu, features=hvgs)
  seu <- RunPCA(seu, features=hvgs)
  seu <- FindNeighbors(seu, dims=1:10)
  seu <- FindClusters(seu, resolution=0.5)
  seu <- RunUMAP(seu, dims=1:10)

      # Manually build a pure SCE to avoid multi-layer issues
    counts_mat <- as.matrix(seu@assays$RNA@counts)
    logc_mat   <- as.matrix(seu@assays$RNA@data)
    sce <- SingleCellExperiment(
      assays = list(
        counts = counts_mat,
        logcounts = logc_mat
      ),
      colData = seu@meta.data
    )
    reducedDim(sce, "PCA")  <- Embeddings(seu, "pca")
    reducedDim(sce, "UMAP") <- Embeddings(seu, "umap")
    saveRDS(sce, file.path(data_dir, paste0(prefix, "_sce.rds")))

  # Milo pipeline
  milo_obj <- Milo(sce)
  milo_obj <- buildGraph(milo_obj, k=30, d=30, reduced.dim="PCA")
  milo_obj <- makeNhoods(milo_obj, prop=0.5, k=10, d=30,
                         refined=TRUE, refinement_scheme="graph")
  milo_obj <- countCells(milo_obj,
                         meta.data=as.data.frame(colData(milo_obj)),
                         samples="sample_id")

  # Design matrix
design_df <- as.data.frame(colData(milo_obj))[, c("sample_id", group_col, annot_col)]
dr <- distinct(design_df)
  rownames(dr) <- dr$sample_id
  dr <- dr[colnames(nhoodCounts(milo_obj)), ]

  # Differential testing
  # contrast: second level vs first level
  lvls <- unique(dr[[group_col]])
  contrast <- paste0(group_col, lvls[2], " - ", group_col, lvls[1])
  da_res <- testNhoods(milo_obj,
                       design=as.formula(paste0("~0+", group_col)),
                       design.df=dr,
                       model.contrasts=contrast,
                       fdr.weighting="graph-overlap")
  write.csv2(da_res, file.path(data_dir, paste0(prefix, "_DA_results.csv")))

  # Visualization: p-value histogram
  jpeg(file.path(vis_dir, paste0(prefix, "_pvalue_hist.jpg")), width=8, height=6, res=300)
  print(ggplot(da_res, aes(PValue)) + geom_histogram(bins=50))
  dev.off()

  # Annotate & group results
  da_annot <- annotateNhoods(milo_obj, da_res, coldata_col=annot_col)
  da_annot$cellType <- ifelse(
    da_annot[[paste0(annot_col, "_fraction")]] < 0.7,
    "Mixed", da_annot[[annot_col]]
  )
  write.csv2(da_annot, file.path(data_dir, paste0(prefix, "_DA_grouped.csv")))

  # Beeswarm plot
  jpeg(file.path(vis_dir, paste0(prefix, "_beeswarm.jpg")), width=6, height=6, res=300)
  plotDAbeeswarm(da_annot, group.by=annot_col)
  dev.off()

  # Neighborhood graph on UMAP
  milo_obj <- buildNhoodGraph(milo_obj)
  jpeg(file.path(vis_dir, paste0(prefix, "_umap_nhood.jpg")), width=10, height=10, res=300)
  p1 <- plotReducedDim(milo_obj, dimred="UMAP", colour_by=group_col) + guides(fill="none")
  p2 <- plotNhoodGraphDA(milo_obj, da_annot, layout="UMAP", alpha=0.1)
  print(p1 + p2)
  dev.off()

  # Save Milo object
  saveRDS(milo_obj, file.path(data_dir, paste0(prefix, "_milo.rds")))
}

# Run in parallel over RDS files
mclapply(rds_files, process_rds, mc.cores = ncores <- total_cores)

message("All comparisons processed.")

