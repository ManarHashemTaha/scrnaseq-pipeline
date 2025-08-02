#!/usr/bin/env Rscript
# automated_cellchat_analysis.R

# -------------------------------------------------------------------------
# Parallel CellChat analysis pipeline with comprehensive visualizations (JPEG)
# -------------------------------------------------------------------------

# 1. Install & load dependencies
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_github("sqjin/CellChat", upgrade="never")
cran_pkgs <- c("Seurat","patchwork","future","future.apply","dplyr","ggplot2","ComplexHeatmap")
for (pkg in cran_pkgs) if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(tools); library(parallel); library(future); library(future.apply)
  library(Seurat); library(SeuratObject); library(CellChat)
  library(ComplexHeatmap); library(dplyr); library(patchwork); library(ggplot2); library(grid)
})

# 2. Parse arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) stop("Usage: Rscript automated_cellchat_analysis.R <rds_dir> [group_col] [annot_col]")
rds_dir   <- args[1]
group_col <- ifelse(length(args)>=2, args[2], "stage")
annot_col <- ifelse(length(args)>=3, args[3], "predicted.id")

# 3. Output directories
analysis_name <- basename(rds_dir)
out_dir  <- file.path(getwd(), paste0(analysis_name, "_cellchat_analysis"))
vis_dir  <- file.path(out_dir, "Visualization")
data_dir <- file.path(out_dir, "Data")
for (d in c(out_dir, vis_dir, data_dir)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

# 4. Parallel config
n_cores <- max(1, detectCores()-1)
plan(multisession, workers=n_cores)
options(future.globals.maxSize=8*1024^3)
message("Using ", n_cores, " workers")

# 5. Helper: build CellChat from flat matrix
create_cc_obj <- function(path) {
  sample <- file_path_sans_ext(basename(path))
  message("[", sample, "] loading...")
  seu <- readRDS(path)

    # Flatten RNA assay layers into a single layer
  DefaultAssay(seu) <- "RNA"
  # Merge counts, data, and scale.data layers into one
  seu <- JoinLayers(seu, assay = "RNA", layers = c("counts", "data", "scale.data"))
  # Extract combined counts matrix
  expr_mat <- as.matrix(LayerData(seu, assay = "RNA", layer = "counts"))

  # Clean metadata
  seu@meta.data[[annot_col]] <- droplevels(factor(seu@meta.data[[annot_col]]))

  # Build CellChat from matrix
  cc <- createCellChat(
    object   = expr_mat,
    meta     = seu@meta.data,
    group.by = group_col
  )
  cc <- setIdent(cc, ident.use = annot_col)
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling", key = "annotation")
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  saveRDS(cc, file.path(data_dir, paste0(sample, "_cc.rds")))
  return(cc)
}

# 6. Process files Process files Process files
rds_files <- list.files(rds_dir, pattern="\\.rds$", full.names=TRUE)
cc_list <- future_lapply(rds_files, create_cc_obj, future.seed=TRUE)
names(cc_list) <- file_path_sans_ext(basename(rds_files))

# 7. Merge
# Filter out any failed CellChat objects
valid_idx <- sapply(cc_list, function(x) inherits(x, "CellChat"))
if (any(!valid_idx)) {
  warning("Dropping ", sum(!valid_idx), " failed CellChat objects before merge.")
  cc_list <- cc_list[valid_idx]
}
if (length(cc_list) < 2) stop("Need at least two valid CellChat objects to merge.")
common <- Reduce(intersect, lapply(cc_list, function(z) levels(z@idents)))
cc_list <- lapply(cc_list, function(z){ z@idents <- factor(z@idents, levels=common); z })
merged_cc <- mergeCellChat(cc_list, add.names=names(cc_list), cell.prefix=TRUE, merge.data=FALSE)
saveRDS(cc_list, file.path(data_dir, "cc_list.rds"))
saveRDS(merged_cc, file.path(data_dir, "cc_merged.rds")))

# 8. Visualizations (JPEG) ------------------------------------------------------------
# 8.1 compareInteractions
jpeg(file.path(vis_dir,"compareInteractions.jpg"), width=8, height=6, units="in", res=300)
print(compareInteractions(merged_cc, group=c(1,2)))
print(compareInteractions(merged_cc, group=c(1,2), measure="weight"))
dev.off()

# 8.2 differential interactions
jpeg(file.path(vis_dir,"diffInteractions_circle.jpg"), width=15, height=10, units="in", res=300)
par(mfrow=c(1,2), xpd=TRUE)
netVisual_diffInteraction(merged_cc, weight.scale=TRUE)
netVisual_diffInteraction(merged_cc, weight.scale=TRUE, measure="weight")
dev.off()

jpeg(file.path(vis_dir,"diffInteractions_heatmap.jpg"), width=10, height=7, units="in", res=300)
g1<-netVisual_heatmap(merged_cc)
g2<-netVisual_heatmap(merged_cc, measure="weight")
print(g1+g2)
dev.off()

# 8.3 circle plot per dataset
jpeg(file.path(vis_dir,"splitted_diffCircle.jpg"), width=10, height=7, units="in", res=300)
wmax<-getMaxWeight(cc_list, attribute="count")
par(mfrow=c(1,2), xpd=TRUE)
for(i in seq_along(cc_list)){
  netVisual_circle(cc_list[[i]]@net$count, weight.scale=TRUE, edge.width.max=12,
                   title.name=names(cc_list)[i])
}
dev.off()

# 8.4 signaling role scatter
ggs <- lapply(cc_list, netAnalysis_computeCentrality)
gg <- lapply(seq_along(ggs), function(i) netAnalysis_signalingRole_scatter(ggs[[i]], title=names(ggs)[i]))
jpeg(file.path(vis_dir,"signalingRole_scatter.jpg"), width=10, height=6, units="in", res=300)
print(wrap_plots(gg))
dev.off()

# 8.5 rankNet
jpeg(file.path(vis_dir,"rankNet.jpg"), width=8, height=6, units="in", res=300)
g1 <- rankNet(merged_cc, mode="comparison", measure="weight", stacked=TRUE, do.stat=TRUE)
g2 <- rankNet(merged_cc, mode="comparison", measure="weight", stacked=FALSE, do.stat=TRUE)
print(g1+g2)
dev.off()

# 8.6 signaling role heatmap
i1<-netAnalysis_signalingRole_heatmap(ggs[[1]], pattern="outgoing", signaling=union(ggs[[1]]@netP$pathways,ggs[[2]]@netP$pathways))
i2<-netAnalysis_signalingRole_heatmap(ggs[[2]], pattern="outgoing", signaling=union(ggs[[1]]@netP$pathways,ggs[[2]]@netP$pathways))
jpeg(file.path(vis_dir,"role_outgoing_heatmap.jpg"), width=10, height=10, units="in", res=300)
ComplexHeatmap::draw(i1 + i2, ht_gap=unit(0.5,'cm'))
dev.off()

# 8.7 netVisual_bubble for all idents
i<-1; targets<-2:length(cc_list)
jpeg(file.path(vis_dir,'bubble_overall.jpg'), width=8, height=6, units="in", res=300)
g1<-netVisual_bubble(merged_cc, sources.use=i, targets.use=targets, comparison=c(1,2), angle.x=45)
g2<-netVisual_bubble(merged_cc, sources.use=i, targets.use=targets, comparison=c(1,2), max.dataset=2, title.name='Incr', angle.x=45)
print(g1+g2)
dev.off()

# 8.8 gene expression violin per pathway
datasets<-names(cc_list)
merged_cc@meta$datasets <- factor(merged_cc@meta$datasets, levels=datasets)
for(pw in c("CCL","MIF","GALECTIN","BAFF","MK","IL16","FLT3")){
  jpeg(file.path(vis_dir,paste0('expr_',pw,'.jpg')), width=8, height=6, units="in", res=300)
  print(plotGeneExpression(merged_cc, signaling=pw, split.by='datasets', colors.ggplot=TRUE, type='violin'))
  dev.off()
}

message("✅ Done: visuals in ", vis_dir)

