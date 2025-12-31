#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)
library(tidyverse)
library(logger)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

log_info("Get files and parameters from Snakemake")
inputs <- snakemake@input
params <- snakemake@params

case_names <- params$case_names
control_names <- params$control_names
samples <- c(case_names, control_names)

read_data <- function(mtx, genes, barcodes) {
  ReadMtx(
    mtx = mtx,
    features = genes,
    cells = barcodes
  )
}

log_info("Reading count matrices")
cases <- mapply(
  read_data,
  inputs$cases_mtx,
  inputs$cases_genes,
  inputs$cases_barcodes,
  SIMPLIFY = FALSE
)
names(cases) <- case_names

controls <- mapply(
  read_data,
  inputs$controls_mtx,
  inputs$controls_genes,
  inputs$controls_barcodes,
  SIMPLIFY = FALSE
)
names(controls) <- control_names

data <- c(cases, controls)



log_info("Creating Seurat objects")
seurat_objs <- lapply(names(data), function(sample) {
  CreateSeuratObject(
    counts = data[[sample]],
    project = sample,
    min.cells = snakemake@params$min_cells,
    min.features = snakemake@params$min_features
  )
})
names(seurat_objs) <- names(data)

log_info("Adding group information")
for (sample in names(seurat_objs)) {
  if (sample %in% control_names) {
    seurat_objs[[sample]]$group <- "control"
  } else {
    seurat_objs[[sample]]$group <- "case"
  }
}

log_info("Normalizing and running PCA")
for (sample in names(seurat_objs)) {
  seurat_objs[[sample]] <- NormalizeData(
    seurat_objs[[sample]],
    normalization.method = snakemake@params$normalization_method,
    scale.factor = snakemake@params$scale_factor
  )

  seurat_objs[[sample]] <- FindVariableFeatures(
    seurat_objs[[sample]],
    selection.method = snakemake@params$selection_method,
    nfeatures = snakemake@params$nfeatures
  )
  seurat_objs[[sample]] <- ScaleData(seurat_objs[[sample]])
  seurat_objs[[sample]] <- RunPCA(
    seurat_objs[[sample]],
    features = VariableFeatures(seurat_objs[[sample]])
  )
}

log_info("Integrating datasets")
features <- SelectIntegrationFeatures(seurat_objs)
anchors <- FindIntegrationAnchors(
  object.list = seurat_objs,
  anchor.features = features
)
combined <- IntegrateData(anchors)

log_info("Get number of cells")
total_cells <- ncol(combined)

a = data.frame(table(combined@meta.data$orig.ident))
colnames(a) = c("orig.ident","n_cell")
cell_summary = left_join(a,unique(combined@meta.data[,c("orig.ident","group")]))

cell_summary <- cell_summary %>%
  mutate(
    percent_total = round(n_cell / total_cells * 100, 2),
    total_cells = total_cells
  ) %>%
  rename(
    sample = orig.ident
  )

log_info("Save cell summary")
write.csv(
  cell_summary,
  snakemake@output[["cell_counts"]],
  row.names = FALSE
)

log_info("Saving RDS file")
saveRDS(combined, snakemake@output$rds)
log_info("Integration finished successfully")
