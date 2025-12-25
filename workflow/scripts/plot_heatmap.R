
#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(logger)

log_info("Reading combined reduction RDS file")
combined2 = readRDS(snakemake@input$rds)

log_info("Finding markers for heatmap")
markers <- FindAllMarkers(
  combined2,
  only.pos = TRUE,
  min.pct = snakemake@params[["min_pct"]],
  logfc.threshold = snakemake@params[["logfc_threshold"]]
)

log_info(paste0("Generating Heatmap with top ", snakemake@params[["n_markers"]], " markers per cluster"))
top <- markers %>%
  group_by(cluster) %>%
  top_n(snakemake@params[["n_markers"]], avg_log2FC)

log_info("Generating Heatmap plot")
p = DoHeatmap(
  combined2,
  features = top$gene
) + NoLegend()

log_info("Save plot")
ggsave(
  filename = snakemake@output[["heatmap"]],
  plot = p,
  width = 10,
  height = 8
)
