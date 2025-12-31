#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(logger)


# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")


combined = readRDS(snakemake@input$rds)

combined <- ScaleData(combined)
combined <- RunPCA(combined)

p = ElbowPlot(combined, ndims = snakemake@params[["npcs"]])

ggsave(
  filename = snakemake@output[["elbowplot"]],
  plot = p,
  width = 6,
  height = 4
)


sd <- combined@reductions$pca@stdev  
diffs <- abs(diff(sd))
components = which.min(diffs[1:30])



set.seed(snakemake@params[["seed"]])

combined2 <- RunTSNE(combined, dims = 1:components)
combined2 <- FindNeighbors(combined2, dims = 1:components)
combined2 <- FindClusters(combined2, 
    resolution = snakemake@params[["resolution"]], 
    algorithm = snakemake@params[["algorithm"]])

combined2 <- RunUMAP(
  combined2,
  dims = 1:components,
  reduction = "pca",
  n.neighbors = snakemake@params[["n_neighbors"]],
  min.dist = snakemake@params[["min_dist"]]
)

markers <- FindAllMarkers(
  combined2,
  only.pos = TRUE,
  min.pct = snakemake@params[["min_pct"]],
  logfc.threshold = snakemake@params[["logfc_threshold"]],
  test.use = snakemake@params[["test_use"]]
)

labels <- markers %>% 
  group_by(cluster) %>%
  top_n(1, avg_log2FC) %>% 
  arrange(cluster) %>% 
  mutate(name = paste(cluster,gene,sep = ": "))

ids = labels$name
names(ids) = labels$cluster


write.csv(labels[,c("cluster","gene")],snakemake@output$markers,
          row.names = FALSE)


metadata = combined2@meta.data %>%
  mutate(long_name = ids[seurat_clusters])

combined2@meta.data$long_name = metadata$long_name

Idents(combined2) = metadata$seurat_clusters

saveRDS(combined2, file = snakemake@output$reduction)
