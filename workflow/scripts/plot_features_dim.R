
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

log_info("Reading files")
combined2 = readRDS(snakemake@input$rds)
genes = read.csv(snakemake@input$genes, header = TRUE) %>% pull(gene)

log_info("Generating FeaturePlot")
p = FeaturePlot(
  combined2,
  features = genes,
  ncol = 6,
)

log_info("Save plot")
ggsave(
  filename = snakemake@output[["featureplot"]],
  plot = p,
  width = 16,
  height = 12
)