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

log_info("Reading combined reduction RDS file")
combined2 = readRDS(snakemake@input$rds)

log_info(
  paste0("Generating DimPlot:", snakemake@params[['reduction']], ";", snakemake@params[['legend']])
  )

p = DimPlot(
  combined2,
  reduction = snakemake@params[["reduction"]],
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  group.by = snakemake@params[["legend"]] 
)

log_info("Save plot")
ggsave(
  filename = snakemake@output[["dimplot"]],
  plot = p,
  width = 6,
  height = 5
)