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

log_info("Read RDS")
combined = readRDS(snakemake@input$rds)
DefaultAssay(combined) = "RNA"

log_info("Find top features")
top<- head(VariableFeatures(combined), snakemake@params$top_features)

log_info("Plot most variable features")
p1 = VariableFeaturePlot(combined) 
p2 = LabelPoints(plot = p1, points = top, repel = TRUE)

log_info("Save plot")
ggsave(snakemake@output$plot,p2)
