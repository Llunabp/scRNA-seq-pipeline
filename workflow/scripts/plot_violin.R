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

log_info("Extract data")
df <- FetchData(
  combined,
  vars = c("nFeature_RNA", "nCount_RNA", "group")
)

log_info("Plot nFeature")
p1 = ggplot(df, aes(x = group, y = nFeature_RNA, fill = group)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  theme_classic()

log_info("Plot nCount")
p2 = ggplot(df, aes(x = group, y = nCount_RNA, fill = group)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  theme_classic()

log_info("Save plots")
ggsave(snakemake@output$plot_feature,p1)
ggsave(snakemake@output$plot_count,p2)