#!/usr/bin/env Rscript

library(tidyverse)
library(logger)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

df = read.csv(snakemake@input$genes, header = TRUE)
colnames(df) = c("cluster","marker")
markers_cells = read.csv(snakemake@input$markers, header = TRUE)

df = left_join(df,markers_cells[,c("cell","marker")]) %>%
  group_by(cluster,marker) %>%
  summarise(cells = paste0(unique(cell), collapse = ";"))

write.csv(df, snakemake@output$annotated_genes, row.names = FALSE)
