#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(logger)
library(clusterProfiler)


# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")
set.seed(snakemake@params[["seed"]])

run_gsea_cluster <- function(
    de_table,        
    cluster_name,    
    gene_set,        
    min_genes = snakemake@params[["min_genes"]]
) {
  
  # filter by min genes
  if (nrow(de_table) < min_genes) return(NULL)
  
  # ranking each cluster
  gene_rank <- de_table$avg_log2FC
  names(gene_rank) <- de_table$gene
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  
  # GSEA
  gsea <- tryCatch(
    GSEA(
      geneList = gene_rank,
      TERM2GENE = data.frame(
        term = "Signature_genes",
        gene = gene_set
      ),
      pvalueCutoff = 1,
      verbose = FALSE,
      pAdjustMethod = snakemake@params[["GSEA_padj_method"]],
      seed = snakemake@params[["seed"]]
    ),
    error = function(e) NULL
  )
  
  if (is.null(gsea) || nrow(gsea@result) == 0) return(NULL)
  
  
  data.frame(
    cluster = cluster_name,
    ES = gsea@result$enrichmentScore,
    NES = gsea@result$NES,
    P.value = gsea@result$pvalue,
    FDR = gsea@result$p.adjust
  )
}


combined2 = readRDS(snakemake@input$rds)
clusters <- levels(Idents(combined2))
a = combined2@meta.data %>% select(seurat_clusters,long_name) %>% unique() %>% 
  rename(cluster = seurat_clusters,
         Cluster = long_name)

DE_summary = data.frame(cluster = clusters)
DE_summary = left_join(DE_summary,a, by = "cluster")

DE_by_cluster <- list()
for (cl in clusters) {
  log_info(paste0("Cluster ", cl))
  cells_cl <- WhichCells(combined2, idents = cl)
  meta_cl <- combined2@meta.data[cells_cl, ]
  
  # comprobar que hay ambos grupos
  if (length(unique(meta_cl$group)) < 2) {
    log_info("  -> Skipping cluster, only one group present")
    next
  }
  
  de <- FindMarkers(
    combined2,
    ident.1 = "case",
    ident.2 = "control",
    group.by = "group",
    subset.ident = cl,
    test.use = snakemake@params[["test_use"]],
    logfc.threshold = snakemake@params[["logfc_threshold"]],
    min.pct = snakemake@params[["min_pct"]]
  )
  
  de$gene <- rownames(de)
  de$cluster <- cl
  
  de_sig <- de %>%
    filter(p_val_adj < snakemake@params[["FDR"]] & abs(avg_log2FC) > snakemake@params[["FC"]])
  
  DE_by_cluster[[cl]] <- de_sig
}


DE_all <- bind_rows(DE_by_cluster)
DE_all_sum <- DE_all %>%
  group_by(cluster) %>%
  summarise(
    DE_genes = n(),
    Up = sum(avg_log2FC > 0),
    Down = sum(avg_log2FC < 0)
  )

DE_summary = left_join(DE_summary,DE_all_sum, by = "cluster")


cell_counts <- combined2@meta.data %>%
  count(orig.ident, seurat_clusters, group)

total_cells_sample <- combined2@meta.data %>%
  count(orig.ident) %>%
  rename(total = n)

cell_percent_sample <- cell_counts %>%
  left_join(total_cells_sample, by = "orig.ident") %>%
  mutate(percent = n / total * 100)


cell_percent_summary <- cell_percent_sample %>%
  group_by(seurat_clusters, group) %>%
  summarise(
    mean_percent = mean(percent),
    sd_percent   = sd(percent),
    .groups = "drop"
  )

cell_percent_wide <- cell_percent_summary %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = c(mean_percent, sd_percent)
  )

cell_percent_wide <- cell_percent_wide %>%
  mutate(
    FoldChange = log2(mean_percent_case / mean_percent_control),
    cluster = seurat_clusters,
    "Control_%_cells" = paste0(mean_percent_control, " ± ", sd_percent_control),
    "Case_%_cells" = paste0(mean_percent_case, " ± ", sd_percent_case)
  ) %>% select(cluster,"Control_%_cells","Case_%_cells",FoldChange)

DE_summary = left_join(DE_summary, cell_percent_wide, by = "cluster")

pvals_cluster <- cell_percent_sample %>%
  group_by(seurat_clusters) %>%
  summarise(
    p_ttest = tryCatch(
      t.test(
        percent[group == "case"],
        percent[group == "control"]
      )$p.value,
      error = function(e) NA
    ),
    p_wilcox = tryCatch(
      wilcox.test(
        percent[group == "case"],
        percent[group == "control"]
      )$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  )

colnames(pvals_cluster) = c("cluster","P-value (t.test)", "P-value (wilcox)")

DE_summary = left_join(DE_summary, pvals_cluster, by = "cluster") 

DE_summary = DE_summary %>% select(!cluster)

write.csv(DE_summary, snakemake@output[["de_summary"]], row.names = FALSE)


DE_genes <- scan(
  snakemake@input[["gene_set"]],
  what = "character"
)


GSEA_list <- list()

for (cl in names(DE_by_cluster)) {
  
  log_info(paste0("Running GSEA for cluster ", cl))

  res <- run_gsea_cluster(
    de_table = DE_by_cluster[[cl]],
    cluster_name = cl,
    gene_set = DE_genes
  )
  
  if (!is.null(res)) {
    GSEA_list[[cl]] <- res
  }
}

GSEA_table = data.frame(cluster = clusters)
GSEA_table = left_join(GSEA_table, a, by = "cluster")

GSEA_table = left_join(GSEA_table,bind_rows(GSEA_list))
GSEA_table = GSEA_table %>%
  mutate(
    passes_filter = ifelse(P.value < snakemake@params[["pvalue_cutoff"]] & FDR < snakemake@params[["FDR_cutoff"]], "Yes", "No")
  )


write.csv(GSEA_table, snakemake@output[["gsea_table"]], row.names = FALSE)