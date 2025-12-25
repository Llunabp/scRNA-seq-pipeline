rule dif_expression:
    conda:
        "../envs/r_env.yaml"
    params:
        min_genes = config["GSEA"]["MIN_GENES"],
        pvalue_cutoff = config["GSEA"]["PVALUE"],
        logfc_threshold = config["MARKERS"]["LOGFC_THRESHOLD"],
        min_pct = config["MARKERS"]["MIN_PCT"],
        FDR = config["DE"]["FDR"],
        FC = config["DE"]["FC"],
        seed = config["SEED"]
    input:
        rds = OUTDIR / "combined_reduction.rds",
        gene_set = config["GENE_SET_FILE"]
    output:
        de_summary = REPORT_DIR_TABLES / "DE_summary.csv",
        gsea_table = REPORT_DIR_TABLES / "GSEA_results.csv"
    log:
        LOGDIR / "dif_expression" / "log.txt"
    script:
        "../scripts/DE_GSEA.R"