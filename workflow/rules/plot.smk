rule dim_plot:
    conda:
        "../envs/r_env.yaml"
    params:
        reduction = lambda wildcards: wildcards.reduction,
        legend = lambda wildcards: wildcards.legend,
    input:
        rds = OUTDIR / "combined_reduction.rds"
    output:
        dimplot = REPORT_DIR_PLOTS / "{reduction}" / "dim_plot_{legend}.png"
    log:
        LOGDIR / "dim_plot_{reduction}" / "log_{legend}.txt"
    script:
        "../scripts/plot_dim.R"

rule plot_heatmap:
    conda:
        "../envs/r_env.yaml"
    params:
        n_markers = config["HEATMAP"]["N_MARKERS"],
        logfc_threshold = config["MARKERS"]["LOGFC_THRESHOLD"], 
        min_pct = config["MARKERS"]["MIN_PCT"],
        seed = config["SEED"]
    input:
        rds = OUTDIR / "combined_reduction.rds"
    output:
        heatmap = REPORT_DIR_PLOTS / "heatmap.png"
    log:
        LOGDIR / "plot_heatmap" / "log.txt"
    script:
        "../scripts/plot_heatmap.R"

rule plot_features_dim:
    conda:
        "../envs/r_env.yaml"
    input:
        rds = OUTDIR / "combined_reduction.rds",
        genes = OUTDIR / "genes.csv"
    output:
        featureplot = REPORT_DIR_PLOTS / "feature_plot.png"
    log:
        LOGDIR / "plot_features_dim" / "log.txt"
    script:
        "../scripts/plot_features_dim.R"

rule plot_HGVc:
    conda:
        "../envs/r_env.yaml"
    params:
        top_features = config["HGV"]["N_MARKERS"]
    input:
        rds = OUTDIR / "combined.rds",
    output:
        plot = REPORT_DIR_PLOTS / "HGVs_plot.png",
    log:
        LOGDIR / "plot_HGVs" / "log.txt"
    script:
        "../scripts/plot_HGVs.R"

rule plot_violin:
    conda:
        "../envs/r_env.yaml"
    input:
        rds = OUTDIR / "combined.rds",
    output:
        plot_feature = REPORT_DIR_PLOTS / "violin_features.png",
        plot_count = REPORT_DIR_PLOTS / "violin_counts.png",
    log:
        LOGDIR / "plot_violin" / "log.txt"
    script:
        "../scripts/plot_violin.R"