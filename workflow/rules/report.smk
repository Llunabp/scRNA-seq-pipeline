rule report:
    conda:
        "../envs/quarto.yaml"
    params:
        legends = LEGENDS,
        reductions = REDUCTIONS,
        dim_plot_path = FILE_DIM,
        case = config["GEO_ACCESION"]["Cases"],
        control = config["GEO_ACCESION"]["Controls"],
        gsea_pval = config["GSEA"]["PVALUE"],
        gsea_mingenes = config["GSEA"]["MIN_GENES"]
        de_fdr = config["DE"]["FDR"],
        de_fc = config["DE"]["FC"],
        sc_min_cells = config["READ_DATA"]["MIN_CELLS_PER_GENE"],
        sc_min_features = config["READ_DATA"]["MIN_FEATURES_PER_CELL"],
        sc_scale_factor = config["READ_DATA"]["SCALE_FACTOR"],
        sc_normalization_method = config["READ_DATA"]["NORMALIZATION_METHOD"],
        sc_selection_method = config["READ_DATA"]["SELECTION_METHOD"],
        sc_nfeatures = config["READ_DATA"]["NFEATURES"],

    input:
        qmd = "template.qmd",
        dim_plots = expand(REPORT_DIR_PLOTS / "{reduction}" / "dim_plot_{legend}.png", reduction=config["DIM_REPRESENTATION"]["REDUCTION_METHOD"], legend=config["DIM_REPRESENTATION"]["LEGEND"]),
        heatmap = REPORT_DIR_PLOTS / "heatmap.png",
        feature_plot = REPORT_DIR_PLOTS / "feature_plot.png",
        de_summary = REPORT_DIR_TABLES / "DE_summary.csv",
        gsea_table = REPORT_DIR_TABLES / "GSEA_results.csv",
        annotated_genes = REPORT_DIR_TABLES / "annotated_genes.csv",
    output:
        OUTDIR / "report" / "report.html"
    log:
        LOGDIR / "report" / "log.txt"
    shell:
        "set +o pipefail; "
            "Rscript -e 'library(quarto)' "
            "-e \"quarto_render("
                "input = '{input.qmd}', "
                "execute_params=list("
                    "heatmap='{input.heatmap}', "
                    "feature_plot='{input.feature_plot}', "
                    "de_summary='{input.de_summary}', "
                    "gsea_table='{input.gsea_table}', "
                    "annotated_genes='{input.annotated_genes}', "
                    "legends='{params.legends}', "
                    "reductions='{params.reductions}', "
                    "dim_plot_path='{params.dim_plot_path}', "
                    "gsea_pval='{params.gsea_pval}', "
                    "gsea_mingenes='{params.gsea_mingenes}', "
                    "de_fdr='{params.de_fdr}', "
                    "de_fc='{params.de_fc}', "
                    "sc_min_cells='{params.sc_min_cells}', "
                    "sc_min_features='{params.sc_min_features}', "
                    "sc_scale_factor='{params.sc_scale_factor}', "
                    "sc_normalization_method='{params.sc_normalization_method}', "
                    "sc_selection_method='{params.sc_selection_method}', "
                    "sc_nfeatures='{params.sc_nfeatures}'))\" "
            ">{log} 2>&1 && "
            'mv "$(dirname {input.qmd:q})/report.html" {output.html:q}'