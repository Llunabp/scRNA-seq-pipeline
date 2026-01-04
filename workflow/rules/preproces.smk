rule load_counts:
    conda:
        "../envs/r_env.yaml"
    params:
        case_names = list(config["CASES"].keys()),
        control_names = list(config["CONTROLS"].keys()),
        min_cells = config["READ_DATA"]["MIN_CELLS_PER_GENE"],
        min_features = config["READ_DATA"]["MIN_FEATURES_PER_CELL"],
        scale_factor = config["READ_DATA"]["SCALE_FACTOR"],
        normalization_method = config["READ_DATA"]["NORMALIZATION_METHOD"],
        selection_method = config["READ_DATA"]["SELECTION_METHOD"],
        nfeatures = config["READ_DATA"]["NFEATURES"],
    input:
        cases_mtx = [config["CASES"][s]["mtx"] for s in config["CASES"]],
        cases_genes = [config["CASES"][s]["genes"] for s in config["CASES"]],
        cases_barcodes = [config["CASES"][s]["barcodes"] for s in config["CASES"]],
        controls_mtx = [config["CONTROLS"][s]["mtx"] for s in config["CONTROLS"]],
        controls_genes = [config["CONTROLS"][s]["genes"] for s in config["CONTROLS"]],
        controls_barcodes = [config["CONTROLS"][s]["barcodes"] for s in config["CONTROLS"]],

    output:
        rds = OUTDIR / "combined.rds",
        cell_counts = REPORT_DIR_TABLES / "cell_counts.csv"
    log:
        LOGDIR / "load_counts" / "log.txt"
    script:
        "../scripts/read_data.R"

rule reduce_data:
    conda:
        "../envs/r_env.yaml"
    params:
        npcs = config["PCA"]["N_PCS"],
        resolution = config["CLUSTERING"]["RESOLUTION"],
        algorithm = config["CLUSTERING"]["ALGORITHM"],
        n_neighbors = config["UMAP"]["N_NEIGHBORS"],
        min_dist = config["UMAP"]["MIN_DIST"],
        min_pct = config["MARKERS"]["MIN_PCT"],
        test_use = config["MARKERS"]["TEST"],
        logfc_threshold = config["MARKERS"]["LOGFC_THRESHOLD"],
        seed = config["SEED"],
    input:
        rds = OUTDIR / "combined.rds"
    output:
        reduction = OUTDIR / "combined_reduction.rds",
        markers = OUTDIR / "genes.csv",
        elbowplot = REPORT_DIR_PLOTS / "elbow_plot.png"
    log:
        LOGDIR / "reduce_data" / "log.txt"
    script:
        "../scripts/reduce_data.R"

rule annotated_genes:
    conda:
        "../envs/r_env.yaml"
    input:
        genes = OUTDIR / "genes.csv",
        markers = config["MARKERS_CELLS"]
    output:
        annotated_genes = REPORT_DIR_TABLES / "annotated_genes.csv"
    log:
        LOGDIR / "annotate_genes" / "log.txt"
    script:
        "../scripts/annotate_genes.R"