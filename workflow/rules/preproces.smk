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
        rds = OUTDIR / "combined.rds"
    log:
        LOGDIR / "load_counts" / "log.txt"
    script:
        "../scripts/read_data.R"