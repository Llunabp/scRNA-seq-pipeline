rule write_config_json:
    params:
        content=config,
    output:
        OUTDIR / "config.json",
    priority: 50
    run:
        import json

        with open(output[0], "w") as fw:
            json.dump(params.content, fw, indent=2)

rule check_legend:
    conda:
        "../envs/r_env.yaml"
    params:
        possibilities = ["long_name", "seurat_clusters", "group", "orig.ident"],
        legends = LEGENDS
    output:
        touch(OUTDIR / "checks_legend.done")
    log:
        LOGDIR / "check_legend" / "log.txt"
    run:
        legends_list = params.legends.split(";")
        if not all(legend in params.possibilities for legend in legends_list):
            raise(f"Error: One or more legends in {params.legends} are not among the allowed options: {', '.join(params.possibilities)}")

rule check_reduction:
    conda:
        "../envs/r_env.yaml"
    params:
        possibilities = ["umap", "tsne", "pca"],
        reductions = REDUCTIONS
    output:
        touch(OUTDIR / "checks_reduction.done")
    log:
        LOGDIR / "check_reduction" / "log.txt"
    run:
        reductions_list = params.reductions.split(";")
        if not all(reduction in params.possibilities for reduction in reductions_list):
            raise(f"Error: One or more reductions in {params.reductions} are not among the allowed options: {', '.join(params.possibilities)}")