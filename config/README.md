# Instructions

To run the pipeline, you will need an environment with `snakemake`
(check [the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)).

## Inputs and outputs

The workflow execution parameters are set in three configuration files in YAML format. Workflow configuration are set via [`config.yaml`](/config/config.yaml) and study configuration are set via [`target.yaml`](/config/target.yaml)

## Workflow configuration variables

All of the following variables are pre-defined in the [`config.yaml`](/config/config.yaml) file:
- `READ_DATA`: preprocessing and normalization configuration:
  - `MIN_CELLS_PER_GENE` (Default: **3**): minimum number of cells in which a gene must be detected to be retained.
  - `MIN_FEATURES_PER_CELL` (Default: **200**): minimum number of detected genes required for a cell to be retained.
  - `SCALE_FACTOR` (Default: **10000**): scaling factor used during normalization to correct for sequencing depth.
  - `NORMALIZATION_METHOD` (Default: **LogNormalize**): method used to normalize gene expression values.
  - `SELECTION_METHOD` (Default: **vst**): method used to identify highly variable genes. Options: vst, mean.var.plt, dispersion
  - `NFEATURES` (Default: **2000**): number of highly variable genes selected for downstream analyses.
- `PCA`: principal component analysis configuration:
  - `N_PCS` (Default: **50**): number of principal components explored during PCA.
- `CLUSTERING`: graph-based clustering configuration:
  - `RESOLUTION` (Default: **0.4**): resolution parameter controlling cluster granularity.
  - `ALGORITHM` (Default: **3**): clustering algorithm used. Options: 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm.
- `UMAP`: UMAP dimensionality reduction configuration:
  - `N_NEIGHBORS` (Default: **30**): number of nearest neighbors used for UMAP graph construction.
  - `MIN_DIST` (Default: **0.3**): minimum distance between points in the UMAP embedding.
- `MARKERS`: marker gene detection parameters:
  - `LOGFC_THRESHOLD` (Default: **0.25**): minimum log2 fold-change required for marker detection.
  - `MIN_PCT` (Default: **0.25**): minimum fraction of cells expressing a gene within a cluster.
  - `TEST` (Default: **wilcox**): statistical test used to identify marker genes. Options: 
- `SEED` (Default: **123**): random seed used to ensure reproducibility.
- `DIM_REPRESENTATION`: dimensionality reduction visualization settings:
  - `REDUCTION_METHOD` (Default: **tsne, umap**): dimensionality reduction methods used for visualization. Options: pca, tsne, umap
  - `LEGEND` (Default: **long_name, group, orig.ident**): metadata variables used to color cells in embeddings. Options: long_name, group, orig.ident, seurat_cluster
- `HEATMAP`: heatmap visualization parameters:
  - `N_MARKERS` (Default: **5**): number of top marker genes per cluster displayed in heatmaps.
- `HGV`: highly variable gene visualization parameters:
  - `N_MARKERS` (Default: **10**): number of highly variable genes highlighted in variability plots
- `GSEA`: Gene Set Enrichment Analysis configuration:
  - `MIN_GENES` (Default: **10**): minimum number of genes required to perform GSEA for a cluster.
  - `PVALUE` (Default: **0.05**): nominal p-value threshold for significance.
  - `FDR` (Default: **0.1**): false discovery rate threshold for multiple testing correction.
  - `P_ADJ_METHOD` (Default: **BH**): method used for p-value adjustment.
- `DE`: differential expression configuration:
  - `FDR` (Default: **0.05**): false discovery rate threshold to define significantly differentially expressed genes.
  - `FC` (Default: **0.4**): minimum absolute log2 fold-change required to retain differentially expressed genes.

All of the following variables are pre-defined in the [`target.yaml`](/config/target.yaml) file:
- `GEO_ACCESION`: GEO datasets used in the analysis:
  - `Cases`: GEO accession containing case samples.
  - `Controls`: GEO accession containing control samples.

- `CASES`: case sample input files:
  - `Sample`:
    - `mtx`: path to the count matrix file.
    - `genes`: path to the gene annotation file.
    - `barcodes`: path to the cell barcode file.
  
- `CONTROLS`: control sample input files:
  - `Sample`:
    - `mtx`: path to the count matrix file.
    - `genes`: path to the gene annotation file.
    - `barcodes`: path to the cell barcode file.
- `OUTPUT_DIRECTORY`: directory where all pipeline results and intermediate files are stored.
- `GENE_SET_FILE`: file containing the predefined gene set used for Gene Set Enrichment Analysis (GSEA).
- `MARKERS_CELLS`: file containing reference marker genes used for cell type annotation.


An example "config" file could look like this:

```yaml
GEO_ACCESION:
  Cases: "GSEXXXXXX" 
  Controls: "GSEXXXXXX" 
  
CASES:
  SAMPLE_1:
    mtx: "path/to/case_sample_1/matrix.mtx"
    genes: "path/to/case_sample_1/genes.tsv"
    barcodes: "path/to/case_sample_1/barcodes.tsv"
  SAMPLE_2:
    mtx: "path/to/case_sample_2/matrix.mtx"
    genes: "path/to/case_sample_2/genes.tsv"
    barcodes: "path/to/case_sample_2/barcodes.tsv"
  ...

CONTROLS:
  SAMPLE_1:
    mtx: "path/to/control_sample_1/matrix.mtx"
    genes: "path/to/control_sample_1/genes.tsv"
    barcodes: "path/to/control_sample_1/barcodes.tsv"
  SAMPLE_2:
    mtx: "path/to/control_sample_2/matrix.mtx"
    genes: "path/to/control_sample_2/genes.tsv"
    barcodes: "path/to/control_sample_2/barcodes.tsv"
  ...

OUTPUT_DIRECTORY: "path/to/output_dir"
GENE_SET_FILE: "path/to/gene_set.txt"
MARKERS_CELLS: "path/to/markers_cells.csv"
```
Configuration variables may also be set through the `--config` option.

## Run modes

To run the analysis with the default configuration, run the following command
(change the `-c/--cores` argument to use a different number of CPUs (n) and study for the study name that you are using):

```shell
snakemake --use-conda -c n 
```

## Workflow graphs

To generate a simplified rule graph, run:

```shell
snakemake --rulegraph | dot -Tpng > .rulegraph.png
```
![Snakemake rule graph from workflow.](/.rulegraph.png)


To generate the directed acyclic graph (DAG) of all rules
to be executed, run:

```shell
snakemake --forceall --dag | dot -Tpng > .dag.png
```
![Snakemake DAG from workflow.](/.dag.png)






