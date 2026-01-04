# Pipeline Instructions
To run the workflow, download the input data using [`00_download_data.sh`](/00_download_data.sh). It includes case and control studies count matrix with its own barcode and gene files. Also include enfisema gene set file and marker associated to cell identities file. 

Then, execute [`01_run_default.sh`](/01_run_default.sh). You can previous observe what is going to be executed using [`00_dry_run.sh`](/00_dry_run.sh). The default settings are already established. To modify them, check [`config`](/config/README.md) info.

Example of use:

```shell
snakemake --use-conda -c 1 
```
