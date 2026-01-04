#!/usr/bin/env bash
#SBATCH --job-name=bioinsilico
#SBATCH --output=slurm-%x-%j.out
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=4
#SBATCH --qos=medium

set -e

logthis () {
    echo $(date) "|" $@
}

NCPU=4


logthis "Running Pipeline"
snakemake --use-conda -c ${NCPU} 

logthis "Done!"
