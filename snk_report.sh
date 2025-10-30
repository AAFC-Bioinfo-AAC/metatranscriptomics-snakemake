#!/bin/bash
#SBATCH --job-name=full_test_run_report
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --cluster=gpsc8
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=00:30:00

source path/to/source/conda/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.9.0
export PATH="$PWD/bin:$PATH"

#Can just run this from the head node.
snakemake --report path/to/where/you/whant/the/report/metatranscriptomics_report.html
