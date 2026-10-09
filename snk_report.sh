#!/bin/bash
#SBATCH --job-name=snakemake_report
#SBATCH --output=snakemake_report_%j.out
#SBATCH --error=snakemake_report_%j.err
#SBATCH --clusters=<CLUSTER_NAME>
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=<MEMORY_MB>
#SBATCH --time=<HH:MM:SS>

set -euo pipefail

# Activate the same Snakemake environment used for the analysis.
source "/absolute/path/to/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake-9.20.0

PIPELINE_DIR="/absolute/path/to/code/metatranscriptomics-snakemake"
cd "$PIPELINE_DIR"

export TMPDIR="/absolute/path/to/scratch/${USER}/tmpdir"
mkdir -p "$TMPDIR"

# Generate the report after the workflow has completed.
REPORT_DIR="$PIPELINE_DIR/results"
mkdir -p "$REPORT_DIR"

snakemake \
    --snakefile "$PIPELINE_DIR/workflow/Snakefile" \
    --configfile "$PIPELINE_DIR/config/config.yaml" \
    --profile "$PIPELINE_DIR/profiles/slurm" \
    --conda-prefix "/absolute/path/to/shared/conda/metatranscriptomic-conda-env" \
    --report "$REPORT_DIR/metatranscriptomics_report.html"
