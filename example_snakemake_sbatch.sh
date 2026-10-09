#!/bin/bash
#SBATCH --job-name=run_snakemake
#SBATCH --output=run_snakemake_%j.out
#SBATCH --error=run_snakemake_%j.err
#SBATCH --clusters=<CLUSTER_NAME>
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=<MEMORY_MB>
#SBATCH --time=<HH:MM:SS>

set -euo pipefail

# Activate the environment containing Snakemake and the SLURM executor plugin.
source "/absolute/path/to/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake-9.20.0

# Run from the repository root.
PIPELINE_DIR="/absolute/path/to/code/metatranscriptomics-snakemake"
cd "$PIPELINE_DIR"

export PATH="$PIPELINE_DIR/bin:$PATH"

# This directory must be writable on the compute nodes.
export TMPDIR="/absolute/path/to/scratch/${USER}/tmpdir"
mkdir -p "$TMPDIR"

snakemake \
    --snakefile "$PIPELINE_DIR/workflow/Snakefile" \
    --profile "$PIPELINE_DIR/profiles/slurm" \
    --configfile "$PIPELINE_DIR/config/config.yaml" \
    --conda-prefix "/absolute/path/to/shared/conda/metatranscriptomic-conda-env" \
    --printshellcmds \
    --latency-wait 120 \
    --keep-going
