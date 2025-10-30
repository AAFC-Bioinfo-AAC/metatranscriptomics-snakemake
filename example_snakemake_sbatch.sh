#!/bin/bash
#SBATCH --job-name=run_snakemake.sh
#SBATCH --output=run_snakemake_%j.out 
#SBATCH --error=run_snakemake_%j.err 
#SBATCH --cluster=<CLUSTER_NAME>
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --mem=<MEMORY_MB>         # e.g., 2000
#SBATCH --time=<HH:MM:SS>         # Must be long enough for completion of workflow 

source abs/path/to/source/conda/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.9.0
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile abs/path/to/profile/directory/code/metatranscriptomics-snakemake/profiles/slurm \
    --configfile abs/path/to/the/file/including/file/name/code/metatranscriptomics-snakemake/config/config_GPSC.yaml \
    --conda-prefix abs/path/to/where/you/created/the/conda/environment/code/metatranscriptomics-snakemake/metatranscriptomic-conda-env \
    --printshellcmds \
    --latency-wait 120 \
    --keep-going
    
