#!/bin/bash
#SBATCH --job-name=profile_config_test
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --cluster=gpsc8 
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=72:00:00

source path/to/source/conda/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.9.0
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile abs/path/to/profile/directory/code/metatranscriptomics-snakemake/profiles/slurm \
    --configfile abs/path/to/the/file/including/file/name/code/metatranscriptomics-snakemake/config/config_GPSC.yaml \
    --conda-prefix abs/path/to/where/you/created/the/conda/environment/code/metatranscriptomics-snakemake/metatranscriptomic-conda-env \
    --printshellcmds \
    --latency-wait 120 \
    --keep-going
    
