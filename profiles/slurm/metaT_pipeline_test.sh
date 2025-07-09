#!/bin/bash
#SBATCH --job-name=metaT_full_test
#SBATCH --output=metaT_full_test_%j.out 
#SBATCH --error=metaT_full_test_%j.err 
#SBATCH --cluster=gpsc8 
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=8:00:00

source /gpfs/fs7/aafc/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.6.0
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile /home/kjg000/J-003165_abcc_rcba/code/metatranscriptomics-pipeline-snakemake/profiles/slurm \
    --configfile /gpfs/fs7/aafc/projects/J-003165_abcc_rcba/code/metatranscriptomics-pipeline-snakemake/config/config.yaml \
    --printshellcmds 
