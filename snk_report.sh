#!/bin/bash
#SBATCH --job-name=dry_run
#SBATCH --output=dry_run_%j.out
#SBATCH --error=dry_run_%j.err
#SBATCH --cluster=gpsc8
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=00:30:00

source /gpfs/fs7/aafc/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.6.0
export PATH="$PWD/bin:$PATH"

snakemake \
    --profile /gpfs/fs7/aafc/projects/J-003165_abcc_rcba/code/metatranscriptomics-snakemake/profiles/slurm \
    --configfile /gpfs/fs7/aafc/projects/J-003165_abcc_rcba/code/metatranscriptomics-snakemake/config/config.yaml \
    --conda-prefix /gpfs/fs7/aafc/scratch/kjg000/tmpdir/metatranscriptomics_conda \
    --printshellcmds \
    --keep-going \
    --report /gpfs/fs7/aafc/projects/J-003165_abcc_rcba/code/metatranscriptomics-snakemake/metatranscriptomics_report.html 
