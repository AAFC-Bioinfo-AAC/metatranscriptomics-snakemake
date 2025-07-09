# Snakemake Profile for running on SLURM

## Directory layout 
my_project/
├── Snakefile
├── config/
│   └── config.yaml          # ← workflow variables for your code/rules
└── profiles/
    └── slurm/
        ├── config.yaml      # ← profile settings for snakemake/SLURM and per-rule resource config 
      

## Main profile
- Controls how Snakemake runs overall (max cores, jobs, conda, default-resources, latency-wait, etc)
- Per rule resources need to go here. Snakemake no longer take a cluster_config.ymal file
- Is used when you launch Snakemake with --profile <dir>.
### Example config.yaml
```bash
### How Snakemake will run on SBATCH ###
cores: 64 # total number of cores Snakemake can request at any time
jobs: 20 # max amount of jobs Snakemake runs at once
latency-wait: 60 # gives time for containers to start, and for other I/O delays
rerun-incomplete: true
quiet: false # Makes Snakemake output more verbose about its operations
retries: 2               # so jobs killed by IO hiccups are auto retried
max-jobs-per-second: 10  # slowing down submission rate (tune for your site!)

### SLURM ###
slurm: true
default-resources:
  - mem_mb=8000             # default memory per job: 8GB (change as desired)
  - slurm_partition=standard  # SLURM partition to use by default
  - slurm_account=aafc_aac # SLURM account

### Env modules ###
# use-envmodules: false # use-envmodules: true only if: Your system disables Conda/containers and expects you to use module load bioinfo-tool for each step or you have a properly configured profile and know which modules are needed for every rule.

### Conda ###
use-conda: true
conda-frontend: mamba   # if you have mamba installed

### Resource scopes ###
# Affects how Snakemake defines resource limits
# Note that cores and threads are always considered local 
# set-resource-scopes:
    cores: local 

### Containerization if using containers like Docker ###
# use-singularity: true 
# singularity-prefix-dir: 
# singularity-args:
# cleanup-containers: false 

### Rule-specific configs ###
# These go into the cluster_config.yaml
## Per rule resources
set-resources:
  fastp_pe:
    cpus: 4
    mem_mb: 40000
    time: "00:40:00"
    partition: standard
    account: aafc_aac

  bowtie2_align:
    cpus: 24
    mem_mb: 48000
    time: "00:30:00"
    partition: standard
    account: aafc_aac

  extract_unmapped_fastq:
    cpus: 60
    mem_mb: 64000
    time: "00:30:00"
    partition: standard
    account: aafc_aac

  sortmerna_pe:
    cpus: 48
    mem_mb: 32000
    time: "06:00:00"
    partition: standard
    account: aafc_aac

  kraken2:
    cpus: 2
    mem_mb: 600000
    time: "00:30:00"
    partition: large
    account: aafc_aac__large

  bracken:
    cpus: 2
    mem_mb: 4000
    time: "00:10:00"
    partition: standard
    account: aafc_aac

  combine_bracken_outputs:
    cpus: 1
    mem_mb: 2000
    time: "00:20:00"
    partition: standard
    account: aafc_aac

  bracken_extract:
    cpus: 1
    mem_mb: 2000
    time: "00:10:00"
    partition: standard
    account: aafc_aac

  rgi_reload_database:
    cpus: 1
    mem_mb: 2000
    time: "00:30:00"
    partition: standard
    account: aafc_aac

  rgi_bwt:
    cpus: 20
    mem_mb: 64000
    time: "01:00:00"
    partition: standard
    account: aafc_aac

  rna_spades:             
    cpus: 48
    mem_mb: 64000
    time: "04:00:00"
    partition: standard
    account: aafc_aac   

  rnaquast_busco:
    cpus: 4
    mem_mb: 16000
    time: "00:10:00"
    partition: standard
    account: aafc_aac  

  megahit_coassembly:
    cpus: 60
    mem_mb: 512000
    time: "04:00:00"
    partition: standard
    account: aafc_aac 

  index_coassembly:
    cpus: 8
    mem_mb: 16000
    time: "02:00:00"
    partition: standard
    account: aafc_aac

  bowtie2_map_transcripts:
    cpus: 16
    mem_mb: 32000
    time: "00:30:00"
    partition: standard
    account: aafc_aac

  assembly_stats_depth:
    cpus: 2
    mem_mb: 2000
    time: "00:30:00"
    partition: standard
    account: aafc_aac

  prodigal_genes:
    cpus: 1
    mem_mb: 2000
    time: "01:00:00"
    partition: standard
    account: aafc_aac

  featurecounts:
    cpus: 4
    mem_mb: 8000
    time: "00:10:00"
    partition: standard
    account: aafc_aac
  ```

### Example
```bash
#!/bin/bash
#SBATCH --job-name=metaT_full_test
#SBATCH --output=metaT_full_test_%j.out 
#SBATCH --cluster=gpsc8 
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=1:00:00

source ~/.bashrc
mamba activate snakemake-9.6.0
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile path/to/profiles/slurm/ \
    --configfile my_project/config/cluster_config.yaml \
    --printshellcmds \
    --forceall ## keep for testing. Remove for actual runs
```     
