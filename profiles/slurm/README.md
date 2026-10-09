# Snakemake Profile for Running on SLURM

This profile defines how Snakemake submits jobs to SLURM and assigns job resources.

## Directory layout

Paths below are relative to the repository root, `metatranscriptomics-snakemake/`.

| Path | Contents |
|---|---|
| `workflow/Snakefile` | Main workflow definition. |
| `workflow/rules/` | Snakemake rule modules. |
| `workflow/scripts/` | Scripts called by workflow rules. |
| `workflow/envs/` | Rule-specific Conda environment definitions. |
| `resources/` | Reference files or database setup information. |
| `profiles/slurm/config.yaml` | SLURM execution profile. |
| `config/config.yaml` | Workflow paths, parameters and thread settings. |
| `config/samplesheet.csv` | Sample IDs and paired fastq filenames. |
| `run_snakemake.sh` | SLURM launcher script. |
| `.env` | Project paths and environment variables. |

Tool log locations are configured through `log_files` in `config/config.yaml`.

## Main profile

Save the profile as `profiles/slurm/config.yaml`.

An editable template is provided as `profiles/slurm/example_config.yaml`. Copy it to `config.yaml` and update the settings for your cluster.

```bash
cp profiles/slurm/example_config.yaml profiles/slurm/config.yaml
```

Invoke Snakemake from the repository root using:

```bash
snakemake \
    --profile /absolute/path/to/metatranscriptomics-snakemake/profiles/slurm
```

The profile controls execution settings. Sample information, database paths and analysis parameters belong in `config/config.yaml`.

## Temporary working directory

The profile's `envvars` setting lists environment variables to pass to jobs; it does not assign their values.

Define a writable temporary directory in the launcher before invoking Snakemake:

```bash
export TMPDIR="/absolute/path/to/scratch/${USER}/tmpdir"
mkdir -p "$TMPDIR"
```

This path must be writable on the compute nodes. The revised rules create private temporary directories beneath it and remove them when they finish.

Final outputs and intermediate files passed between jobs must remain on storage accessible to the relevant compute nodes.

## Example profile

The following example is intended for Snakemake 9.20.0 with the SLURM executor plugin installed in the Snakemake environment.

Replace angle-bracket placeholders with valid settings for your cluster. Remove the `clusters` or `qos` entries if those options are not required.

Runtime values are in minutes. `mem_mb` specifies total memory per job.

The resource allocations are examples and should be adjusted using observed memory consumption and execution times. Assign a suitable partition to any rule whose requirements exceed the standard partition's limits.

```yaml
# SLURM execution profile for Snakemake 9.20.0.
# Replace angle-bracket placeholders with valid cluster settings.

executor: slurm

# Maximum thread allocation available to an individual remote job.
cores: 60

# Maximum number of active cluster jobs.
jobs: 10

# CPU limit for rules executed locally by the controller.
local-cores: 1

latency-wait: 60
rerun-incomplete: true
retries: 2
max-jobs-per-second: 2

# Reevaluate outputs when inputs, parameters, code or environments change.
rerun-triggers: [mtime, input, params, software-env, code]

# Values must be defined before Snakemake starts.
envvars:
  - TMPDIR

use-conda: true

# Cluster defaults.
default-resources:
  slurm_account: "<STANDARD_ACCOUNT_NAME>"
  slurm_partition: "<STANDARD_PARTITION_NAME>"
  clusters: "<CLUSTER_NAME>"
  qos: "<QOS_LEVEL>"
  runtime: 60
  mem_mb: 4000

# Example starting allocations.
set-resources:
  fastp_pe:
    mem_mb: 4000
    runtime: 40

  bowtie2_align:
    mem_mb: 24000
    runtime: 90

  extract_unmapped_fastq:
    mem_mb: 18000
    runtime: 60

  sortmerna_pe:
    mem_mb: 32000
    runtime: 360

  kraken2:
    slurm_account: "<LARGE_ACCOUNT_NAME>"
    slurm_partition: "<LARGE_PARTITION_NAME>"
    mem_mb: 840000
    runtime: 120

  bracken:
    mem_mb: 4000
    runtime: 10

  combine_bracken_outputs:
    mem_mb: 2000
    runtime: 20

  bracken_extract:
    mem_mb: 2000
    runtime: 10

  rgi_validate_database:
    mem_mb: 2000
    runtime: 30

  rgi_bwt:
    mem_mb: 64000
    runtime: 60

  rna_spades:
    mem_mb: 64000
    runtime: 240

  rnaquast_busco:
    mem_mb: 16000
    runtime: 120

  megahit_coassembly:
    mem_mb: 256000
    runtime: 720

  index_coassembly:
    mem_mb: 16000
    runtime: 120

  bowtie2_map_transcripts:
    mem_mb: 32000
    runtime: 720

  assembly_stats_depth:
    mem_mb: 2000
    runtime: 30

  prodigal_genes:
    mem_mb: 2000
    runtime: 60

  featurecounts:
    mem_mb: 8000
    runtime: 10

  prepare_cazyme_proteins:
    mem_mb: 2000
    runtime: 30

  cazyme_annotation:
    mem_mb: 16000
    runtime: 240

  cazyme_rna_counts:
    mem_mb: 8000
    runtime: 30
```

## Resource settings

- Rule thread settings determine the CPU request for each job.
- `jobs: 10` limits concurrent cluster jobs.
- `cores: 60` does not impose an aggregate limit of 60 CPUs across all remote jobs.
- `local-cores: 1` limits CPU use by local rules.
- The launcher's `#SBATCH` directives allocate resources to the Snakemake controller. Rule jobs receive their own allocations through the executor.

With `rna_spades.memory_gb: 60`, retain sufficient memory in the `rna_spades` allocation. The revised rule checks the SPAdes memory setting against the job's allocated `mem_mb`.

The `default-resources` settings apply unless overridden by a rule or a rule-specific profile setting. The Kraken2 entry above overrides the standard account and partition.

CAZyme analysis is part of the default workflow. Its three rule allocations are included above.

## Conda environments

The profile enables rule-specific Conda environments with `use-conda: true`.

Supply the environment cache path through `--conda-prefix` in the launcher. Use the same prefix when creating environments and running the workflow.

For example:

```bash
snakemake \
    --profile /absolute/path/to/metatranscriptomics-snakemake/profiles/slurm \
    --conda-prefix /absolute/path/to/metatranscriptomic-conda-env
```

## References

- [Snakemake 9.20.0 command-line documentation](https://snakemake.readthedocs.io/en/v9.20.0/executing/cli.html)
- [SLURM executor plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
