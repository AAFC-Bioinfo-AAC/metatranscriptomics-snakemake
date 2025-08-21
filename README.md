# Metatranscriptomics Snakemake Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## About

The Metatranscriptomics Snakemake Pipeline uses paired-end FASTQ files from Illumina shotgun metatranscriptomic sequencing as input. The pipeline can be broken down into four main stages: sample read processing, sort read analysis, individual sample assembly, and co-assembly. Sample processing consists of fastp, Bowtie2, and SortMeRNA to perform quality filtering, remove host and PhiX contamination, and ribosomal (r)RNA depletion. These cleaned reads are used for the sort read analysis consisting of taxonomic classification with Kraken2 using GTDB and antimicrobial gene profiling with RGI using CARD. Individual samples are assembled into presumptive messenger (m)RNA transcripts using RNA SPAdes. Assembly quality is evaluated with rnaQUAST. The co-assembly stage prepares the data for gene expression analysis. All cleaned reads are co-assembled with MEGAHIT, and the resulting co-assembly is indexed with Bowtie2. The cleaned sample reads are then mapped back to the co-assembly, and SAMtools is used generate assembly statistics, mapping summaries, and sequencing depth across the co-assembly. With Prodigal the protein and nucleotide coding regions of the co-assembly are predicted. FeatureCounts quantifies the predicted coding regions and generates a table for gene expression analysis. If metagenomic sequencing was done for these samples then the trimmed and host/PhiX removed metagenomic reads should be used for the co-assembly step.

### Future Additions

Contemplated future additions to this pipeline include the tool CoverM to map the metatransriptomic reads to the assembled metagenomes, and a CAZyme analysis module.

---

## Table of Contents

- [Metatranscriptomics Snakemake Pipeline](#metatranscriptomics-snakemake-pipeline)
  - [About](#about)
    - [Future Additions](#future-additions)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
    - [Workflow diagram](#workflow-diagram)
    - [Snakemake rules](#snakemake-rules)
    - [Module  `preprocessing.smk` contains these rules:](#module--preprocessingsmk-contains-these-rules)
    - [Module  `sortmerna.smk` contains this rule:](#module--sortmernasmk-contains-this-rule)
    - [Module  `taxonomy.smk` contains these rules:](#module--taxonomysmk-contains-these-rules)
    - [Module  `amr_short_reads.smk` contains these rules:](#module--amr_short_readssmk-contains-these-rules)
    - [Module  `sample_assembly.smk` contains these rules:](#module--sample_assemblysmk-contains-these-rules)
    - [Module  `coassebly_annotation.smk` contains these rules:](#module--coassebly_annotationsmk-contains-these-rules)
  - [Data](#data)
  - [Parameters](#parameters)
  - [Usage](#usage)
    - [Pre-requisites](#pre-requisites)
      - [Software](#software)
      - [Databases](#databases)
    - [Setup Instructions](#setup-instructions)
      - [1. Installation](#1-installation)
      - [2. SLURM Profile](#2-slurm-profile)
        - [2.1. SLURM Profile Directory Structure](#21-slurm-profile-directory-structure)
        - [2.2. Profile Configuration](#22-profile-configuration)
      - [3. Configuration](#3-configuration)
        - [3.1. config/config.yaml](#31-configconfigyaml)
        - [3.2. Environment file](#32-environment-file)
        - [3.3. Sample list](#33-sample-list)
      - [4. Running the pipeline](#4-running-the-pipeline)
        - [4.1. Conda environments](#41-conda-environments)
        - [4.2. SLURM launcher](#42-slurm-launcher)
    - [Notes](#notes)
      - [Warnings](#warnings)
      - [Current issues](#current-issues)
      - [Resource usage](#resource-usage)
  - [OUTPUT](#output)
  - [Credits](#credits)
  - [Contribution](#contribution)
  - [License](#license)
  - [References](#references)
    - [Publications](#publications)
    - [Resources](#resources)
    - [Tools/Software](#toolssoftware)
  - [Citation](#citation)

---

## Overview

### Workflow diagram

 ```mermaid
 flowchart TD
    %% Title
    %% Metatranscriptome Assembly and Analysis Pipeline

    subgraph PREPROC [Pre-processing]
        A[Paired Reads] -->|QC & trim| B{fastp}
        B --> C[Trimmed Reads - temp]
        B --> L((Fastp QC Report))
        C -->|Host/PhiX removal| D{Bowtie2}
        D --> E[Non-host, Non-PhiX Reads]
    end

    subgraph DEPLETION [rRNA Depletion]
        E --> F{SortMeRNA}
        F --> G[rRNA-depleted Reads]
    end

    subgraph COASSEMBLY [Co-Assembly]
        %% megahit_coassembly rule
        G --> MH{MEGAHIT Co-assembly}
        MH --> FA((Co-assembled Transcripts))

        %% index_coassembly rule
        FA --> IB{Index Bowtie2}
        IB --> IB1((Coassembly index))

        %% bowtie2_map_transcripts rule (per-sample mapping)
        IB1 --> BT2{Bowtie2 Map Transcripts}
        G --> BT2
        BT2 --> BAM((Sample BAMs))

        %% assembly_stats_depth rule (per-sample)
        BAM --> STATS{Assembly Stats+Depth}
        STATS --> STATSOUT((Stats/Depth/IdxStats))

        %% prodigal_genes rule
        FA --> PG{Prodigal Genes}
        PG --> PROD_OUTS1((Predicted protein sequences
        Predicted nucleotide sequences))
        PG --> PROD_OUTS2((Gene annotation file
                        Simplified Annotation Format))

        %% featurecounts rule (per-sample)
        PROD_OUTS2 --> FC{featureCounts}
        BAM --> FC
        FC --> FCT((Sample Counts.txt))
    end

    subgraph ASSEMBLY [Sample Assembly]
        G --> H{rnaSPAdes}
        H --> I((Sample Transcripts))
    end

    subgraph QC_REPORTS [Reports and Files for Downstream Analysis]
        I --> S{rnaQUAST}
        S --> U((Assembly QC Report))
        G --> M{Kraken2}
        M --> N{Bracken}
        N --> O((Taxonomic Profile))
        G --> W{RGI}
        W --> Q((AMR Profile))
    end

    %% TEMP FILE STYLING
    style C fill:#f2f2f2,stroke-dasharray: 5 5
    style L fill:#f2f2f2,stroke-dasharray: 5 5
```

### Snakemake rules

The pipeline is modularized, with each module located in the `metatranscriptomics-snakemake/workflow/rules` directory. The modules are `preprocessing.smk`, `sortmerna.smk`, `taxonomy.smk`,`amr_short_reads.smk`, and `coassebly_annotation.smk`. 

---   
### Module  `preprocessing.smk` contains these rules:
🔹 **`rule fastp_pe` *Quality Control & Trimming***

- **Purpose:** Performs adapter trimming, quality trimming, and filtering of paired-end reads.
- **Inputs:** `samplesheet.csv` defines sample IDs and corresponding read pairs.
- **Outputs:**
  - Trimmed paired reads: `sample_r1.fastq.gz`, `sample_r2.fastq.gz`
  
- **Notes:**
  - Parameters are defined in **`config/config.ymal`** for `fastp`.
  - These files are marked as temporary in the rule: `sample_u1.fastq.gz`, `sample_r2.fastq.gz`,`sample.fastp.html`, and `sample.fastp.json`. If these are required the temporary() flag on the output files in the rule can be removed.


🔹 **`rule bowtie2_align` *Alignment to Host/Phix***
- **Purpose:** Aligns trimmed reads to a user created reference (Host/PhiX) that has been indexed by Bowtie2 index.
- **Inputs:** 
  - Trimmed paired reads: `*_r1.fastq.gz`, `*_r2.fastq.gz`
  - Bowtie2 index files with the suffix `.bt2`
- **Outputs:**
  - None. The reference-aligned `BAM` file is marked as temporary and automatically cleaned up once no longer needed by the pipeline.

- **Notes:**
  - Uses **default parameters** from `Bowtie2`.
  - This file is marked as temporary in the rule: `sample.bam`. If it is required the temporary() flag on the output file in the rule can be removed.
- **Performance Notes:**
  > **Wall time:**  
  > - 60 cores (bowtie2: 44, SAMtools view: 4, SAMtools sort: 12): 12m 18s  
  > - 24 cores (bowtie2: 16, SAMtools view: 4, SAMtools sort: 8): ??


🔹 **`rule extract_unmapped_fastq` *Decontamination***
- **Purpose:** extracts the reads that did not align into paired-end FASTQ files depleted of host and PhiX reads
- **Inputs:**
  - Sorted BAM file: `sample.bam`
- **Outputs:**
  - Clean read pairs: `sample_trimmed_clean_R1.fastq.gz`/`sample_trimmed_clean_R2.fastq.gz` 
- **Notes:**
  - Uses **default parameters** from `Bowtie2`.
- **Performance Notes:**
  >  **Wall time:**  
  > - 60 cores, no splitting: 17m 57s  
  > - Optimized run with core splitting at a 80:20 ratio between SAMtools and pigz (SAMtools: 48, pigz: 12) : time???

---  
### Module  `sortmerna.smk` contains this rule:
🔹 **`rule sortmerna` *rRNA Removal***
- **Purpose:** Align the clean read pairs to an rRNA database and outputs the rRNA-depleted reads
- **Inputs:**
  - Clean read pairs: `sample_trimmed_clean_R1.fastq.gz`/`sample_trimmed_clean_R2.fastq.gz` 
- **Outputs:**
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Notes:**
  - Uses **default parameters** from `SortMeRNA`
  - The database used for testing the pipeline was `smr_v4.3_default_db.fasta`, available from the Reference RNA databases (database.tar.gz) file at [sortmerna release v4.3.3](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3)

- **Performance Notes:**
  >  **Wall time:**  
  > - 60 cores wall time: 4h 48m 56s
  > - 48 cores wall time: ??

--- 
### Module  `taxonomy.smk` contains these rules:

🔹 **`rule kraken2` *Assign Taxonomy***
- **Purpose:** Assign taxonomy to the clean reads using a Kraken2-formatted GTDB
- **Inputs:**
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**
  - Kraken and report for each sample: `sample.kraken` and `sample.report.txt`
- **Notes:**
  - Uses confidence threshold of 0.5 and default parameters from `Kraken2`
  - New Kraken2 database in GPSC common has not been tested yet.
  - Devin Holman's database with additional references has been tested and worked when memory was increased to 840 GB
  - Must use **Large compute node**

- **Performance Notes:**
  >  **Wall time:**  
  > - Large compute node with 600 GB. With 16 CUPs wall time was 7m 56s
  > - Large compute node with 600 GB. With 2 CUPs wall time was 19m 13s

🔹 **`rule bracken` *Abundance Estimation***
  - **Purpose:** Refines Kraken classification to provide abundance estimates at the species, genus and phylum level for each sample.
  - **Inputs:** Kraken report: `sample.report.txt`
  - **Outputs:**  
  - Bracken reports at:
    - Species level: `sample_bracken.species.report.txt`
    - Genus level: `sample_bracken.genus.report.txt`
    - Phylum level: `sample_bracken.phylum.report.txt`
 - **Notes:**
  - Outputs are used as **intermediate files** for downstream rule: `combine_bracken_outputs`
  - his rule is also making `sample.report_bracken_species.txt` at each level in the `kraken2` directory. At some point see if we can either place these into a directory called `reports` or have them cleaned up in the shell block.

- **Performance Notes:**
  >  **Wall time:**  
  > - 10 threads the wall time was 9s.
  > - 2 threads ??

🔹 **`rule combine_bracken_outputs` *Merging Abundance Tables***
- **Inputs:**  
  - Bracken reports at species, genus, and phylum levels from `rule bracken`
- **Outputs:**  
  - Combined abundance tables for:
    - Species level: `merged_abundance_species.txt`
    - Genus level: `merged_abundance_genus.txt`
    - Phylum level: `merged_abundance_pylum.txt`

🔹 **`rule bracken_extract` *Relative Abundance Tables***
- **Purpose:** generate tables for the raw and relative abundance for each taxonomic level for all samples
- **Inputs:**
  - Combined abundance tables for:
    - Species level: `merged_abundance_species.txt`
    - Genus level: `merged_abundance_genus.txt`
    - Phylum level: `merged_abundance_pylum.txt`
- **Outputs:**
  - Combined relative and raw abundance tables for
    - Species level: `Bracken_species_raw_abundance.csv` and `Bracken_species_relative_abundance.csv`
    - Genus level: `Bracken_genus_raw_abundance.csv` and `Bracken_genus_relative_abundance.csv`
    - Phylum level: `Bracken_phylum_raw_abundance.csv` and `Bracken_genus_relative_abundance.csv`
---
### Module  `amr_short_reads.smk` contains these rules:
🔹 **`rule rgi_reload_database` *Load CARD DB***
- **Purpose:** Checks if the CARD Database has been loaded from a common directory or user specific directory
- **Inputs:** 
  - `card_reference.fasta`
  - `card.json`
- **Outputs:**
  - Done marker `rgi_reload_db.done` to prevent the rule from re-running every time the pipeline is called.

🔹 **`symlink_rgi_card` *Symlink CARD to the working directory***
- **Purpose:** Prevent the re-loading of the CARD DB

🔹 **`rule rgi_bwt` *Antimicrobial Resistance Gene Profiling***
- **Purpose:** performs antimicrobial resistance gene profiling on the cleaned reads using *k*-mer alignment (kma)
- **Inputs:** 
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**    
  - `sample_paired.allele_mapping_data.txt` – Text-formatted allele mapping  
  - `sample_paired.artifacts_mapping_stats.txt` – Statistics on mapping artifacts  
  - `sample_paired.gene_mapping_data.txt` – Per-gene alignment details  
  - `sample_paired.overall_mapping_stats.txt` – Summary statistics across all mappings  
  - `sample_paired.reference_mapping_stats.txt` – Reference-specific mapping stats  

- **Notes:**
  - Uses default RGI BWT parameters.
  - These files are marked as temporary in the rule: `sample_paired.allele_mapping_data.json`, `sample_paired.sorted.length_100.bam`, and `sample_paired.sorted.length_100.bam.bai`. If these are required the temporary() flag on the output files in the rule can be removed.

- **Performance Notes:**
  >  **Wall time:**  
  > - 40 cores wall time: 18m 7s
  > - 4 cores wall time: ??If time does not increase much further reduce cores.

---   
- **`rule coverm`**
   > **Note to self:** Add in the option of running CoverM. This should not be part of the main pipeline but an option if MAGs from metagenomic sequencing of the same samples are available.

- **`rules Cazymes`**
     > **Note to self:** Do we want to include this in the pipeline or use transcripts in existing Bash pipeline made by Arun.
---
### Module  `sample_assembly.smk` contains these rules:
🔹 **`rule rna_spades` *Assemble transcripts***
- **Purpose:** The rRNA-depleted reads are assembled into presumptive mRNA transcripts
- **Inputs:**
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs**
  - Presumptive transcripts: `sample.fasta`

- **Notes:**
  - Uses default SPAdes parameters with the -rna flag.

- **Performance Notes:**

  >  **Wall time:** 
  > -  One sample with 60 cores ran for 2h 15m 53s.
  > -  Try 48 and 32 cores to see if wall time is similar. Reduce cores if so.
  
🔹 **`rule rnaquast_busco` *QC for transcripts***
- **Purpose:** Reports the number of transcripts, transcripts over 500 bp, transcripts over 1000 bp and the BUSCO completeness.
- **Inputs:**
  - Presumptive transcripts: `sample.fasta`
  - Busco lineage: `bacteria_odb12` and `archaea_odb12`
- **Outputs:**
  - QUAST report in `sample_bacteria` and `sample_archaea` directories 

- **Notes:**
  - The software is not intended for metatranscriptomics. Use caution when interpreting the results. For instance the BUSCO completeness cannot be interpreted as the percentage of assembly quality but instead it is a representation of the core functions from the bacteria_odb12 and archaea_odb12 lineages.

--- 
### Module  `coassebly_annotation.smk` contains these rules:
 🔹 **`megahit_coassembly` *Co-assembly of all samples***
- **Purpose:** rRNA-depleted reads are co-assembled with MEGAHIT
- **Inputs:**
  - Cleaned sample reads from all samples: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**
  - Presumptive transcripts from the coassembly: `final.contigs.fa`

- **Notes:**
  - If Metagenomic sequencing was done the co-assembly of those reads would be a better choice.
  - The Co-assembly is used as a index to produce sorted BAM files for each assembly. These sorted BAM files can then be used in featureCounts and downstream expression analysis.

- **Performance Notes:**
  >  **Wall time:** 
  > - For the co-assembly of three samples was 52 min 28 sec with 60 cores. 
  > - For large co-assemblies a large mem node will need to be used (will need to be tested at some point)
  
 🔹 **`index_coassembly` *Create index***
- **Purpose:** Bowtie2 is used to make an index that can be used to map the reads to the co-assembly
- **Inputs:**
  - Presumptive transcripts from the coassembly: `final.contigs.fa`
- **Outputs:**
  - Bowtie2 index `coassembly.1.bt2`, `coassembly.2.bt2`, `coassembly.3.bt2`, `coassembly.4.bt2`, `coassembly.rev.1.bt2`, and `coassembly.rev.2.bt2`

- **Performance Notes:**
  >  **Wall time:** 
  > - Wall time with 8 threads was 1m 5s for a 185 kb assembly. 

🔹 **`bowtie2_map_transcripts` *Map samples to co-assembly***
- **Purpose:** Map the rRNA depleted cleaned reads from each sample to the co-assembly
- **Inputs:**
  - Bowtie2 index `coassembly.1.bt2`, `coassembly.2.bt2`, `coassembly.3.bt2`, `coassembly.4.bt2`, `coassembly.rev.1.bt2`, and `coassembly.rev.2.bt2`
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**
  - BAM file for each sample: `sample.coassembly.sorted.bam`

- **Performance Notes:**
  > **Wall time:**
  > - For one sample using 40 cores was 9m 53s. Reduced cores to 16 and will check wall time.

🔹 **`assembly_stats_depth` *QC and coverage for mapped reads***
- **Purpose:** `samtools flagstat` provides alignment statisitics that include the total reads, reads that mapped to the co-assembly, properly paired reads and duplicates. The flagstat is used to check how each sample aligns to the co-assembly. `samtools depth` computes the per-base sequencing depth across the co-assembly to evaluate sequncing depth and uniformity of coverage. `samtools idxstats` provides sequnces level mapping statistics with the sample contig name that is used to identify contigs that are over or under represented.
- **Inputs:**
  - BAM file for each sample: `sample.coassembly.sorted.bam`
**Outputs:**
  - Alignment statistics: `sample.flagstat.txt`
  - Sequencing depth: `sample.coverage.txt.gz`
  - Mapping statiscics: `sample.idxstats.txt.gz` 

🔹 **`rule prodigal_genes` *Gene prediction***
- **Purpose:** Predict the protein and nucleotide sequnces in the co-assembly. Generate a simplified annotation formate file that is used by `featurecounts`.
- **Inputs:**
  - Presumptive transcripts from the coassembly: `final.contigs.fa`
- **Outputs:**
  - Predicted protein sequnces: `coassembly.faa`
  - Predicted nuelotide sequnces: `coassembly.fna`
  - Feature formatted annotation file: `coassembly.gff`
  - Simplified annotation formate file: `coassembly.saf`

- **Notes:**
  - Go back and decided if this output should be designated temporary.

🔹 **`rule featurecounts` *Count table***
- **Purpose:** generates a table for each sample that includes the Geneid (unique identifier), the co-assembly contig name, the start and end positions of each gene on the contig, the strand orientation (+ or -), the gene length, and the number of reads mapped to each gene. Since all samples are mapped to the same co-assembly reference, the resulting tables can be combined for downstream analysis of gene expression across samples.
- **Inputs:**
  - Simplified annotation formate file: `coassembly.saf`
  - BAM file for each sample: `sample.coassembly.sorted.bam`
- **Outputs:**
  - Featurecounts table: `sample_counts.txt`
---

## Data

The raw input data must be in the form of paired-end FASTQ files generated from metatranscriptomics experiments.

- Each sample must include both forward (R1) and reverse (R2) read files.

**Example:**

- **Dataset 1 Filename**: Sequencing reads (FASTQ) from beef cattle rumen samples are provided for three samples: `LLC42Nov10C`, `LLC42Sep06CR`, and `LLC82Sep06GR`. There is also a sub-sampled test file, `test_LLC82Sep06GR`, which can be used for all steps dealing with un-assembled data but will fail for the RNA SPAdes step.

---

## Parameters

| Parameter          | Value                                                                                               |
| -------------------- | ----------------------------------------------------------------------------------------------------- |
| *samplesheet.csv* | *The samplesheet is described here: [Sample list](#sample-list)* |
| *preprocessing* | *The papameters are changed in the config/config.yaml. They are `cut_tail: true`, `cut_front: true`
  cut_mean_quality: 20
  cut_window_size: 4
  qualified_quality_phred: 15
  detect_adapter_for_pe: true
  length_required: 100*              |
| *parameter_name_2* | *Description of what the parameter does and the expected value (e.g., boolean, list).*              |

---

## Usage

### Pre-requisites

#### Software

- Snakemake version 9.6.0
- Snakemake-executor-plugin-slurm

#### Databases

- **Bowtie2**  
  Bowtie2 uses an index of reference sequences to align reads. This index must be created before running the pipeline. The index files (with the `.bt2` extension) must be located in the directory you specify in the `config/config.yaml` file. Make sure to update the prefix of these files in the `config.yaml` file.

- **SortMeRNA**  
  SortMeRNA requires a ribosomal (r)RNA database in the `rRNA_DB` directory. Update the `config.yaml` file with the filename of the database used. You can download the database from [SortMeRNA releases](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3). The file `smr_v4.3_default_db.fasta` was used for pipeline testing.

- **Kraken2**  
  Kraken2 requires a Kraken2-formatted GTDB database. The GTDB release tested with this pipeline was 220.

- **RGI BWT/CARD**  RGI BWT requires the CARD (Comprehensive Antibiotic Resistance Database) database. The version tested in this pipeline was 4.0.1. The database can be located on a common drive or in your working directory.  
  Instructions for installing the CARD database are available on [CARD RGI github](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst).  
  Steps copied from the RGI documentation:

  **Download CARD data:**

  ```bash
  wget https://card.mcmaster.ca/latest/data
  tar -xvf data ./card.json

  rgi load --card_json /path/to/card.json --local

  rgi card_annotation -i /path/to/card.json > card_annotation.log 2>&1

  rgi load -i /path/to/card.json --card_annotation card_database_v3.0.1.fasta --local
  ```

  **Note:** the files after loading and annotating card must be called `card.json` and `card_reference.fasta`

- **BUSCO**
rnaQUAST uses the BUSCO bacterial and archaeal lineages. The directory path to these lineages must be provided.

### Setup Instructions

#### 1. Installation

Clone the repository into the directory where you want to run the metatranscriptomics Snakemake pipeline.  
**Note:** This location must be on an HPC (High Performance Computing) cluster with access to a high-memory node (at least 600 GB RAM) and sufficient storage for all metatranscriptomics analyses.

```bash
cd /path/to/code/directory
git clone <repository-url>
```
#### 2. SLURM Profile
##### 2.1. SLURM Profile Directory Structure
```
metatranscriptomics_pipeline/
├── Workflow/
│   └── Snakemake/
│   └── ... 
├── profiles/
│   └── slurm/
│       └── config.yaml         ← profile config
├── config/
│   └── config.yaml             ← workflow data/sample config
|   └── samples.txt
├── run_snakemake.sh            ← your SLURM launcher
├── .env
└── ...                         
```
##### 2.2. Profile Configuration
The SLURM execution settings are configured in profiles/slurm/config.yaml. This file defines resource defaults, cluster submission commands, and job script templates for Snakemake. The pre-rule resources need to be adjusted for the size and number of input samples for each rule.

**Example for profiles/slurm/config.yaml:**
```bash
### How Snakemake assigns resources to rules
cores: 60
jobs: 10 
latency-wait: 60 
rerun-incomplete: true
retries: 2              
max-jobs-per-second: 2 
executor: slurm


### Env Vars ###
envvars:
  TMPDIR: "/gpfs/fs7/aafc/scratch/$USER/tmpdir"

default-resources:
  - slurm_account=aafc_aac
  - slurm_partition=standard
  - slurm_cluster=gpsc8
  - runtime=60       # minutes
  - mem_mb=4000
  - cpus=1

### Env modules ###
# use-envmodules: false 

### Conda ###
use-conda: true
conda-frontend: mamba   

### Resource scopes ###
set-resource-scopes:
  cores: local 

## Per rule resources
set-resources:
  fastp_pe:
    cpus: 2
    mem_mb: 4000
    runtime: 40
    slurm_partition: standard
    slurm_account: aafc_aac
    slurm_cluster: gpsc8

  bowtie2_align:
    cpus: 24
    mem_mb: 48000
    runtime: 30
    slurm_partition: standard
    slurm_account: aafc_aac
    slurm_cluster: gpsc8
```

#### 3. Configuration

The pipeline requires the following configuration files: `config.yaml`, `.env`, and `samples.txt`.

##### 3.1. config/config.yaml

The `config.yaml` file must be located in the `config` directory, which resides in the main Snakemake working directory. This file specifies crucial settings, including:

- Path to the `samples.txt`  
- Input and output directories  
- File paths to required databases 
- Parameters for each rule **NEED TO UPDATE RULES** 

**Note:**  
You must edit `config.yaml` **before** running the pipeline to ensure all paths are correctly set.  
For best practice, use database paths that are in common locations to all users on the HPC.

##### 3.2. Environment file

This file must contain paths to the **PROJECT ROOT**,  **USER SCRATCH**, and **RGI COMMON DATABASE**. Follow these instructions:

- In the main Snakemake directory (where you are running Snakemake from)

```bash
touch .env
```

- Open the .env file and add

```bash
 PROJECT_ROOT = path/to/project/root
 TMPDIR = path/to/temp/on/cluster **Issue with $USER. I had to use my actual username in the .env file**
 RGI_CARD = path/to/card.json and card_reference.fasta
 ```

##### 3.3. Sample list

 `samplesheet.csv` Has the following column names: "sample","fastq_1","fastq_2". For the column 'sample" use the sampleID for the read pair, and for "fastq_1","fastq_2" have the names of the read1 and read2 files as they appear in the raw fastq files directory. The file location of the `samplesheet.csv` must be`config/samplesheet.csv`.

**Example `samplesheet.csv`:**
sample,fastq_1,fastq_2
test_LLC82Nov10GR,test_LLC82Nov10GR_r1.fastq.gz,test_LLC82Nov10GR_r2.fastq.gz
test_LLC82Sep06GR,test_LLC82Sep06GR_r1.fastq.gz,test_LLC82Sep06GR_r2.fastq.gz

#### 4. Running the pipeline

Complete steps **1.Installation**, **2.SLURM Profile**, and **3.Configuration** and ensure database paths have been added to the 'config/config.yaml'. Required databases are described in the [Pre-requisites](#pre-requisites).

##### 4.1. Conda environments

Snakemake can automatically create and load Conda environments for each rule in your workflow. Check to see that you have the following configuration files in the `workflow/envs` directory:

- `bedtools.yaml`
- `bowtie2.yaml`
- `featurecounts.yaml`
- `kraken2.yaml`
- `megahit.yaml`
- `rgi.yaml`
- `rnaquast.yaml`
- `RNAspades.yaml`
- `sortmerna.yaml`

Load the required conda environments for the pipeline with:

```bash
snakemake --use-conda \
  --conda-create-envs-only \
  --conda-prefix path/to/common/lab/folder/conda/metatranscriptomics-snakemake-conda
```
##### 4.2. SLURM launcher
This is the script you use to submit the Snakemake pipeline to SLURM.
- **Before submitting job to SLURM run `export SLURM_CONF="/etc/slurm-llnl/gpsc8.science.gc.ca.conf"`**
- Defines resources for the job scheduler
- Activates the Snakemake environment
- Submits and manages jobs using the Snakemake `--profile` configuration `(profiles/slurm/)`.
- Contains any additional Snakemake arguments (e.g.., `--unlock`, `--dry-run`, `--rerun-incomplete`)
- For a snakemake report with runtime and software versions use --report path/to/metatranscriptomics_report.html after the pipeline has completed

```bash
#!/bin/bash
#SBATCH --job-name=run_snakemake.sh
#SBATCH --output=run_snakemaket_%j.out 
#SBATCH --error=run_snakemake_%j.err 
#SBATCH --cluster=gpsc8 
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --mem=2000
#SBATCH --time=8:00:00

source /gpfs/fs7/aafc/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake-9.6.0
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile absolute/path/to/profiles/slurm \
    --configfile absolute/path/to/config/config.yaml \
    --conda-prefix absolute/path/to/common/conda/metatranscriptomics-snakemake-conda \
    --printshellcmds \
    --keep-going 
  ```
### Notes
- temp folder is set to `/gpfs/fs7/aafc/scratch/$USER/tmpdir` for running on the GPSC.
#### Warnings

- The conda environments will not be created if the conda configuration is `conda config --set channel_priority strict`.
- Set conda to `conda config --set channel_priority flexible` or use libmamba.
- The `.env` file can overwrite the `config/config.yaml` file

#### Current issues

- In the .env file /gpfs/fs7/aafc/scratch/$USER/ was not solving to user so as a temporary fix I put in my user name.
- When poor sample reads are used in the pipeline rna SPAdes cannot make an transcripts.fasta file. A temporary solution is to make a dummy fasta file. This results in failed downstream rules for rnaQUAST.

#### Resource usage

- Kraken2: Large compute node with 600 GB. With 16 CUPs wall time was 7m 56s. With 2 CPUs wall time was 19m 13s.

## OUTPUT

*Provide format, location, and naming of result files, and a brief description.*

*Example Output:*

*Output files include:*

*- results/reports/summary.csv: Key metrics from analysis.*

*- results/logs/pipeline.log: Step-by-step log.*

*- results/plots/visualization.png: Output plot.*

---

## Credits

This repository was written by Katherine James-Gzyl and assisted by Devin Holman and Arun Kommadath.

---

## Contribution

If you would like to contribute to this project, please review the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md) and ensure you adhere to our [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) to foster a respectful and inclusive environment.

---

## License

This project is distributed under the MIT License. For complete details and copyright information, see the [LICENSE](LICENSE) file.

---

## References

*Provide references to key publications and any useful resources for tools/software used. Formal citations of the tools used may also be provided via a CITATIONS.md file.*

*Example References:*

### Publications

The pipeline and analysis associated with it is published here:

- Your published paper title – Journal, Year.

### Resources

- Link to Snakemake Manual
- Link to Tool X Documentation
  
### Tools/Software

References to tools and software used here can be found in the [CITATIONS.md](CITATIONS.md) file.

## Citation

*Provide information on how to cite this repository. Use a CITATION.cff file where required. Citation tools like GitHub and Zenodo will use this file to generate standardized references.*

If you use this project in your work, please cite it using the [CITATION.cff](CITATION.cff) file.
