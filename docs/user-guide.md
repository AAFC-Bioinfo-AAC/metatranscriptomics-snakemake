<!-- omit in toc -->
# METATRANSCRIPTOMICS SNAKEMAKE PIPELINE - USER GUIDE

---

<!-- omit in toc -->
## Table of Contents

- [Overview](#overview)
  - [Workflow diagram](#workflow-diagram)
  - [Snakemake rules](#snakemake-rules)
    - [Module `preprocessing.smk`](#module-preprocessingsmk)
    - [Module `sortmerna.smk`](#module-sortmernasmk)
    - [Module `taxonomy.smk`](#module-taxonomysmk)
    - [Module `amr_short_reads.smk`](#module-amr_short_readssmk)
    - [Module `sample_assembly.smk`](#module-sample_assemblysmk)
    - [Module `coassembly_annotation.smk`](#module-coassembly_annotationsmk)
    - [Module `db_can.smk`](#module-db_cansmk)
    - [Module `env_versions.smk`](#module-env_versionssmk)
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
      - [3.4. Scripts called in rules](#34-scripts-called-in-rules)
    - [4. Running the pipeline](#4-running-the-pipeline)
      - [4.1. Conda environments](#41-conda-environments)
      - [4.2. SLURM launcher](#42-slurm-launcher)
      - [4.3. Submit launcher to SLURM](#43-submit-launcher-to-slurm)
  - [Notes](#notes)
    - [Warnings](#warnings)
    - [Current issues](#current-issues)
    - [Resource usage](#resource-usage)
- [Output](#output)

---

## Overview

### Workflow diagram

```mermaid
flowchart TD
    subgraph PREPROC["READ PROCESSING"]
        A["Paired-end RNA reads"] --> B{"fastp"}
        B --> C["Trimmed paired reads"]
        B --> L["Retained HTML/JSON reports"]
        C --> D{"Bowtie2 host/PhiX alignment"}
        D --> HBAM["Sorted BAM"]
        HBAM --> EXTRACT{"SAMtools and BEDTools"}
        EXTRACT --> E["Pairs with both mates unmapped"]
        E --> F{"SortMeRNA"}
        F --> G["rRNA-filtered paired reads"]
    end

    subgraph SHORTREAD["SHORT-READ ANALYSIS"]
        G --> K{"Kraken2"}
        K --> BR{"Bracken"}
        BR --> TAX["RNA taxonomic profiles"]
        G --> RGI{"RGI BWT"}
        RGI --> ARG["ARG-associated RNA profiles"]
    end

    subgraph ASSEMBLY["INDIVIDUAL SAMPLE ASSEMBLY"]
        G --> SP{"rnaSPAdes"}
        SP --> TRANS["Sample transcript assemblies"]
        TRANS --> RQ{"rnaQUAST with BUSCO"}
        RQ --> QC["Assembly evaluation reports"]
    end

    subgraph SHARED["SHARED REFERENCE AND GENE QUANTIFICATION"]
        G -->|No external reference| MH{"MEGAHIT"}
        MH --> REF["Shared reference FASTA"]
        DNA["Matched metagenomic assembly"] -->|reference_assembly set| REF

        REF --> INDEX{"Bowtie2-build"}
        INDEX --> IDX["Reference index"]
        IDX --> MAP{"Bowtie2 RNA mapping"}
        G --> MAP
        MAP --> BAM["Sorted and indexed sample BAMs"]

        BAM --> STATS{"SAMtools"}
        STATS --> STATSOUT["Mapping and depth statistics"]

        REF --> PG{"Prodigal"}
        PG --> PROTEINS["Proteins and CDS sequences"]
        PG --> ANNOT["GFF and SAF annotations"]

        ANNOT -->|SAF| FC{"featureCounts"}
        BAM --> FC
        FC --> COUNTS["Raw paired-fragment counts"]
    end

    subgraph CAZYME["CAZYME ANNOTATION AND RNA COUNTS"]
        PROTEINS -->|Proteins| PREP{"Match protein and gene IDs"}
        ANNOT -->|GFF| PREP
        PREP --> CAZINPUT["Protein FASTA and ID map"]

        CAZINPUT -->|Protein FASTA| DBCAN{"run_dbcan CAZyme_annotation"}
        DBCAN --> CAZANNOT["CAZyme annotations"]

        CAZANNOT --> SUMMARY{"Filter annotations and summarize counts"}
        CAZINPUT -->|ID map| SUMMARY
        COUNTS --> SUMMARY
        SUMMARY --> MATRICES["CAZyme gene and family count matrices"]
    end

    classDef default fill:#1f2937,stroke:#F8B229,color:#e5e7eb;
    classDef temporary fill:#1f2937,stroke:#22d3ee,stroke-dasharray:5 5,color:#e5e7eb;
    class C,HBAM temporary;
    linkStyle default stroke:#F8B229,stroke-width:1.5px;
```

### Snakemake rules

The pipeline is modularized, with each module located in the `metatranscriptomics-snakemake/workflow/rules` directory. The modules are `preprocessing.smk`, `sortmerna.smk`, `taxonomy.smk`,`amr_short_reads.smk`, `sample_assembly.smk`, `coassembly_annotation.smk`, and `env_versions.smk`.

---

#### Module `preprocessing.smk`

This module performs read-quality control, trimming and removal of reads originating from the host or PhiX control.

**Default configuration settings**

The following values are supplied in `config/config.yaml`. These workflow settings may differ from the default settings used by the individual software packages.

| Configuration setting | Default | Description |
|---|---:|---|
| `fastp: threads` | `2` | Number of threads used by *fastp*. |
| `fastp: cut_tail` | `true` | Enables sliding-window quality trimming from the 3′ end. |
| `fastp: cut_front` | `true` | Enables sliding-window quality trimming from the 5′ end. |
| `fastp: cut_mean_quality` | `20` | Minimum mean Phred quality required within the trimming window. |
| `fastp: cut_window_size` | `4` | Number of bases included in the sliding quality window. |
| `fastp: qualified_quality_phred` | `15` | Minimum Phred score used to define a qualified base. |
| `fastp: detect_adapter_for_pe` | `true` | Enables automatic adapter detection for paired-end reads. |
| `fastp: length_required` | `100` | Minimum read length retained after trimming. |
| `bowtie2_align: threads` | `12` | Total number of threads allocated among *Bowtie2* and *SAMtools*. |
| `extract_unmapped_fastq: threads` | `8` | Total number of threads allocated among *SAMtools*, *BEDTools*, and two *pigz* compressors. |

Quality-filtering parameters should be selected according to the sequencing platform, read length and study objectives.

**Rule: `fastp_pe` Quality Control & Trimming**

- **Purpose:** Uses *fastp* to perform adapter detection, adapter trimming, quality trimming, quality filtering, and length filtering of paired-end reads.
- **Inputs:**
- Paired-end fastq files specified for each sample in `samplesheet.csv`.
- **Outputs:**
  - Trimmed R1 reads: `sample_r1.fastq.gz`
  - Trimmed R2 reads: `sample_r2.fastq.gz`
  - Unpaired R1 reads: `sample_u1.fastq.gz`
  - Unpaired R2 reads: `sample_u2.fastq.gz`
  - HTML quality-control report: `sample.fastp.html`
  - JSON quality-control report: `sample.fastp.json`
- **Notes:**
  - Only the paired trimmed reads are used by subsequent preprocessing rules.
  - The four fastq outputs from this rule are marked with the Snakemake `temp()` function and are removed automatically when they are no longer required.
  - The HTML and JSON quality-control reports are retained automatically in the `fastp` subdirectory beneath the configured log directory.
  - To retain intermediate FASTQ files, remove `temp()` from the corresponding outputs in `workflow/rules/preprocessing.smk`.
  - The processing log is written to `fastp/sample.fastp.log` beneath the configured log directory.


**Rule: `bowtie2_align` — Host and PhiX alignment**

- **Purpose:** Aligns the trimmed paired reads against a user-supplied Bowtie2 index containing the relevant host genome sequence and PhiX reference sequence. The alignments are converted into a coordinate-sorted BAM file using *SAMtools*.
- **Inputs:**
  - Trimmed R1 reads: `sample_r1.fastq.gz`
  - Trimmed R2 reads: `sample_r2.fastq.gz`
  - Complete Bowtie2 index specified by `bowtie2_index` in `config/config.yaml`
- **Outputs:**
  - Coordinate-sorted alignment file: `bam/sample.bam`
- **Notes:**
  - The workflow currently declares the six `.bt2` index files; `.bt2l` large indexes are not currently supported by its input declarations.
  - Bowtie2 default alignment settings are used apart from the configured number of threads and the addition of a read-group ID and sample tag.
  - Available threads are divided among *Bowtie2* and *SAMtools*. The rule requires at least three allocated threads.
  - The BAM file is marked with `temp()` because it is an intermediate file used to recover host/PhiX-depleted paired reads.
  - The BAM contains both mapped and unmapped records.
  - The alignment log is written to `bowtie2/sample.log` beneath the configured log directory.

**Rule: `extract_unmapped_fastq` — Host and PhiX read removal**

- **Purpose:** Extracts paired reads for which neither mate aligned to the combined host and PhiX reference.
- **Inputs:**
  - Coordinate-sorted alignment file: `bam/sample.bam`
- **Outputs:**
  - Host/PhiX-depleted R1 reads: `sample_trimmed_clean_R1.fastq.gz`
  - Host/PhiX-depleted R2 reads: `sample_trimmed_clean_R2.fastq.gz`
- **Notes:**
  - *SAMtools* retains read pairs for which both mates are unmapped and excludes secondary and supplementary alignments.
  - Pairs with one mapped mate are excluded, even if the other mate is unmapped.
  - The retained BAM records are sorted by read name before *BEDTools* converts them back into paired fastq files.
  - The fastq files are compressed using *pigz*.
  - Available threads are divided among *SAMtools*, *BEDTools* and two *pigz* compressors. The rule requires at least five allocated threads.
  - Temporary sorting files are written beneath `TMPDIR` or `/tmp` if `TMPDIR` is unset, and removed when the rule finishes.
  - The host/PhiX-depleted fastq files are retained as standard outputs; this rule does not use `protected()`.
  - These host/PhiX-depleted paired reads are inputs to the SortMeRNA module. Its rRNA-filtered paired outputs are used for downstream taxonomic profiling, ARG profiling, assembly and RNA mapping.
  - The extraction log is written to `bedtools/sample.log` beneath the configured log directory.
---

#### Module `sortmerna.smk`

This module computationally filters rRNA from host/PhiX-depleted paired reads.

**Default configuration settings**

The following value is supplied in `config/config.yaml`.

| Configuration setting | Default | Description |
|---|---:|---|
| `sortmerna_pe: threads` | `12` | Number of threads allocated to SortMeRNA. |

**Rule: `sortmerna_pe` — rRNA filtering**

- **Purpose:** Uses *SortMeRNA* to align host/PhiX-depleted paired reads against a configured rRNA reference database and retain pairs for which neither mate meets the rRNA alignment criteria.
- **Inputs:**
  - Host/PhiX-depleted read pairs: `sample_trimmed_clean_R1.fastq.gz` / `sample_trimmed_clean_R2.fastq.gz`
  - rRNA reference fasta file specified by `sortmerna_DB` in `config/config.yaml`.
- **Outputs:**
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz` / `sample_rRNAdep_R2.fastq.gz`
  - Filtering statistics: `sample_sortmerna_pe.stats`
- **Notes:**
  - The `--paired_in` option excludes both mates from the retained output if either mate matches the rRNA reference. This preserves pairing but can also remove a non-rRNA mate paired with an rRNA read.
  - The `--out2` option writes retained R1 and R2 reads into separate compressed fastq files.
  - The database used for testing the pipeline was `smr_v4.3_default_db.fasta`, available from the Reference RNA databases archive (`database.tar.gz`) at [SortMeRNA release v4.3.3](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3).
  - The Conda environment specifies SortMeRNA version `4.3.6`; the database download release and software version are separate.
  - Each sample uses a unique temporary working directory beneath `TMPDIR`, or `/tmp` if `TMPDIR` is unset. The reference index is built within that directory rather than reused from a shared index.
  - The statistics file contains the SortMeRNA alignment summary. The retained fastq files and statistics file are saved before the temporary working directory is removed.
  - The processing log is written to `sortmerna/samplereads_pe.log` beneath the configured log directory.
  - The retained paired reads are used for downstream taxonomic profiling, ARG profiling, assembly and RNA mapping.

---

#### Module `taxonomy.smk`

**`rule kraken2` *Assign Taxonomy***

- **Purpose:** Assign taxonomy to the clean reads using a Kraken2-formatted GTDB
- **Inputs:**
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**
  - Kraken and report for each sample: `sample.kraken` and `sample.report.txt`
- **Notes:**
  - Must use **Large compute node** with at least 600 GB.

**`rule bracken` *Abundance Estimation***

- **Purpose:** Refines Kraken classification to provide abundance estimates at the species, genus and phylum level for each sample.
- **Inputs:** Kraken report: `sample.report.txt`
- **Outputs:**
  - Bracken reports at:
    - Species level: `sample_bracken.species.report.txt`
    - Genus level: `sample_bracken.genus.report.txt`
    - Phylum level: `sample_bracken.phylum.report.txt`
- **Notes:**
  - Outputs are used as **intermediate files** for downstream rule: `combine_bracken_outputs`
  - This rule is also making `sample.report_bracken_species.txt` at each level in the `kraken2` directory. At some point see if we can either place these into a directory called `reports` or have them cleaned up in the shell block.

**`rule combine_bracken_outputs` *Merging Abundance Tables***

- **Inputs:**
  - Bracken reports at:
    - Species level: `sample_bracken.species.report.txt`
    - Genus level: `sample_bracken.genus.report.txt`
    - Phylum level: `sample_bracken.phylum.report.txt`
- **Outputs:**
  - Combined abundance tables for:
    - Species level: `merged_abundance_species.txt`
    - Genus level: `merged_abundance_genus.txt`
    - Phylum level: `merged_abundance_phylum.txt`

**`rule bracken_extract` *Relative Abundance Tables***

- **Purpose:** generate tables for the raw and relative abundance for each taxonomic level for all samples
- **Inputs:**
  - Combined abundance tables for:
    - Species level: `merged_abundance_species.txt`
    - Genus level: `merged_abundance_genus.txt`
    - Phylum level: `merged_abundance_phylum.txt`
- **Outputs:**
  - Combined relative and raw abundance tables for
    - Species level: `Bracken_species_raw_abundance.csv` and `Bracken_species_relative_abundance.csv`
    - Genus level: `Bracken_genus_raw_abundance.csv` and `Bracken_genus_relative_abundance.csv`
    - Phylum level: `Bracken_phylum_raw_abundance.csv` and `Bracken_genus_relative_abundance.csv`

---

#### Module `amr_short_reads.smk`

**`rule rgi_reload_database` *Load CARD DB***

- **Purpose:** Checks if the CARD Database has been loaded from a common directory or user specific directory
- **Inputs:**
  - `card_reference.fasta`
  - `card.json`
- **Outputs:**
  - Done marker `rgi_reload_db.done` to prevent the rule from re-running every time the pipeline is called.

**`symlink_rgi_card` *Symlink CARD to the working directory***

- **Purpose:** Prevent the re-loading of the CARD DB

**`rule rgi_bwt` *Antimicrobial Resistance Gene Profiling***

- **Purpose:** performs antimicrobial resistance gene profiling on the cleaned reads using *k*-mer alignment (kma)
- **Inputs:**

  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**

  - `sample_paired.allele_mapping_data.txt`
  - `sample_paired.artifacts_mapping_stats.txt`
  - `sample_paired.gene_mapping_data.txt`
  - `sample_paired.overall_mapping_stats.txt`
  - `sample_paired.reference_mapping_stats.txt`
- **Notes:**

  - Uses default RGI BWT parameters.
  - For large sample files the large memory node may be required.
  - These files are marked as temporary in the rule: `sample_paired.allele_mapping_data.json`, `sample_paired.sorted.length_100.bam`, and `sample_paired.sorted.length_100.bam.bai`. If these are required the temporary() flag on the output files in the rule can be removed.

---

#### Module `sample_assembly.smk`

**`rule rna_spades` *Assemble transcripts***

- **Purpose:** The rRNA-depleted reads are assembled into presumptive mRNA transcripts
- **Inputs:**

  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs**

  - Presumptive transcripts: `sample.fasta`
- **Notes:**

  - Poor quality samples that result in no assembly hav an empty sample.fasta file

  **`rule rnaquast_busco` *QC for transcripts***
- **Purpose:** Reports the number of transcripts, transcripts over 500 bp, transcripts over 1000 bp and the BUSCO completeness.
- **Inputs:**

  - Presumptive transcripts: `sample.fasta`
  - Busco lineage: `bacteria_odb12` and `archaea_odb12`
- **Outputs:**

  - QUAST report in `sample_bacteria` and `sample_archaea` directories
- **Notes:**

  - The software is not intended for metatranscriptomics. Use caution when interpreting the results. For instance the BUSCO completeness cannot be interpreted as the percentage of assembly quality but instead it is a representation of the core functions from the bacteria_odb12 and archaea_odb12 lineages.

---

#### Module `coassembly_annotation.smk`

**`megahit_coassembly` *Co-assembly of all samples***

- **Purpose:** rRNA-depleted reads are co-assembled with MEGAHIT
- **Inputs:**

  - Cleaned sample reads from all samples: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**

  - Presumptive transcripts from the coassembly: `final.contigs.fa`
- **Notes:**

  - If Metagenomic sequencing was done the co-assembly of those reads would be a better choice.
  - The Co-assembly is used as a index to produce sorted BAM files for each assembly. These sorted BAM files can then be used in featureCounts and downstream expression analysis.

**`index_coassembly` *Create index***

- **Purpose:** Bowtie2 is used to make an index that can be used to map the reads to the co-assembly
- **Inputs:**
  - Presumptive transcripts from the coassembly: `final.contigs.fa`
- **Outputs:**
  - Bowtie2 index `coassembly.1.bt2`, `coassembly.2.bt2`, `coassembly.3.bt2`, `coassembly.4.bt2`, `coassembly.rev.1.bt2`, and `coassembly.rev.2.bt2`

**`bowtie2_map_transcripts` *Map samples to co-assembly***

- **Purpose:** Map the rRNA depleted cleaned reads from each sample to the co-assembly
- **Inputs:**
  - Bowtie2 index `coassembly.1.bt2`, `coassembly.2.bt2`, `coassembly.3.bt2`, `coassembly.4.bt2`, `coassembly.rev.1.bt2`, and `coassembly.rev.2.bt2`
  - rRNA-depleted reads: `sample_rRNAdep_R1.fastq.gz`/`sample_rRNAdep_R2.fastq.gz`
- **Outputs:**
  - BAM file for each sample: `sample.coassembly.sorted.bam`

**`assembly_stats_depth` *QC and coverage for mapped reads***

- **Purpose:** `samtools flagstat` provides alignment statistics that include the total reads, reads that mapped to the co-assembly, properly paired reads and duplicates. The flagstat is used to check how each sample aligns to the co-assembly. `samtools depth` computes the per-base sequencing depth across the co-assembly to evaluate sequencing depth and uniformity of coverage. `samtools idxstats` provides sequences level mapping statistics with the sample contig name that is used to identify contigs that are over or under represented.
- **Inputs:**
  - BAM file for each sample: `sample.coassembly.sorted.bam`
    **Outputs:**
  - Alignment statistics: `sample.flagstat.txt`
  - Sequencing depth: `sample.coverage.txt.gz`
  - Mapping statistics: `sample.idxstats.txt.gz`

**`rule prodigal_genes` *Gene prediction***

- **Purpose:** Predict the protein and nucleotide sequences in the co-assembly. Generate a simplified annotation format file that is used by `featurecounts`.
- **Inputs:**

  - Presumptive transcripts from the coassembly: `final.contigs.fa`
- **Outputs:**

  - Predicted protein sequences: `coassembly.faa`
  - Predicted nucleotide sequences: `coassembly.fna`
  - Feature formatted annotation file: `coassembly.gff`
  - Simplified annotation format file: `coassembly.saf`
- **Notes:**

  - Go back and decided if this output should be designated temporary.

**`rule featurecounts` *Count table***

- **Purpose:** generates a table for each sample that includes the Geneid (unique identifier), the co-assembly contig name, the start and end positions of each gene on the contig, the strand orientation (+ or -), the gene length, and the number of reads mapped to each gene. Since all samples are mapped to the same co-assembly reference, the resulting tables can be combined for downstream analysis of gene expression across samples.
- **Inputs:**
  - Simplified annotation formate file: `coassembly.saf`
  - BAM file for each sample: `sample.coassembly.sorted.bam`
- **Outputs:**
  - Featurecounts table: `sample_counts.txt`

#### Module `env_versions`

**`software_report` *Print versions of all conda packages***

**`filter_key_bioinformatics_versions` *Print versions of only the core bioinformatic conda packages***

**`filter_key_bioinformatics_html` *Print versions of only the core bioinformatic conda packages in the html report***

---

## Data

The raw input data must be in the form of paired-end FASTQ files generated from metatranscriptomics experiments.

A set of sub-sampled raw FASTQ files are provided for testing (`data/test_LLC82Sep06GR_{r1,r2}*.fastq.gz`). These files can be used for all workflow stages that operate on un-assembled data (e.g., quality control, filtering, mapping, and annotation). For the `RNA-SPAdes` step, the test sample (or any poor/low coverage samples) will produce an empty fasta file to let Snakemake know that the file output has been created. For such samples, the next step in the pipeline (rnaQUAST) will be skipped.

---

## Parameters

The `config/config.yaml` file contains the editable pipeline parameters, thread allocation for rules with more than one core, and the relative file paths for input and output. The prefix of the absolute file path must go in `.env`. Most tools in the pipeline have default parameters. The tools with parameters different from default or that can be edited in the `config/config.yaml` file are listed below.

| Parameter                       | Value                                                                                                                                                                    |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| samplesheet.csv                 | The samplesheet is described here:[Sample list](#33-sample-list)                                                                                                         |
| fastp:*cut_tail*                | If true, trim low quality bases from the 3′ end until a base meets or exceeds the*cut_mean_quality* threshold. If false,disabled.                                       |
| fastp:*cut_front*               | If true, trim low quality bases from the 5′ end until a base meets or exceeds the*cut_mean_quality* threshold. If false,disabled.                                       |
| fastp:*cut_mean_quality*        | A positive integer specifying the minimum average quality score threshold for sliding window trimming.                                                                   |
| fastp:*cut_window_size*         | A positive integer specifying the sliding window size in bp when using*cut_mean_quality*.                                                                                |
| fastp:*qualified_quality_phred* | A positive integer specifying the minimum Phred score that a base needs to be considered qualified.                                                                      |
| fastp:*detect_adapter_for_pe*   | If true, auto adapter detection. If false,disabled.                                                                                                                      |
| fastp:*length_required*         | Reads shorter then this positive integer will be discarded.                                                                                                              |
| kraken2:*conf_threshold*        | Interval between 0 and 1. Higher values require more of a read’s k-mers to match the same taxon before it is classified, increasing precision but reducing sensitivity. |
| bracken:*readlen*               | The read length of your data in bp.                                                                                                                                      |
| rna_spades:*memory*             | Memory limit set in mb.                                                                                                                                                  |

---

## Usage

### Pre-requisites

#### Software

- Snakemake version 9.6.0 *Most rules also tested with Snakemake version 9.9.0.*
- Snakemake-executor-plugin-slurm version 1.6.1. *Earlier version resulted in SLURM communication issues*

#### Databases

- **Bowtie2**
  Bowtie2 uses an index of reference sequences to align reads. This index must be created before running the pipeline. The index files (with the `.bt2` extension) must be located in the directory specified in `config/config.yaml`. Make sure to update the prefix of these files in the `config.yaml` file. For instructions on creating the index please see the [Bowtie2 GitHub repository](https://github.com/BenLangmead/bowtie2).
- **SortMeRNA**
  SortMeRNA requires a ribosomal (r)RNA database in the `rRNA_DB` directory. Update the `config.yaml` file with the filename of the database used. You can download the database from [SortMeRNA releases](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3). The file `smr_v4.3_default_db.fasta` was used for pipeline testing.
- **Kraken2**
  Kraken2 requires a Kraken2-formatted GTDB. The GTDB release tested with this pipeline was 220. Pre-built Kraken2-formatted GTDB are available from [Kraken 2, KrakenUniq and Bracken indexes](https://benlangmead.github.io/aws-indexes/k2), and instructions for building custom Kraken2-formatted GTDBs are available on the [Kraken2 GitHub repository](https://github.com/DerrickWood/kraken2).
- **RGI BWT/CARD**  RGI BWT requires the CARD (Comprehensive Antibiotic Resistance Database) database. The version tested in this pipeline was 4.0.1. The database can be located on a common drive or in your working directory.
  Instructions for installing the CARD database are available on [CARD RGI GitHub repository](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst).
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
  rnaQUAST uses the BUSCO bacterial and archaeal lineages. The directory path to these lineages must be provided in `config/config.yaml`. The [BUSCO lineages](https://busco.ezlab.org/busco_userguide.html#lineage-datasets) are available on the webpage.

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

```bash
metatranscriptomics_pipeline/
├── Workflow/
│   └── Snakefile
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

The SLURM execution settings must be configured in `profiles/slurm/config.yaml.` An editable example is provided in this repository at `profiles/slurm/example_config.yaml` After editing, rename this file to `config.yaml` so that Snakemake recognizes it.

- This configuration file defines resource defaults, cluster submission commands, and job script templates for Snakemake. It should be customized for each specific HPC environment.
- Remember to update the rerun-triggers: [input, params, software-env] setting whenever the pipeline is modified.
- Pre-rule resource allocations should also be adjusted according to the size and number of input samples for each rule.

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

# Prevent rerunning jobs just for Snakefile edits
## flags available [input, mtime, params, software-env, code, resources, none]
rerun-triggers: [input, params, software-env]

### Env Vars ###
envvars:
  TMPDIR: "/path/to/scratch/${USER}/tmpdir"

default-resources:
  - slurm_account=<ACCOUNT_NAME>
  - slurm_partition=<PARTITION_NAME>
  - slurm_cluster=<CLUSTER_NAME>
  - slurm_qos=<QOS_LEVEL>      # e.g., 'low' if jobs are held in queue for long
  - runtime=<RUNTIME_MINUTES>  # e.g., 60
  - mem_mb=<MEMORY_MB>         # e.g., 4000

### Env modules ###
# use-envmodules: false 

### Conda ###
use-conda: true
conda-frontend: mamba   

### Resource scopes ###
set-resource-scopes:
  cores: local 

# Reusable Slurm Blocks (anchors)
# Standard partition/account/cluster used by most rules
_slurm_std: &slurm_std
  slurm_partition: <PARTITION_NAME>
  slurm_account: <ACCOUNT_NAME_standard> # e.g., standard, large memory 
  slurm_cluster: <CLUSTER_NAME>

# Large memory partition/account/cluster used by some rules
_slurm_large: &slurm_large
  slurm_partition: <PARTITION_NAME>
  slurm_account: <ACCOUNT_NAME_large> # e.g., standard, large memory 
  slurm_cluster: <CLUSTER_NAME>

## Per rule resources
set-resources:
  fastp_pe:
    <<: *slurm_std
    mem_mb: 4000
    runtime: 40

  kraken2:
    <<: *slurm_large
    mem_mb: 600000
    runtime: 30
```

#### 3. Configuration

The pipeline requires the following configuration files: `config.yaml`, `.env`, and `samples.txt`.

##### 3.1. config/config.yaml

The `config.yaml` file must be located in the `config` directory, which resides in the main Snakemake working directory. This file specifies crucial settings, including:

- Path to the `samples.txt`
- Input and output directories
- File paths to required databases
- Threads for each rule
- Parameters for software see the [Parameters](#parameters) section.

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
 TMPDIR = path/to/temp/on/cluster
 RGI_CARD = path/to/card.json and card_reference.fasta
```

##### 3.3. Sample list

`samplesheet.csv` Has the following column names: "sample","fastq_1","fastq_2". For the column 'sample" use the sampleID for the read pair, and for "fastq_1","fastq_2" have the names of the read1 and read2 files as they appear in the raw fastq files directory. The file location of the `samplesheet.csv` must be`config/samplesheet.csv`.

**Example `samplesheet.csv`:**
sample,fastq_1,fastq_2
test_LLC82Nov10GR,test_LLC82Nov10GR_r1.fastq.gz,test_LLC82Nov10GR_r2.fastq.gz
test_LLC82Sep06GR,test_LLC82Sep06GR_r1.fastq.gz,test_LLC82Sep06GR_r2.fastq.gz

##### 3.4. Scripts called in rules

The scripts called in the Snakemake pipeline are located in workflow/scripts.

- Module [taxonomy.smk](#module-taxonomysmk) uses the `extract_bracken_columns.py` script in `rule combine_bracken_outputs`.

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

- Defines resources for the job scheduler
- Activates the Snakemake environment
- Submits and manages jobs using the Snakemake `--profile` configuration `(profiles/slurm/)`.
- Contains any additional Snakemake arguments (e.g.., `--unlock`, `--dry-run`, `--rerun-incomplete`)
- For a snakemake report with runtime and software versions use --report path/to/metatranscriptomics_report.html after the pipeline has completed

##### 4.3. Submit launcher to SLURM

- Submit to SLURM compute node with bash terminal with `sbatch name_of_your_script.sh`
- Example below needs to be edited with the headers for the HPC you are using.

```bash
#!/bin/bash
#SBATCH --job-name=run_snakemake.sh
#SBATCH --output=run_snakemake_%j.out 
#SBATCH --error=run_snakemake_%j.err 
#SBATCH --cluster=<CLUSTER_NAME>
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --mem=<MEMORY_MB>         # e.g., 2000
#SBATCH --time=<HH:MM:SS>         # Must be long enough for completion of workflow 

source path/to/source/conda/common/miniforge/miniforge3/etc/profile.d/conda.sh

conda activate snakemake_env
export PATH="$PWD/bin:$PATH"

  snakemake \
    --profile absolute/path/to/profiles/slurm \
    --configfile absolute/path/to/config/config.yaml \
    --conda-prefix absolute/path/to/common/conda/metatranscriptomics-snakemake-conda \
    --printshellcmds \
    --keep-going 
```

### Notes

- The `profile/slurm/config.yaml` has been configured for our SLURM cluster. This will need to be configured for the cluster you are using.
- temp folder is set to `path/to/scratch/${USER}/tmpdir`
- A Snakemake report can be generated from the head node with `snakemake --report path/to/report/report_name.html`

#### Warnings

- The conda environments will not be created if the conda configuration is `conda config --set channel_priority strict`.
- Set conda to `conda config --set channel_priority flexible` or use libmamba.
- The `.env` file can overwrite the `config/config.yaml` file

#### Current issues

None.

#### Resource usage

- Kraken2: Large compute node with 600 GB. With 16 CUPs wall time was 7m 56s. With 2 CPUs wall time was 19m 13s.
- Generate Snakemake report to track walltime

---

## Output

**All output file paths are set in the `config/config.yaml` file and need to be edited prior to running the pipeline.**

The following table includes the key outputs of the metatranscriptomics pipeline. The [Snakemake rules](#snakemake-rules) section provides greater detail on all file outputs.

| Output Type                  | Description                                                                                                    | Filename                                                                                                                                                                                                                |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Processed sample reads       | Processed reads with Host reads and rRNA removed.                                                              | sample_rRNAdep_R1.fastq.gz / sample_rRNAdep_R2.fastq.gz                                                                                                                                                                 |
| Assembled transcripts        | Individual sample assemblies                                                                                   | sample.fasta                                                                                                                                                                                                            |
| Transcripts from Co-assembly | Co-assembly of all samples                                                                                     | final.contigs.fa                                                                                                                                                                                                        |
| Report                       | Kraken taxonomy summary for each sample                                                                        | sample.report.txt                                                                                                                                                                                                       |
| Report                       | Bracken report for the raw, and relative abundance at each taxonomic level                                     | Bracken_species_raw_abundance.csv, Bracken_species_relative_abundance.csv,Bracken_genus_raw_abundance.csv, Bracken_genus_relative_abundance.csv, Bracken_phylum_raw_abundance.csv, Bracken_genus_relative_abundance.csv |
| Report                       | Antimicrobial resistance gene profiling using RGI and the CARD.                                                | sample_paired.allele_mapping_data.txt, sample_paired.artifacts_mapping_stats.txt, sample_paired.gene_mapping_data.txt, sample_paired.overall_mapping_stats.txt, sample_paired.reference_mapping_stats.txt               |
| Report                       | rnaQUAST quality control report for individual sample assemblies using the BUSCO bacteria and archaea lineages | Reports are found in`sample_bacteria/` and `sample_archaea/` directories which contain the short_report files with .pdf, .tsv, and .txt extensions                                                                      |
| Report                       | Alignment statistics of the sample reads to the co-assembly                                                    | sample.flagstat.txt                                                                                                                                                                                                     |
| Report                       | Per-base sequencing depth across the co-assembly                                                               | sample.coverage.txt.gz                                                                                                                                                                                                  |
| Report                       | Sequence level mapping statistics with the sample contig name                                                  | sample.idxstats.txt.gz                                                                                                                                                                                                  |
| Annotation files             | Annotation tables for gene prediction of the co-assembly with protein sequences and nucleotide sequences       | coassembly.faa, coassembly.fna, coassembly.gff, coassembly.saf                                                                                                                                                          |
| Report                       | Feature count table for each sample                                                                            | sample_counts.txt                                                                                                                                                                                                       |

---
