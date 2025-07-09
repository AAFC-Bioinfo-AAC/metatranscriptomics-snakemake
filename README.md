# Metatranscriptomics Snakemake Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## About

The Metatranscriptomics Snakemake Pipeline uses paired-end FASTQ files from Illumina shotgun metatranscriptomic sequencing as input. The first part of the pipeline processes the reads using fastp, Bowtie2, and SortMeRNA to perform quality filtering, deplete host and PhiX reads, and removes ribosomal (r)RNA. These cleaned reads are then used for taxonomic classification with Kraken2 and GTDB, antimicrobial gene profiling with RGI and CARD, and transcriptome assembly into presumptive messenger (m)RNA transcripts using RNA SPAdes. To asses the quality of transcripts rnaQUAST and mapping the cleaned reads back to assembly are used. The SAM files from mapping the reads back to the assembly can be used in further expression studies.

- *Tools that still need to be added to the pipeline are CoverM, and CAZyme analysis*
- *The different ways the code can be configured or customized for specific use cases.*
- *Could include a brief mention of any unique features or benefits of the project.*

---

## Table of Contents

- [Metatranscriptomics Snakemake Pipeline](#metatranscriptomics-snakemake-pipeline)
  - [About](#about)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
    - [Workflow diagram](#workflow-diagram)
    - [Snakemake rules](#snakemake-rules)
  - [Data](#data)
  - [Parameters](#parameters)
  - [Usage](#usage)
    - [Pre-requisites](#pre-requisites)
      - [Software](#software)
      - [Databases](#databases)
    - [Setup Instructions](#setup-instructions)
      - [1. Installation](#1-installation)
      - [2. Configuration](#2-configuration)
        - [2.1. config.yaml](#21-configyaml)
        - [2.2. Environment file](#22-environment-file)
        - [2.3. Sample list](#23-sample-list)
      - [3. Running the pipeline](#3-running-the-pipeline)
        - [3.1.Conda environments](#31conda-environments)
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
        B --> C[Trimmed Reads]
        B --> L((Fastp QC Report))
        C -->|Host/PhiX removal| D{Bowtie2}
        D --> E[Non-host, Non-PhiX Reads]
    end

    subgraph DEPLETION [rRNA Depletion]
        E --> F{SortMeRNA}
        F --> G[rRNA-depleted Reads]
    end

    subgraph ASSEMBLY [Assembly]
        G --> H{rnaSPAdes}
        H --> I((Assembled Transcripts))
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

    %% Optional transition from assembled transcripts back to Bowtie2 mapping (read support check)
    I --> T{Bowtie2}
    T --> V((Read Mapping, 
            Coverage, and
            Alignment Stats))
    T --> X((Sorted BAM files))
```

### Snakemake rules

- **`rule fastp_pe`** reads the sample names from the `samples.txt` file and processes paired-end FASTQ files with the naming conventions `*_R1.fastq.gz` / `*_r1.fastq.gz` and `*_R2.fastq.gz` / `*_r2.fastq.gz`. This step performs adapter trimming, quality trimming, and filtering. The output files are the trimmed paired reads (`*_r1.fastq.gz`/`*_r2.fastq.gz`), unpaired reads (`*_u1.fastq.gz`/`*_u2.fastq.gz`), and QC reports. Default parameters are used. The unpaired reads can be removed after confirming that the majority of reads have a pair that passed QC.

- **`rule bowtie2_align`** uses the trimmed paired-end files from `rule fastp_pe` and aligns them to the specified index. The resulting file is a reference-aligned BAM file. Default parameters are used.
  > **Wall time:** Tests for single sample. When total cores was 60 (bowtie2 44, SAMtools view 4, SAMtools sort 12) time was 12m 18s. Reduced to total cores 24 (bowtie2 16, SAMtools view 4, SAMtools sort 8)

- **`rule extract_unmapped_fastq`** takes the sorted BAM file as input and extracts the reads that did not align into paired-end FASTQ files depleted of host and PhiX reads (`*_trimmed_clean_R1.fastq.gz`/`*_trimmed_clean_R2.fastq.gz`). Default parameters are used.

  > **Wall time:** Tests for single sample. When 60 cores were used and there was no splitting the wall time was 17m 57s. Changed to a 80:20 ratio of splitting between SAMtools and pigz. Total cores is 60 so SAMtools will get 48 and pigz 12. A temp directory was added and wall time logging has been improved.

- **`rule sortmerna`** aligns the host-depleted reads to an rRNA database and outputs the rRNA-depleted reads (`*_rRNAdep_R1.fastq.gz`/`*_rRNAdep_R2.fastq.gz`). These rRNA-depleted reads will be used for downstream analysis. The database used for testing the pipeline was `smr_v4.3_default_db.fasta`, available from the Reference RNA databases (database.tar.gz) file at [sortmerna release v4.3.3](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3). Default parameters were used.

  > **Wall time:** With 60 cores wall time was 4h 48m 56s for one sample. See what wall time is with a reduction to 48 cores and then 32.

- **`rule rna_spades`**

  > **Wall time:** One sample with 60 cores ran for 2h 15m 53s. The wall time logging and log file has been fixed. Try 48 and 32 cores to see if wall time is similar. Reduce cores if so.
  > **Note to self:** The order of this rule should be changed in the Snakemake workflow. Keep for now while testing.
   > **Note to self:** In the shell block copy the spades.log. It contains some assembly stat that are not printed to the current log file. SHOULD BE FIXED. CHECK AFTER NEXT RUN.

  The rRNA-depleted reads are assembled into presumptive mRNA transcripts using the `--rna` flag and default parameters. The transcripts are output to a `*.fasta` file.

- **`rule rnaquast_busco`** uses the transcripts from RNA SPAses to print the number of transcripts, transcripts over 500 bp, transcripts over 1000 bp and the BUSCO completeness. The software is not intended for metatranscriptomics. Use caution when interpreting the results. For instance the BUSCO completeness cannot be interpreted as the percentage of assembly quality but instead it is a representation of the core functions from the bacteria_odb12 and archaea_odb12 lineages.  
  
- **`megahit_coassembly`** here the rRNA-depleted reads are co-assembled with MEGAHIT. If Metagenomic sequencing was done the co-assembly of those reads would be a better choice. The Co-assembly is used as a index to produce sorted BAM files for each assembly. These sorted BAM files can then be used in featureCounts and downstream expression analysis.

  > **Wall time:** For the co-assembly of three samples was 52 min 28 sec with 60 cores. For large co-assemblies a large mem node will need to be used (will need to be tested at some point)

- **`index_coassembly`** use Bowtie2 to make an index that can be used to map the reads to the co-assembly.
  > **Wall time:** Wall time with 8 threads was 1m 5s.

- **`bowtie2_map_transcripts`** Maps the cleaned reads to the co-assembly resulting in a `.coassembly.sorted.bam` for each sample.

  > **Wall time:** For one sample using 40 cores was 9m 53s. Reduced cores to 16 and will check wall time.

- **`assembly_stats_depth`** produces a `flagstat.txt` summary of the reads that mapped back to the assembly, a `coverage.txt.gz` depth file with per-base coverage across the co-assembly for each sample, and a `idxstats.txt.gz` with read counts per transcript for each sample.

- **`rule prodigal_genes`** used to make a FASTA file of predicted protein sequences `coassembly.faa` and the predicted genes `coassembly.fna`, a feature formatted annotation file `coassembly.gff` and a simplified annotation formate file `coassembly.saf` that is used by feature counts.

  > **Temp file:** Go back and decided if this output should be designated temporary
  > **Wall time** was 8m with 1 core. prodigal does not support more than one core.

- **`rule featurecounts`** generates a table for each sample that includes the Geneid (unique identifier), the co-assembly contig name, the start and end positions of each gene on the contig, the strand orientation (+ or -), the gene length, and the number of reads mapped to each gene. Since all samples are mapped to the same co-assembly reference, the resulting tables can be combined for downstream analysis of gene expression across samples.

  > **Wall time** for one sample with 4 cores was 17s.

- **`rule kraken2`** assigns taxonomy to the rRNA-depleted reads using a Kraken2-formatted GTDB. A confidence threshold of 0.5 is used and all other parameters are defaults. The output files are `*.report` and `*.kraken`.
  > **Wall time** Large compute node with 600 GB. With 16 CUPs wall time was 7m 56s. With 2 CPUs wall time was 19m 13s.

- **`rule bracken`** uses the report file from kraken to output a report at the species, genus and phylum level for each sample. These intermediate files are fed into `rule combine_bracken_outputs`.

   > **Wall time** with 10 threads the wall time was 9s. The cores have been reduced to 2.
   > **Note to self:** This rule is also making `.report_bracken_species.txt` at each level in the `06_kraken` directory. At some point see if we can either place these into a directory called `reports` or have them cleaned up in the shell block.

- **`rule combine_bracken_outputs`** combines the reports for all the samples  

- **`rule bracken_extract`** used a python script `scripts/extract_bracken_columns.py` to generate tables for the raw and relative abundance for each taxonomic level used in `rule bracken`. The resulting outputs are `Bracken_[species/genus/phylum]_relative_abundance.csv` and Bracken_[species/genus/phylum]_raw_abundance.csv`.

- **`rule rgi_reload_database`** Loads the CARD database from a common folder to the working directory only if `localDB` has not been previously loaded. After the step is completed, there should be a `localDB` folder in the main Snakemake directory and a `rgi_reload_db.done` file in the logs directory to prevent the rule from re-running every time the pipeline is called.

- **`rule rgi_bwt`** performs antimicrobial resistance gene profiling on the rRNA-depleted reads using *k*-mer alignment (kma) and default parameters. Output files are:  
  - `*_paired.allele_mapping_data.json`
  - `*_paired.allele_mapping_data.txt`
  - `*_paired.artifacts_mapping_stats.txt`
  - `*_paired.gene_mapping_data.txt`
  - `*_paired.overall_mapping_stats.txt`
  - `*_paired.reference_mapping_stats.txt`
  - `*_paired.sorted.length_100.bam`
  - `*_paired.sorted.length_100.bam.bai`  
  Remove files that are not required after this step completes.

  > **Wall time** with 40 cores the job took 18m 7s. Try reducing to 20 cores. If time does not increase much further reduce cores.
- **`rule coverm`**
   > **Note to self:** Add in the option of running CoverM. This should not be part of the main pipeline but an option if MAGs from metagenomic sequencing of the same samples are available.

- **`rules Cazymes`**
     > **Note to self:** Do we want to include this in the pipeline or use transcripts in existing Bash pipeline made by Arun.

---

## Data

The raw input data must be in the form of paired-end FASTQ files generated from metatranscriptomics experiments.

- Each sample should include both forward (R1) and reverse (R2) read files.
- Both uppercase (`R1`/`R2`) and lowercase (`r1`/`r2`) naming formats are accepted (e.g., `sample_R1.fastq.gz`, `sample_r2.fastq.gz`).
- The files must be organized in a directory named `01_raw`.
- The path to the `01_raw` directory must be specified in the `config.yaml` file.

**Example:**

- **Dataset 1 Filename**: Sequencing reads (FASTQ) from beef cattle rumen samples are provided for three samples: `LLC42Nov10C`, `LLC42Sep06CR`, and `LLC82Sep06GR`. There is also a sub-sampled test file, `test_LLC82Sep06GR`, which can be used for all steps dealing with un-assembled data but will fail for the RNA SPAdes step.

---

## Parameters

| Parameter          | Value                                                                                               |
| -------------------- | ----------------------------------------------------------------------------------------------------- |
| *parameter_name_1* | *Description of what the parameter does and the expected value (e.g., integer, string, file path).* |
| *parameter_name_2* | *Description of what the parameter does and the expected value (e.g., boolean, list).*              |

---

## Usage

### Pre-requisites

#### Software

- Snakemake version 9.6.0

#### Databases

- **Bowtie2**  
  Bowtie2 uses an index of reference sequences to align reads. This index must be created before running the pipeline. The index files (with the `.bt2` extension) must be located in the `index` directory. Make sure to update the prefix of these files in the `config.yaml` file.

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

- The conda environments used in the pipeline need to be generated before running with `snakemake --use-conda --conda-create-envs-only`

### Setup Instructions

#### 1. Installation

Clone the repository into the directory where you want to run the metatranscriptomics Snakemake pipeline.  
**Note:** This location must be on an HPC (High Performance Computing) cluster with access to a high-memory node (at least 600 GB RAM) and sufficient storage for all metatranscriptomics analyses.

```bash
cd /path/to/metatranscriptomics
git clone <repository-url>
```

#### 2. Configuration

The pipeline requires the following configuration files: `config.yaml`, `.env`, and `sample.txt`.

##### 2.1. config.yaml

The `config.yaml` file must be located in the `config` directory, which resides in the main Snakemake working directory. This file specifies crucial settings, including:

- Path to the `sample.txt`  
- Input and output directories  
- File paths to required databases  

**Note:**  
You must edit `config.yaml` **before** running the pipeline to ensure all paths are correctly set.  
For best practice, use database paths that are in common locations to all users on the HPC.

##### 2.2. Environment file

This file must contain paths to the TMPDIR file and the RGI database. Follow these instructions:

- In the main Snakemake directory (where you are running Snakemake from)

```bash
touch .env
```

- Open the .env file and add

```bash
 TMPDIR = path/to/temp/on/cluster
 RGI_CARD = path/to/card.json and card_reference.fasta
 ```

##### 2.3. Sample list

 Sample names must be stored in a plain text file named `samples.txt`. Each sample name should appear on a separate line, with no additional formatting or headers.

**Example `samples.txt`:**
`LLC42Nov10CR`
`LLC42Sep06CR`
`LLC82Sep06GR`

#### 3. Running the pipeline

Complete steps **1.Installation** and **2.Configuration** and ensure database paths have been added to the 'config.yaml'. Required databases are described in the [Pre-requisites](#pre-requisites).

##### 3.1.Conda environments

Snakemake can automatically create and load Conda environments for each rule in your workflow. Check to see that you have the following configuration files in the `envs` directory:

- `bedtools.yaml`
- `bowtie2.yaml`
- `kraken2.yaml`
- `rgi.yaml`
- `RNAspades.yaml`
- `sortmerna.yaml`

Load the required conda environments for the pipeline with:

```bash
snakemake --conda-create-envs-only
```
  
### Notes

#### Warnings

- The `.env` file can overwrite the `config/config.yaml` file

#### Current issues

- The TMPDIR variable set in the .env file is not working. For a temp fix the TMPDIR has been set in the main Snakemake workflow

- temp folder is set to `/gpfs/fs7/aafc/scratch/$USER/tmpdir` for running on the GPSC.

- `$USER` is changed to actual user name. This was done for testing. I wasn't sure if the `$USER` was being expanded correctly.

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
