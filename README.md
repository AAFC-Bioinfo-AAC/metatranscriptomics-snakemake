# Metatranscriptomics Snakemake Pipeline

[![Repository Template](https://img.shields.io/badge/repository-template-blue)](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-template-repository)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## ABOUT

*Summary of the project:*

The Metatranscriptomics Snakemake Pipeline uses paired-end FASTQ files from Illumina shotgun metatranscriptomic sequencing as input. The first part of the pipeline processes the reads using Fastp, Bowtie2, and SortMeRNA to perform quality filtering, deplete host and PhiX reads, and removes ribosomal (r)RNA. These cleaned reads are then used for taxonomic classification with Kraken2 and GTDB, antimicrobial gene profiling with RGI and CARD, and transcriptome assembly into presumptive messenger (m)RNA transcripts using RNA SPAdes. To asses the quality of transcripts rnaQUAST and mapping the cleaned reads back to assembly are used. The SAM files from mapping the reads back to the assembly can be used in further expression studies. 

- *Tools that still need to be added to the pipeline are CoverM, and CAZyme analysis*
- *The different ways the code can be configured or customized for specific use cases.*
- *Could include a brief mention of any unique features or benefits of the project.*

---

## TABLE OF CONTENTS


| **Section**                                                                | **Description**                                                                                                   |
| ---------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| [ABOUT](#about)                                                            | *A summary of the project, may include its origin, purpose, and functionality, along with configuration options.* |
| [OVERVIEW](#overview)                                                      | *A summary of the project's processes, supported by a visual representation (e.g., a pipeline diagram).*          |
| [DATA](#data)                                                              | *Details of the data files used in the project.*                                                                  |
| [PARAMETERS](#parameters)                                                  | *A table describing configurable parameters, their expected values, and their impact on the output.*              |
| [USAGE](#usage)                                                            | *Detailed guidance on how to use the project, including pre-requisites, instructions, and optional notes.*        |
| &nbsp;&nbsp;&nbsp;&nbsp;[Pre-requisites](#pre-requisites)                  | *Dependencies and hardware/software requirements.*                                                                |
| &nbsp;&nbsp;&nbsp;&nbsp;[Instructions](#instructions)                      | *Step-by-step directions for running the code, including examples and links to related resources.*                |
| &nbsp;&nbsp;&nbsp;&nbsp;[Notes](#notes)                                    | *Additional optional details, tips, or alternative methods.*                                                      |
| [OUTPUT](#output)                                                          | *Details of the output files generated, which may include formats, locations, and naming conventions.*            |
| [CREDITS](#credits)                                                        | *Acknowledgment of contributors, teams, and organizations that supported the project.*                            |
| [CONTRIBUTION](#contribution)                                              | *Guidelines for contributing to the repository, with a link to the `CONTRIBUTING.md` file.*                       |
| [COPYRIGHT](#copyright)                                                    | *Ownership details*                                                                                               |
| [LICENSE](#license)                                                        | *Information about the license, including a link to the `LICENSE` file.*                                          |
| [PUBLICATIONS & ADDITIONAL RESOURCES](#publications--additional-resources) | *Links to publications, articles, or other resources related to the project.*                                     |
| [CITATION](#citation)                                                      | *Instructions for citing the project, with references to the `CITATION.cff` and `CITATIONS.md` files.*            |

---

## OVERVIEW

*Provide a summary of the steps or processes the code performs. Include:*

- **`rule fastp_pe`** reads the sample names from the `samples.txt` file and processes paired-end FASTQ files with the naming conventions `*_R1.fastq.gz` / `*_r1.fastq.gz` and `*_R2.fastq.gz` / `*_r2.fastq.gz`. This step performs adapter trimming, quality trimming, and filtering. The output files are the trimmed paired reads (`*_r1.fastq.gz`/`*_r2.fastq.gz`), unpaired reads (`*_u1.fastq.gz`/`*_u2.fastq.gz`), and QC reports. Default parameters are used. The unpaired reads can be removed after confirming that the majority of reads have a pair that passed QC.

- **`rule bowtie2_align`** uses the trimmed paired-end files from `rule fastp_pe` and aligns them to the specified index. The resulting file is a reference-aligned BAM file. Default parameters are used.

- **`rule extract_unmapped_fastq`** takes the sorted BAM file as input and extracts the reads that did not align into paired-end FASTQ files depleted of host and PhiX reads (`*_trimmed_clean_R1.fastq.gz`/`*_trimmed_clean_R2.fastq.gz`). Default parameters are used.

- **`rule sortmerna`** aligns the host-depleted reads to an rRNA database and outputs the rRNA-depleted reads (`*_rRNAdep_R1.fastq.gz`/`*_rRNAdep_R2.fastq.gz`). These rRNA-depleted reads will be used for downstream analysis. The database used for testing the pipeline was `smr_v4.3_default_db.fasta`, available [here](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3). Default parameters were used.

- **`rule rna_spades`**  
  > **Note to self:** The order of this rule should be changed in the Snakemake workflow. Keep for now while testing.

   > **Note to self:** In the shell block copy the spades.log. It contains some assembly stat that are not printed to the current log file. 

  The rRNA-depleted reads are assembled into presumptive mRNA transcripts using the `--rna` flag and default parameters. The transcripts are output to a `*.fasta` file. 

- **`rule rnaquast_busco`** uses the transcripts from RNA SPAses to print the number of transcripts, transcripts over 500 bp, transcripts over 1000 bp and the BUSCO completeness. The software is not intended for metatranscriptomics. Use caution when interpreting the results. For instance the BUSCO completeness cannot be interpreted as the percentage of assembly quality but instead it is a representation of the core functions from the bacteria_odb12 and archaea_odb12 lineages.  
  

- **`rule bowtie_build_transcript`, `rule bowtie2_map_transcripts` and `assembly_stats_depth`** An index of the transcripts for each sample is built, the rRNA depleted reads are mapped to the transcript file and sorted BAM files are outputted. These BAM files can be used for downstream expression studies. Stats from the SAM files for each sample are printed with samtools flagstat.
  > **Note to self:** determine which files should be marked temporary.  

-**`rule Salmon/kallisto/RSEM` Not currently in pipeline but can be added. 


- **`rule kraken2`** assigns taxonomy to the rRNA-depleted reads using a Kraken2-formatted GTDB. A confidence threshold of 0.5 is used and all other parameters are defaults. The output files are `*.report` and `*.kraken`.

- **`rule bracken`** uses the report file from kraken to output a report at the species, genus and phylum level for each sample. These intermediate files are fed into `rule combine_bracken_outputs`.
   > **Note to self:** This rule is also making `.report_bracken_species.txt` at each level in the `06_kraken` directory. At some point see if we can either place these into a direcorory called `reports` or have them cleaned up in the shell block.

- **`rule combine_bracken_outputs`** combines the reports for all the samples   
  > **Note to self:**  currently using a sub set sample list that will need to be removed after testing. The `sample=SAMPLES_SUBSET` will need to be changed to `SAMPLES`. The subset lines will need to be removed from the `config/config.yaml` file and the top of the `Snakemake`.

- **`rule bracken_extract`** used a python script `scripts/extract_bracken_columns.py` to generate tables for the raw and relative abundance for each taxonomic level used in `rule bracken`. The resulting outputs are `Bracken_[species/genus/phylum]_relative_abundance.csv` and Bracken_[species/genus/phylum]_raw_abundance.csv`.

- **`rule rgi_reload_database`**  
  > **Note to self:** Include a function that allows the rule to not fail if a local CARD database is used.

  Loads the CARD database from a common folder to the working directory. After the step is completed, there should be a `localDB` folder in the main Snakemake directory and a `rgi_reload_db.done` file in the logs directory to prevent the rule from re-running every time the pipeline is called. 

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

- **`rule coverm`**
   > **Note to self:** Add in the option of running CoverM. This should not be part of the main pipeline but an option if MAGs from metagenomic sequencing of the same samples are avalible.

- **`rules Cazymes`**
     > **Note to self:** Do we want to include this in the pipeline or use transcripts in exisiting Bash pipeline made by Arun. 
 
- *A diagram or workflow visual that illustrates the main steps or processes.*

**Example**:

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
    T --> V((Read Mapping Stats))
    T --> X((Sorted BAM files))
```
---

## DATA

The raw input data must be in the form of paired-end FASTQ files generated from metatranscriptomics experiments.

- Each sample should include both forward (R1) and reverse (R2) read files.
- Both uppercase (`R1`/`R2`) and lowercase (`r1`/`r2`) naming formats are accepted (e.g., `sample_R1.fastq.gz`, `sample_r2.fastq.gz`).
- The files must be organized in a directory named `01_raw`.
- The path to the `01_raw` directory must be specified in the `config.yaml` file.

**Example:**

- **Dataset 1 Filename**: Sequencing reads (FASTQ) from beef cattle rumen samples are provided for three samples: `LLC42Nov10C`, `LLC42Sep06CR`, and `LLC82Sep06GR`. There is also a subsampled test file, `test_LLC82Sep06GR`, which can be used for all steps dealing with unassembled data but will fail for the RNA SPAdes step.


*For large external datasets, provide links or instructions to download or utilize them.*

**Example:**

To download the raw data:

```
curl -O https://example.com/path/to/dataset1.tar.gz
```

---

## PARAMETERS


| Parameter          | Value                                                                                               |
| -------------------- | ----------------------------------------------------------------------------------------------------- |
| *parameter_name_1* | *Description of what the parameter does and the expected value (e.g., integer, string, file path).* |
| *parameter_name_2* | *Description of what the parameter does and the expected value (e.g., boolean, list).*              |

---

## USAGE

### Warnings
- **Enviroment varable function**
  The `.env` file can overwirte the `config/config.ymal` file

- **Current Issues**
  * The TMPDIR varable set in the .env file is not working. For a temp fix the TMPDIR has been set in the main Snakemake workflow

  * temp folder is set to `/gpfs/fs7/aafc/scratch/$USER/tmpdir` for running on the GPSC.

  * `$USER` is changed to actual user name. This was done for testing. I wasn't sure if the `$USER` was being expanded correctly. 
### Pre-requisites

#### Software
- Snakemake version 9.6.0

#### Additional Data Files and Databases for Software

- **Bowtie2**  
  Bowtie2 uses an index of reference sequences to align reads. This index must be created before running the pipeline. The index files (with the `.bt2` extension) must be located in the `index` directory. Make sure to update the prefix of these files in the `config.yaml` file.

- **SortMeRNA**  
  SortMeRNA requires a ribosomal (r)RNA database in the `rRNA_DB` directory. Update the `config.yaml` file with the filename of the database used. You can download the database from [SortMeRNA releases](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3). The file `smr_v4.3_default_db.fasta` was used for pipeline testing.

- **Kraken2**  
  Kraken2 requires a Kraken2-formatted GTDB database. The GTDB release tested with this pipeline was 220.

- **RGI BWT/CARD**  
  - * PLACEHOLDER: reminder to fix RGI CARD steps so that it does not fail if a user loaded database is used. 

  RGI BWT requires the CARD (Comprehensive Antibiotic Resistance Database) database. The version tested in this pipeline was 4.0.1. The database can be located on a common drive or in your working directory.  
  Instructions for installing the CARD database are available [here](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst).  
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

- The conda enviroments used in the pipeline need to be generated before running with `snakemake --use-conda --conda-create-envs-only`
- *Specific programming languages, libraries, or frameworks (e.g., Python 3.9, NumPy).*
- *Installation instructions for dependencies (e.g., pip install, conda environments).*
- *Hardware requirements, if any (e.g., CPU/GPU specifications, memory, specs used when running with SLURM).*

### Instructions

#### 1. Installation
Clone the repository into the directory where you want to run the metatranscriptomics Snakemake pipeline.  
**Note:** This location must be on an HPCC (High Performance Computing Cluster) with access to a high-memory node (at least 600 GB RAM) and sufficient storage for all metatranscriptomics analyses.

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
For best practice, use database paths that are in common locations to all users on the HPCC.

##### 2.2. Environment file
This file must contain paths to the TMPDIR file and the RGI databse. Follow these instructions:
- In the main Snakemake direcory (where you are running Snakemake from)
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

##### 3.1.Conda Environments
Snakemake can automatically create and load Conda environments for each rule in your workflow. Check to see that you have the following configuration files in the `envs` directory:
- `bedtools.yaml`
- `bowtie2.yaml`
- `kraken2.yaml`
- `rgi.yaml`
- `RNAspades.yaml`
- `sortmerna.yaml`

Load the required conda enviroments for the pipeline with:
```bash
snakemake --conda-create-envs-only
``` 
  

*DETAILED Step-by-step guide to running the code. Can include:*

- *Command-line examples or scripts to execute.*
- *Screenshots, images, or videos illustrating usage.*
- *Links to detailed documentation or tutorials.*
- *Diagrams showing data flow or system behavior.*

### Notes
- Kraken2: Large compute node with 600 GB. With 16 CUPs wall time was 7m 56s. With 2 CPUs wall time was 19m 13s. 
*IF APPLICABLE: Any information, such as tips, warnings, or alternative ways to run the code.*

*OTHERWISE: Write N/A*

---

## OUTPUT

*Describe the expected outputs of the code. Include:*

- *File types (e.g., `.csv`, `.txt`, `.bam`).*
- *Location of the files.*
- *File naming conventions.*
- *Examples of output files or links to them, if applicable.*

---

## CREDITS

**Example:**
"This repository was written by [Your Name/Team Name]."
"We thank the following people and teams for their assistance in the development of this project:"

- [Contributor 1]
- [Contributor 2]
- [Acknowledged Organizations or Teams]

---

## CONTRIBUTION

If you would like to contribute to this project, please consult [CONTRIBUTING.md](CONTRIBUTING.md)

---

## COPYRIGHT

Government of Canada, Agriculture & Agri-Food Canada

---

## LICENSE

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## PUBLICATIONS & ADDITIONAL RESOURCES

*IF APPLICABLE: Include any publications, articles, or additional resources that are related to the project.*

- *Links to related papers or articles.*
- *References for bioinformatics tools or methods used in the code.*

#### ABCC_RCBA_Guide

Guidelines (under development) for additional context and supplementary materials that align with this project.

---

## CITATION

If you use this repository for your analysis, please cite it using the [CITATION.cff](CITATION.cff) file. An extensive list of references for the tools used can be found in the [CITATIONS.md](CITATIONS.md) file.
