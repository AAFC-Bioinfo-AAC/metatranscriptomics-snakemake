<!-- omit in toc -->
# METATRANSCRIPTOMICS SNAKEMAKE PIPELINE

[![FR](https://img.shields.io/badge/lang-FR-yellow.svg)](README_FR.md)
[![EN](https://img.shields.io/badge/lang-EN-blue.svg)](README.md)

---

<!-- omit in toc -->
## Table of Contents

- [About](#about)
- [Documentation](#documentation)
- [Acknowledgements](#acknowledgements)
- [Security](#security)
- [License](#license)

---

## About

The **Metatranscriptomics Snakemake** Pipeline uses paired-end FASTQ files from Illumina shotgun metatranscriptomic sequencing as input. The pipeline can be broken down into four main stages: sample read processing, sort read analysis, individual sample assembly, and co-assembly. **Sample processing** consists of *fastp*, *Bowtie2*, and *SortMeRNA* to perform quality filtering, remove host and PhiX contamination, and ribosomal (r)RNA depletion. These cleaned reads are used for the **sort read analysis** consisting of taxonomic classification with *Kraken2* using *GTDB* and antimicrobial gene profiling with *RGI* using *CARD*. **Individual samples are assembled** into presumptive messenger (m)RNA transcripts using *RNA SPAdes*. Assembly quality is evaluated with *rnaQUAST*. The **co-assembly stage** prepares the data for gene expression analysis. All cleaned reads are co-assembled with *MEGAHIT*, and the resulting co-assembly is indexed with *Bowtie2*. The cleaned sample reads are then mapped back to the co-assembly, and *SAMtools* is used generate assembly statistics, mapping summaries, and sequencing depth across the co-assembly. With *Prodigal* the protein and nucleotide coding regions of the co-assembly are predicted. *FeatureCounts* quantifies the predicted coding regions and generates a table for gene expression analysis. If metagenomic sequencing was done for these samples then the trimmed and host/PhiX removed metagenomic reads should be used for the co-assembly step.

Planned future enhancements to this pipeline include integrating *CoverM* to map metatranscriptomic reads onto the assembled metagenomes, as well as adding a *CAZyme analysis module* for functional annotation of carbohydrate-active enzymes.

---

## Documentation

For technical details, including installation and usage instructions, please see the [**`User Guide`**](./docs/user-guide.md).

---

## Acknowledgements

- **Credits**: This project was developed at the Lacombe Research and Development Centre, Agriculture & Agri-Food Canada by **Katherine James-Gzyl** and assisted by **Devin Holman** and **Arun Kommadath**.

- **Citation**: To cite this project, click the **`Cite this repository`** button on the right-hand sidebar

- **Contributing**: Contributions are welcome! Please review the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md) and ensure you adhere to our [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) to foster a respectful and inclusive environment.

- **References**: For a list of key resources used here, see [REFERENCES.md](REFERENCES.md)

---

## Security  

⚠️ Do not post any security issues on the public repository! Please report them as described in [SECURITY.md](SECURITY.md)

---

## License

See the [LICENSE](LICENSE) file for details. Visit [LicenseHub](https://licensehub.org) or [tl;drLegal](https://www.tldrlegal.com/) to view a plain-language summary of this license.

**Copyright ©** His Majesty the King in Right of Canada, as represented by the Minister of Agriculture and Agri-Food, 2025.

---
