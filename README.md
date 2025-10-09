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

The **Metatranscriptomics Snakemake Pipeline** is a reproducible, modular workflow designed to process, assemble, and analyze **Illumina paired-end shotgun metatranscriptomic data**. It automates the full analysis—from raw read processing to gene expression quantification—using widely adopted bioinformatics tools integrated through **Snakemake**. The pipeline produces high-quality assemblies, taxonomic and antimicrobial gene profiles, and quantitative gene expression tables suitable for downstream statistical and functional analyses.

The pipeline consists of **four main stages**:

- **Sample Read Processing** — Quality filtering, removal of host and PhiX contamination, and ribosomal RNA (rRNA) depletion using *fastp*, *Bowtie2*, and *SortMeRNA*.
- **Short-Read Analysis** — Taxonomic classification with *Kraken2* (using *GTDB*) and antimicrobial resistance profiling with *RGI* (using *CARD*).
- **Individual Sample Assembly** — Transcript assembly with *RNA SPAdes* and quality assessment via *rnaQUAST*.
- **Co-assembly and Expression Quantification** — Global co-assembly with *MEGAHIT*, followed by read mapping (*Bowtie2*), coverage assessment (*SAMtools*), gene prediction (*Prodigal*), and feature quantification (*FeatureCounts*).

  💡 *If metagenomic sequencing data are available for the same samples, trimmed and host/PhiX-filtered metagenomic reads should be used in the co-assembly stage.*

Some **future enhancements** planned for this pipeline include:

- Integration of *CoverM* for mapping metatranscriptomic reads to assembled metagenomes.
- Addition of a *CAZyme analysis module* for functional annotation of carbohydrate-active enzymes.

---

## Documentation

For technical details, including installation and usage instructions, please see the [**`User Guide`**](./docs/user-guide.md).

---

## Acknowledgements

- **Credits**: This project was developed at the *Lacombe Research and Development Centre, Agriculture & Agri-Food Canada (AAFC)* by **Katherine James-Gzyl** and assisted by **Devin Holman** and **Arun Kommadath**.

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
