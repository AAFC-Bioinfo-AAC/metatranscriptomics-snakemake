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

The **Metatranscriptomics Snakemake Pipeline** is a modular workflow designed to process, assemble, and analyze **Illumina paired-end shotgun metatranscriptomic data**. It coordinates analysis from raw read processing to gene-level quantification and CAZyme annotation using bioinformatics tools integrated through **Snakemake**. The workflow produces transcript assemblies, taxonomic profiles of RNA reads, antimicrobial resistance-associated transcript profiles, CAZyme annotations, and raw gene-count tables and CAZy-family count matrices for downstream analyses.

The pipeline consists of **five main stages**:

- **Sample read processing** — Quality filtering, removal of host and PhiX reads, and computational rRNA filtering using *fastp*, *Bowtie2*, and *SortMeRNA*.
- **Short-read analysis** — Taxonomic classification of filtered RNA reads with *Kraken2* using a GTDB-based database, and read-based profiling of antimicrobial resistance-associated transcripts with *RGI* using CARD.
- **Individual sample assembly** — Transcript assembly with *rnaSPAdes* and assembly evaluation with *rnaQUAST*.
- **Co-assembly and gene quantification** — Co-assembly of filtered reads across samples with *MEGAHIT*, prediction of prokaryotic coding sequences with *Prodigal*, mapping of filtered metatranscriptomic reads with *Bowtie2*, mapping and sequencing-depth statistics with *SAMtools*, and gene-level counting with *FeatureCounts*.
- **CAZyme annotation and transcript quantification** — Annotation of predicted proteins from the shared reference using *dbCAN*, followed by generation of raw RNA count matrices for CAZyme genes and CAZy families from the FeatureCounts results.

  💡 When matched metagenomic data are available, a suitable assembly generated from cleaned metagenomic reads can provide a common reference for mapping and quantifying metatranscriptomic reads and annotating CAZymes. Setting reference_assembly to this assembly’s FASTA path bypasses RNA co-assembly for the shared reference.

  Taxonomic profiles derived from RNA reads reflect the representation of classified transcripts rather than directly measuring microbial cell abundance. Gene-count tables require appropriate downstream processing and statistical analysis to     assess differential expression. Family counts can overlap because genes assigned to multiple CAZy families contribute their counts to each assigned family.

Some **future enhancements** planned for this workflow include:

- Integration of *CoverM* for mapping metatranscriptomic reads to metagenomic references and summarizing coverage.
- KEGG annotation module to assign KEGG Orthology (KO) identifiers to predicted proteins and support functional analyses using KEGG pathways and modules.

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
