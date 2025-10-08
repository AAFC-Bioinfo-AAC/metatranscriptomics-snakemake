<!-- omit in toc -->
# METATRANSCRIPTOMICS SNAKEMAKE PIPELINE

[![FR](https://img.shields.io/badge/lang-FR-yellow.svg)](README_FR.md)
[![EN](https://img.shields.io/badge/lang-EN-blue.svg)](README.md)

---

<!-- omit in toc -->
## Table des matières

- [À propos](#à-propos)
- [Documentation](#documentation)
- [Remerciements](#remerciements)
- [Sécurité](#sécurité)
- [Licence](#licence)

---

## À propos

Le pipeline **Metatranscriptomics Snakemake** utilise comme entrée des fichiers FASTQ appariés issus d’un séquençage métatranscriptomique shotgun Illumina. Le pipeline se divise en quatre grandes étapes : le traitement des lectures d’échantillons, l’analyse des lectures triées, l’assemblage individuel des échantillons et la co-assemblage. Le **traitement des échantillons** comprend l’utilisation de *fastp*, *Bowtie2* et *SortMeRNA* pour effectuer le filtrage de qualité, éliminer la contamination par l’hôte et PhiX, ainsi que pour la déplétion de l’ARN ribosomique (ARNr). Les **lectures nettoyées** servent ensuite à l’**analyse des lectures triées**, qui comprend la classification taxonomique avec *Kraken2* utilisant la base de données GTDB et le profilage des gènes antimicrobiens avec *RGI* utilisant la base *CARD*. Les **échantillons individuels** sont assemblés en transcrits d’ARN messager (ARNm) présumés à l’aide de *RNA SPAdes*. La qualité de l’assemblage est évaluée avec *rnaQUAST*. L’étape de **co-assemblage** prépare les données pour l’analyse de l’expression génique. Toutes les lectures nettoyées sont co-assemblées avec *MEGAHIT*, et le co-assemblage obtenu est indexé avec *Bowtie2*. Les lectures nettoyées des échantillons sont ensuite réalignées sur le co-assemblage, et *SAMtools* est utilisé pour générer des statistiques d’assemblage, des résumés d’alignement et la profondeur de séquençage à travers le co-assemblage. Avec *Prodigal*, les régions codantes protéiques et nucléotidiques du co-assemblage sont prédites. *FeatureCounts* quantifie ces régions codantes prédites et génère un tableau pour l’analyse de l’expression génique. Si un séquençage métagénomique a également été effectué pour ces échantillons, les lectures métagénomiques nettoyées (après retrait des séquences PhiX et de l’hôte) doivent être utilisées pour l’étape de co-assemblage.

Les améliorations futures prévues pour ce pipeline incluent l’intégration de *CoverM* afin de mapper les lectures métatranscriptomiques sur les métagénomes assemblés, ainsi que l’ajout d’un *module d’analyse des CAZymes* pour l’annotation fonctionnelle des enzymes actives sur les glucides.

---

## Documentation

Pour les détails techniques, y compris les instructions d’installation et d’utilisation, veuillez consulter le [Guide de l’utilisateur](/docs/user-guide.md).

---

## Remerciements

- **Crédits** : Ce projet a été réalisé au Centre de recherche et de développement de Lacombe, Agriculture et Agroalimentaire Canada, par **Katherine James-Gzyl**, avec l’aide de **Devin Holman** et d’**Arun Kommadath**.

- **Citation** : Pour citer ce projet, cliquez sur le bouton **`Cite this repository`** dans la barre latérale de droite.

- **Contribution** : Les contributions sont les bienvenues ! Veuillez consulter les lignes directrices dans [CONTRIBUTING.md](CONTRIBUTING.md) et vous assurer de respecter notre [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) afin de favoriser un environnement respectueux et inclusif.

- **Références** : Pour une liste des principales ressources utilisées ici, voir [REFERENCES.md](REFERENCES.md)

---

## Sécurité

⚠️ Ne publiez aucun problème de sécurité sur le répertoire public ! Veuillez les signaler comme décrit dans [SECURITY.md](SECURITY.md).

---

## Licence

Voir le fichier [LICENSE](LICENSE) pour plus de détails. Visitez [LicenseHub](https://licensehub.org/fr) ou [tl;drLegal](https://www.tldrlegal.com/) pour consulter un résumé en langage clair de cette licence.

**Droit d’auteur ©** Sa Majesté le Roi du chef du Canada, représenté par le ministre de l’Agriculture et de l’Agroalimentaire, 2025.

---
