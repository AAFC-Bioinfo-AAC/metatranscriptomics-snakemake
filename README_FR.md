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

Le **pipeline Metatranscriptomics Snakemake** est un flux de travail modulaire et reproductible conçu pour traiter, assembler et analyser des **données métatranscriptomiques shotgun à lectures appariées Illumina**. Il automatise l’ensemble de l’analyse — du traitement des lectures brutes à la quantification de l’expression génique — en utilisant des outils de bioinformatique largement reconnus et intégrés via **Snakemake**.  Ce pipeline génère des assemblages de haute qualité, des profils taxonomiques et de gènes antimicrobiens, ainsi que des tableaux quantitatifs d’expression génique adaptés aux analyses statistiques et fonctionnelles en aval.

Le pipeline comprend **quatre étapes principales** :

- **Traitement des lectures d’échantillons** — Filtrage de qualité, suppression de la contamination par l’hôte et PhiX, et déplétion de l’ARN ribosomique (ARNr) à l’aide de *fastp*, *Bowtie2* et *SortMeRNA*.
- **Analyse des lectures courtes** — Classification taxonomique avec *Kraken2* (en utilisant *GTDB*) et profilage de la résistance antimicrobienne avec *RGI* (en utilisant *CARD*).
- **Assemblage individuel des échantillons** — Assemblage des transcrits avec *RNA SPAdes* et évaluation de la qualité via *rnaQUAST*.
- **Co-assemblage et quantification de l’expression** — Co-assemblage global avec *MEGAHIT*, suivi de l’alignement des lectures (*Bowtie2*), de l’évaluation de la couverture (*SAMtools*), de la prédiction des gènes (*Prodigal*) et de la quantification des régions codantes (*FeatureCounts*).

  💡 *Si des données de séquençage métagénomique sont disponibles pour les mêmes échantillons, les lectures métagénomiques nettoyées (après suppression des séquences PhiX et de l’hôte) doivent être utilisées à l’étape de co-assemblage.*

Des **améliorations futures** sont prévues pour ce pipeline, notamment :

- L’intégration de *CoverM* pour le mappage des lectures métatranscriptomiques sur les métagénomes assemblés.
- L’ajout d’un *module d’analyse des CAZymes* pour l’annotation fonctionnelle des enzymes actives sur les glucides.

---

## Documentation

Pour les détails techniques, y compris les instructions d’installation et d’utilisation, veuillez consulter le [Guide de l’utilisateur](/docs/user-guide.md).

---

## Remerciements

- **Crédits** : Ce projet a été réalisé au *Centre de recherche et de développement de Lacombe, Agriculture et Agroalimentaire Canada (AAC)*, par **Katherine James-Gzyl**, avec l’aide de **Devin Holman** et d’**Arun Kommadath**.

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
