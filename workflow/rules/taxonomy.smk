'''
    Filename: taxonomy.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/16
    Updated: 2026/10/09
    Snakemake version: 9.20.0
'''

import csv
import math
import os
from pathlib import Path


def _bracken_integer(key, default, minimum):
    value = config.get("bracken", {}).get(key, default)

    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < minimum
    ):
        raise ValueError(
            f"bracken.{key} must be an integer of at least {minimum}."
        )

    return value


BRACKEN_READLEN = _bracken_integer("readlen", 150, 1)

BRACKEN_THRESHOLDS = {
    "species": _bracken_integer("threshold_species", 10, 0),
    "genus": _bracken_integer("threshold_genus", 10, 0),
    "phylum": _bracken_integer("threshold_phylum", 10, 0),
}

_confidence = config.get("kraken2", {}).get("conf_threshold", 0.5)

if isinstance(_confidence, bool):
    raise ValueError(
        "kraken2.conf_threshold must be a number between 0 and 1."
    )

try:
    KRAKEN_CONFIDENCE = float(_confidence)
except (TypeError, ValueError):
    raise ValueError(
        "kraken2.conf_threshold must be a number between 0 and 1."
    )

if (
    not math.isfinite(KRAKEN_CONFIDENCE)
    or not 0 <= KRAKEN_CONFIDENCE <= 1
):
    raise ValueError(
        "kraken2.conf_threshold must be a number between 0 and 1."
    )


rule kraken2:
    wildcard_constraints:
        sample = "[^/]+"
    input:
        hash = f"{TAXONOMY_DB}/hash.k2d",
        opts = f"{TAXONOMY_DB}/opts.k2d",
        taxo = f"{TAXONOMY_DB}/taxo.k2d",
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt",
        kraken = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.kraken"
    log:
        f"{LOG_DIR}/kraken2/{{sample}}.log"
    conda:
        "../envs/kraken2.yaml"
    threads:
        config.get("kraken2", {}).get("threads", 2)
    params:
        db = os.path.abspath(TAXONOMY_DB),
        conf_threshold = KRAKEN_CONFIDENCE
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" \
            "$(dirname {output.report:q})"
        : > {log:q}

        # Decompression failures may not propagate through Kraken2.
        gzip -t {input.R1:q} {input.R2:q} >> {log:q} 2>&1

        kraken2 --use-names \
            --threads {threads} \
            --db {params.db:q} \
            --gzip-compressed \
            --confidence {params.conf_threshold} \
            --report-zero-counts \
            --paired {input.R1:q} {input.R2:q} \
            --report {output.report:q} \
            --output {output.kraken:q} \
            >> {log:q} 2>&1

        if [[ ! -s {output.report:q} || ! -f {output.kraken:q} ]]; then
            echo "ERROR: Kraken2 did not produce its report and classification file." >> {log:q}
            exit 1
        fi
        """


rule bracken:
    wildcard_constraints:
        sample = "[^/]+"
    input:
        distribution = f"{TAXONOMY_DB}/database{BRACKEN_READLEN}mers.kmer_distrib",
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt"
    output:
        species = f"{BRACKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt",
        genus = f"{BRACKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt",
        phylum = f"{BRACKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt",
        species_kreport = f"{BRACKEN_OUTPUT_DIR}/reports/{{sample}}_bracken.species.kreport.txt",
        genus_kreport = f"{BRACKEN_OUTPUT_DIR}/reports/{{sample}}_bracken.genus.kreport.txt",
        phylum_kreport = f"{BRACKEN_OUTPUT_DIR}/reports/{{sample}}_bracken.phylum.kreport.txt"
    log:
        f"{LOG_DIR}/bracken/{{sample}}.log"
    conda:
        "../envs/kraken2.yaml"
    threads:
        1
    params:
        readlen = BRACKEN_READLEN,
        threshold_species = BRACKEN_THRESHOLDS["species"],
        threshold_genus = BRACKEN_THRESHOLDS["genus"],
        threshold_phylum = BRACKEN_THRESHOLDS["phylum"],
        distribution_path = lambda wc, input: os.path.abspath(
            input.distribution
        ),
        report_path = lambda wc, input: os.path.abspath(
            input.report
        )
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" \
            "$(dirname {output.species:q})" \
            "$(dirname {output.genus:q})" \
            "$(dirname {output.phylum:q})" \
            "$(dirname {output.species_kreport:q})"
        : > {log:q}

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"
        WORKDIR=$(mktemp -d "$tmpbase/bracken.XXXXXX")
        trap 'rm -rf -- "$WORKDIR"' EXIT

        # Bracken 2.7's wrapper leaves paths unquoted internally.
        # Use simple relative aliases in a private working directory.
        mkdir -p "$WORKDIR/database"

        ln -s {params.distribution_path:q} \
            "$WORKDIR/database/database{params.readlen}mers.kmer_distrib"

        ln -s {params.report_path:q} \
            "$WORKDIR/kraken.report.txt"

        (
            cd "$WORKDIR"

            bracken -d database -i kraken.report.txt \
                -r {params.readlen} -l S \
                -t {params.threshold_species} \
                -o species.tsv -w species.kreport.txt

            bracken -d database -i kraken.report.txt \
                -r {params.readlen} -l G \
                -t {params.threshold_genus} \
                -o genus.tsv -w genus.kreport.txt

            bracken -d database -i kraken.report.txt \
                -r {params.readlen} -l P \
                -t {params.threshold_phylum} \
                -o phylum.tsv -w phylum.kreport.txt

            # The wrapper can exit successfully without producing results.
            for rank in species genus phylum; do
                if [[ ! -s "$rank.tsv" || ! -s "$rank.kreport.txt" ]]; then
                    echo "ERROR: Bracken did not produce both outputs for $rank."
                    exit 1
                fi
            done
        ) >> {log:q} 2>&1

        cp -- "$WORKDIR/species.tsv" {output.species:q}
        cp -- "$WORKDIR/genus.tsv" {output.genus:q}
        cp -- "$WORKDIR/phylum.tsv" {output.phylum:q}

        cp -- "$WORKDIR/species.kreport.txt" {output.species_kreport:q}
        cp -- "$WORKDIR/genus.kreport.txt" {output.genus_kreport:q}
        cp -- "$WORKDIR/phylum.kreport.txt" {output.phylum_kreport:q}
        """


rule combine_bracken_outputs:
    input:
        species = expand(
            f"{BRACKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt",
            sample=SAMPLES
        ),
        genus = expand(
            f"{BRACKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt",
            sample=SAMPLES
        ),
        phylum = expand(
            f"{BRACKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt",
            sample=SAMPLES
        )
    output:
        species = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_species.txt",
        genus = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_genus.txt",
        phylum = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_phylum.txt"
    log:
        f"{LOG_DIR}/bracken/combine_bracken_outputs.log"
    threads:
        1
    run:
        # Merge by taxonomy ID using Python's standard library.
        fields = [
            "name",
            "taxonomy_id",
            "taxonomy_lvl",
            "kraken_assigned_reads",
            "added_reads",
            "new_est_reads",
            "fraction_total_reads",
        ]

        ranks = (
            ("species", "S"),
            ("genus", "G"),
            ("phylum", "P"),
        )

        with open(log[0], "w") as log_handle:
            for level, rank in ranks:
                paths = list(input[level])
                labels = [Path(path).name for path in paths]

                if not paths or len(labels) != len(set(labels)):
                    raise ValueError(
                        f"Missing inputs or duplicate sample filenames "
                        f"at {level} level."
                    )

                counts_by_sample = []
                metadata = {}

                for path in paths:
                    sample_counts = {}

                    with open(path, newline="") as handle:
                        reader = csv.DictReader(
                            handle,
                            delimiter="\t"
                        )

                        if reader.fieldnames != fields:
                            raise ValueError(
                                f"Unexpected Bracken columns in {path}."
                            )

                        for row in reader:
                            if (
                                None in row
                                or any(
                                    value is None
                                    for value in row.values()
                                )
                            ):
                                raise ValueError(
                                    f"Incomplete Bracken row in {path}."
                                )

                            taxid = row["taxonomy_id"]
                            name = row["name"]

                            if (
                                not name
                                or not taxid.isdigit()
                                or int(taxid) < 1
                                or row["taxonomy_lvl"] != rank
                            ):
                                raise ValueError(
                                    f"Invalid taxon or rank in {path}."
                                )

                            if taxid in sample_counts:
                                raise ValueError(
                                    f"Duplicate taxonomy ID {taxid} "
                                    f"in {path}."
                                )

                            assigned, added, estimated = (
                                int(row[field])
                                for field in (
                                    "kraken_assigned_reads",
                                    "added_reads",
                                    "new_est_reads",
                                )
                            )

                            fraction = float(
                                row["fraction_total_reads"]
                            )

                            if (
                                min(assigned, added, estimated) < 0
                                or assigned + added != estimated
                            ):
                                raise ValueError(
                                    f"Invalid estimated counts in {path}."
                                )

                            if (
                                not math.isfinite(fraction)
                                or not 0 <= fraction <= 1
                            ):
                                raise ValueError(
                                    f"Invalid fraction in {path}."
                                )

                            if (
                                taxid in metadata
                                and metadata[taxid] != name
                            ):
                                raise ValueError(
                                    f"Conflicting names for taxonomy "
                                    f"ID {taxid}."
                                )

                            metadata[taxid] = name
                            sample_counts[taxid] = estimated

                    counts_by_sample.append(sample_counts)

                totals = [
                    sum(counts.values())
                    for counts in counts_by_sample
                ]

                Path(output[level]).parent.mkdir(
                    parents=True,
                    exist_ok=True
                )

                with open(
                    output[level],
                    "w",
                    newline=""
                ) as handle:
                    writer = csv.writer(
                        handle,
                        delimiter="\t",
                        lineterminator="\n"
                    )

                    columns = [
                        "name",
                        "taxonomy_id",
                        "taxonomy_lvl",
                    ]

                    for label in labels:
                        columns.extend([
                            label + "_num",
                            label + "_frac",
                        ])

                    writer.writerow(columns)

                    for taxid in sorted(metadata, key=int):
                        row = [metadata[taxid], taxid, rank]

                        for counts, total in zip(
                            counts_by_sample,
                            totals
                        ):
                            estimated = counts.get(taxid, 0)
                            fraction = (
                                estimated / total
                                if total else 0.0
                            )

                            row.extend([
                                estimated,
                                f"{fraction:.5f}",
                            ])

                        writer.writerow(row)

                print(
                    f"Merged {len(paths)} samples and "
                    f"{len(metadata)} taxa at {level} level.",
                    file=log_handle
                )


rule bracken_extract:
    input:
        species_table = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_species.txt",
        genus_table = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_genus.txt",
        phylum_table = f"{BRACKEN_OUTPUT_DIR}/merged_abundance_phylum.txt"
    output:
        species_raw = f"{BRACKEN_OUTPUT_DIR}/Bracken_species_raw_abundance.csv",
        species_rel = f"{BRACKEN_OUTPUT_DIR}/Bracken_species_relative_abundance.csv",
        genus_raw = f"{BRACKEN_OUTPUT_DIR}/Bracken_genus_raw_abundance.csv",
        genus_rel = f"{BRACKEN_OUTPUT_DIR}/Bracken_genus_relative_abundance.csv",
        phylum_raw = f"{BRACKEN_OUTPUT_DIR}/Bracken_phylum_raw_abundance.csv",
        phylum_rel = f"{BRACKEN_OUTPUT_DIR}/Bracken_phylum_relative_abundance.csv"
    conda:
        "../envs/kraken2.yaml"
    script:
        "../scripts/extract_bracken_columns.py"
