'''
    Filename: amr_short_reads.smk
    Original author: Katherine James-Gzyl
    Date created: 2025/07/16
    Updated: 2026/10/08
    Workflow controller version: Snakemake 9.20.0
'''

import os


def project_abspath(path):
    """Resolve relative paths against PROJECT_ROOT."""
    path = os.path.expandvars(os.path.expanduser(os.fspath(path)))
    if not os.path.isabs(path):
        path = os.path.join(PROJECT_ROOT, path)
    return os.path.abspath(path)


def resolve_card_db():
    """Use RGI_CARD when set; otherwise use config.card_latest."""
    card_db = os.getenv("RGI_CARD", "").strip() or config.get("card_latest")

    if not isinstance(card_db, str) or not card_db.strip():
        raise ValueError(
            "Set RGI_CARD or card_latest to the preloaded CARD "
            "database directory."
        )

    return project_abspath(card_db.strip())


def card_db_file(relative_path):
    """Resolve database inputs only when an AMR target is requested."""
    return lambda wildcards: os.path.join(
        resolve_card_db(), relative_path
    )


# Shared database loading and indexing are performed once during setup,
# before sample jobs are submitted.
rule rgi_validate_database:
    input:
        card_json = card_db_file("card.json"),
        card_fasta = card_db_file("card_reference.fasta"),
        metadata = card_db_file("loaded_databases.json"),
        kma_comp = card_db_file("bwt/card_reference/kma.comp.b"),
        kma_length = card_db_file("bwt/card_reference/kma.length.b"),
        kma_name = card_db_file("bwt/card_reference/kma.name"),
        kma_seq = card_db_file("bwt/card_reference/kma.seq.b")
    output:
        marker = f"{LOG_DIR}/rgi_card_db.validated"
    params:
        card_db = lambda wildcards: resolve_card_db()
    log:
        f"{LOG_DIR}/rgi/rgi_validate_database.log"
    conda:
        "../envs/rgi.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" \
                 "$(dirname {output.marker:q})"
        : > {log:q}

        for required_file in {input:q}; do
            if [[ ! -s "$required_file" ]]; then
                echo "ERROR: Required database file is missing or empty: $required_file" >> {log:q}
                exit 1
            fi
        done

        # Confirm that the JSON files are readable and their CARD versions agree.
        python -c '
import json, sys

with open(sys.argv[1]) as stream:
    card = json.load(stream)

with open(sys.argv[2]) as stream:
    metadata = json.load(stream)

version = card.get("_version")
if not version or version == "N/A":
    raise SystemExit("ERROR: card.json has no valid CARD version.")

for key in ("card_json", "card_canonical"):
    if metadata.get(key, {{}}).get("data_version") != version:
        raise SystemExit("ERROR: CARD metadata versions do not match card.json.")
' \
            {input.card_json:q} \
            {input.metadata:q} >> {log:q} 2>&1

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"

        validation_dir="$(mktemp -d "$tmpbase/rgi_validate.XXXXXX")"
        trap 'rm -rf -- "$validation_dir"' EXIT

        ln -s {params.card_db:q} "$validation_dir/localDB"

        (
            cd "$validation_dir"
            rgi database --version --local
        ) >> {log:q} 2>&1

        touch {output.marker:q}
        """


rule rgi_bwt:
    wildcard_constraints:
        sample = "[^/]+"
    input:
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz",
        database_marker = f"{LOG_DIR}/rgi_card_db.validated"
    output:
        json = temp(
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.allele_mapping_data.json"
        ),
        bam = temp(
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.sorted.length_100.bam"
        ),
        bai = temp(
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.sorted.length_100.bam.bai"
        ),
        allele = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.allele_mapping_data.txt"
        ),
        gene = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.gene_mapping_data.txt"
        ),
        artifacts_stats = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.artifacts_mapping_stats.txt"
        ),
        overall_stats = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.overall_mapping_stats.txt"
        ),
        reference_stats = (
            f"{CARD_RGI_OUTPUT_DIR}/{{sample}}/"
            f"{{sample}}_paired.reference_mapping_stats.txt"
        )
    params:
        card_db = lambda wildcards: resolve_card_db(),
        read1 = lambda wildcards, input: os.path.abspath(input.R1),
        read2 = lambda wildcards, input: os.path.abspath(input.R2),
        outprefix = lambda wildcards: project_abspath(
            f"{CARD_RGI_OUTPUT_DIR}/{wildcards.sample}/"
            f"{wildcards.sample}_paired"
        ),
        log_file = lambda wildcards: os.path.abspath(
            f"{LOG_DIR}/rgi/bwt_{wildcards.sample}.log"
        )
    log:
        f"{LOG_DIR}/rgi/bwt_{{sample}}.log"
    threads:
        config.get("rgi_bwt", {}).get("threads", 20)
    conda:
        "../envs/rgi.yaml"
    shell:
        r"""
        set -euo pipefail

        outprefix={params.outprefix:q}
        log_file={params.log_file:q}

        mkdir -p "$(dirname "$outprefix")" \
                 "$(dirname "$log_file")"
        : > "$log_file"

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"

        run_dir="$(mktemp -d "$tmpbase/rgi_bwt.XXXXXX")"
        run_dir="$(cd "$run_dir" && pwd)"
        trap 'rm -rf -- "$run_dir"' EXIT

        # RGI 6.0.4 constructs some internal shell commands without quoting.
        # Use simple local filenames and a shell-safe working-directory path.
        if [[ "$run_dir" =~ [^a-zA-Z0-9_./-] ]]; then
            echo "ERROR: TMPDIR must resolve to a path containing only letters, digits, underscores, dots, slashes or hyphens for RGI 6.0.4." >> "$log_file"
            exit 1
        fi

        ln -s {params.card_db:q} "$run_dir/localDB"
        ln -s {params.read1:q} "$run_dir/R1.fastq.gz"
        ln -s {params.read2:q} "$run_dir/R2.fastq.gz"

        (
            cd "$run_dir"

            rgi bwt \
                -1 R1.fastq.gz \
                -2 R2.fastq.gz \
                -a kma \
                -n {threads} \
                -o result \
                --local \
                --clean
        ) >> "$log_file" 2>&1

        # Catch RGI errors that may be logged despite a zero exit status.
        if grep -Eq '^(ERROR|CRITICAL)([[:space:]]|:)' "$log_file"; then
            echo "ERROR: RGI logged an error; inspect the complete log." >> "$log_file"
            exit 1
        fi

        suffixes=(
            allele_mapping_data.json
            sorted.length_100.bam
            sorted.length_100.bam.bai
            allele_mapping_data.txt
            gene_mapping_data.txt
            artifacts_mapping_stats.txt
            overall_mapping_stats.txt
            reference_mapping_stats.txt
        )

        # Validate fresh results before copying them to the output directory.
        # Reports containing only headers are allowed when no matches are found.
        for suffix in "${{suffixes[@]}}"; do
            if [[ ! -s "$run_dir/result.$suffix" ]]; then
                echo "ERROR: Expected RGI output is missing or empty: result.$suffix" >> "$log_file"
                exit 1
            fi
        done

        samtools quickcheck \
            "$run_dir/result.sorted.length_100.bam" \
            >> "$log_file" 2>&1

        for suffix in "${{suffixes[@]}}"; do
            cp -- "$run_dir/result.$suffix" "$outprefix.$suffix"
        done
        """
