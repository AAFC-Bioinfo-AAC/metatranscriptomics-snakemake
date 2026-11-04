"""Shared RNA-reference CAZyme analysis using the metagenomic dbCAN strategy."""
import hashlib

CAZYME_SETTINGS = config.get("cazyme", {})
if not isinstance(CAZYME_SETTINGS, dict):
    raise ValueError("cazyme must be a YAML mapping.")
CAZYME_MIN_TOOLS = CAZYME_SETTINGS.get("min_tools", 2)
if isinstance(CAZYME_MIN_TOOLS, bool) or not isinstance(CAZYME_MIN_TOOLS, int) or CAZYME_MIN_TOOLS not in (2, 3):
    raise ValueError("cazyme.min_tools must be 2 or 3.")
if not isinstance(CAZYME_SETTINGS.get("database_release", "unspecified"), str) or not CAZYME_SETTINGS.get("database_release", "unspecified").strip():
    raise ValueError("cazyme.database_release must be a nonempty string.")
CAZYME_DIR = relpath("dbcan_output_dir")
DBCAN_SHARED = os.path.join(CAZYME_DIR, "coassembly")
DBCAN_REFERENCE = COASSEMBLY_FASTA
DBCAN_ANNOTATION = os.path.join(DBCAN_SHARED, "annotation")
DBCAN_QUANTIFICATION = os.path.join(DBCAN_SHARED, "quantification_annotation")
DBCAN_BWA_PREFIX = os.path.join(DBCAN_SHARED, "bwa_index", "coassembly")
DBCAN_BWA_INDEX = multiext(DBCAN_BWA_PREFIX, ".amb", ".ann", ".bwt", ".pac", ".sa")
CAZYME_BAMS = expand(f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam", sample=SAMPLE_NAMES)
CAZYME_TARGETS = [DBCAN_ANNOTATION] + [
    f"{DBCAN_SHARED}/samples/{sample}/{sample}_abund/{name}"
    for sample in SAMPLE_NAMES
    for name in ("fam_abund.out", "subfam_abund.out", "EC_abund.out",
                 "fam_substrate_abund.out", "CGC_abund.out",
                 "CGC_substrate_PUL_homology.out", "CGC_substrate_majority_voting.out")
]
CAZYME_TARGETS += [f"{DBCAN_SHARED}/{name}" for name in (
    "cazyme_family_RPM.tsv", "cazyme_subfamily_RPM.tsv",
    "cazyme_substrate_RPM.tsv", "count_summary.json")]


def cazyme_database_inputs(wildcards):
    # Resolve lazily: unrelated targets can run without a CAZyme database.
    if not config.get("dbcan_DB_path"):
        raise ValueError("Set dbcan_DB_path to a prepared dbCAN database directory for CAZyme targets.")
    directory = relpath("dbcan_DB_path")
    return [os.path.join(directory, name) for name in
            ("CAZy.dmnd", "dbCAN.hmm", "dbCAN-sub.hmm", "fam-substrate-mapping.tsv",
             "TCDB.dmnd", "TF.hmm", "TF.dmnd", "STP.hmm", "PUL.dmnd", "dbCAN-PUL.xlsx", "peptidase_db.dmnd", "sulfatlas_db.dmnd")]


def cazyme_script_checksum(name):
    # Shell-invoked helper changes must also trigger reruns under the SLURM profile.
    with open(os.path.join(PIPELINE_DIR, "workflow", "scripts", name), "rb") as stream:
        return hashlib.sha256(stream.read()).hexdigest()


rule cazyme_all:
    input:
        CAZYME_TARGETS


rule merge_cazyme_rpm:
    input:
        files = expand(
            f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/{{name}}",
            sample=SAMPLE_NAMES,
            name=["fam_abund.out", "subfam_abund.out", "fam_substrate_abund.out"]
        )
    output:
        families = f"{DBCAN_SHARED}/cazyme_family_RPM.tsv",
        subfamilies = f"{DBCAN_SHARED}/cazyme_subfamily_RPM.tsv",
        substrates = f"{DBCAN_SHARED}/cazyme_substrate_RPM.tsv",
        summary = f"{DBCAN_SHARED}/count_summary.json"
    params:
        script = f"{PIPELINE_DIR}/workflow/scripts/merge_cazyme_rpm.py",
        script_checksum = cazyme_script_checksum("merge_cazyme_rpm.py"),
        samples = SAMPLE_NAMES,
        directories = [f"{DBCAN_SHARED}/samples/{s}/{s}_abund" for s in SAMPLE_NAMES],
        output_dir = DBCAN_SHARED
    conda:
        "../envs/cazyme_tables.yaml"
    log:
        f"{LOG_DIR}/dbcan/merge_cazyme_rpm.log"
    shell:
        """
        python {params.script:q} --samples {params.samples:q} --directories {params.directories:q} \
          --output {params.output_dir:q} --summary {output.summary:q} > {log:q} 2>&1
        """


rule prepare_cazyme_proteins:
    input:
        proteins = f"{PRODIGAL_DIR}/coassembly.faa",
        gff = f"{PRODIGAL_DIR}/coassembly.gff",
        saf = f"{PRODIGAL_DIR}/coassembly.saf"
    output:
        proteins = f"{DBCAN_SHARED}/reference_proteins.faa",
        mapping = f"{DBCAN_SHARED}/protein_gene_ids.tsv"
    params:
        script = f"{PIPELINE_DIR}/workflow/scripts/prepare_cazyme_proteins.py",
        script_checksum = cazyme_script_checksum("prepare_cazyme_proteins.py")
    threads: 1
    resources:
        mem_mb = 2000
    conda:
        "../envs/cazyme_tables.yaml"
    log:
        f"{LOG_DIR}/dbcan/prepare_proteins.log"
    shell:
        """
        python {params.script:q} --proteins {input.proteins:q} --gff {input.gff:q} \
          --saf {input.saf:q} \
          --output {output.proteins:q} --mapping {output.mapping:q} > {log:q} 2>&1
        """


rule cazyme_annotation:
    input:
        proteins = f"{DBCAN_SHARED}/reference_proteins.faa",
        gff = f"{PRODIGAL_DIR}/coassembly.gff",
        databases = cazyme_database_inputs
    output:
        annotation = directory(f"{DBCAN_SHARED}/annotation")
    params:
        script = f"{PIPELINE_DIR}/workflow/scripts/run_cazyme_annotation.py",
        script_checksum = cazyme_script_checksum("run_cazyme_annotation.py"),
        db_dir = lambda wc, input: os.path.dirname(input.databases[0]),
        database_release = CAZYME_SETTINGS.get("database_release", "unspecified")
    threads: CAZYME_SETTINGS.get("threads", 8)
    resources:
        mem_mb = 64000
    conda:
        "../envs/dbcan.yaml"
    log:
        f"{LOG_DIR}/dbcan/cazyme_annotation.log"
    shell:
        """
        python {params.script:q} --proteins {input.proteins:q} --gff {input.gff:q} --db-dir {params.db_dir:q} \
          --output {output.annotation:q} --log {log:q} --threads {threads} \
          --database-release {params.database_release:q}
        """



def bwa_mapping_threads():
    value = config.get("bwa_mem_mapping", {}).get("threads", 12)
    if isinstance(value, bool) or not isinstance(value, int) or value < 2:
        raise ValueError("bwa_mem_mapping.threads must be an integer of at least 2.")
    return value


rule prepare_cazyme_quantification:
    input:
        annotation = DBCAN_ANNOTATION,
        mapping = f"{DBCAN_SHARED}/protein_gene_ids.tsv"
    output:
        directory(DBCAN_QUANTIFICATION)
    params:
        script = f"{PIPELINE_DIR}/workflow/scripts/prepare_cazyme_quantification.py",
        script_checksum = cazyme_script_checksum("prepare_cazyme_quantification.py"),
        min_tools = CAZYME_MIN_TOOLS
    conda:
        "../envs/cazyme_tables.yaml"
    log:
        f"{LOG_DIR}/dbcan/prepare_quantification.log"
    shell:
        """
        python {params.script:q} --annotation {input.annotation:q} --mapping {input.mapping:q} \
          --output {output:q} --min-tools {params.min_tools} > {log:q} 2>&1
        """


rule bwa_index_coassembly:
    input:
        assembly = DBCAN_REFERENCE
    output:
        DBCAN_BWA_INDEX
    params:
        prefix = DBCAN_BWA_PREFIX
    threads: 1
    log:
        f"{LOG_DIR}/dbcan/bwa_index_coassembly.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log:q})" "$(dirname {params.prefix:q})"
        bwa index -p {params.prefix:q} {input.assembly:q} > {log:q} 2>&1
        """


rule bwa_mem_mapping:
    input:
        assembly = DBCAN_REFERENCE,
        index = DBCAN_BWA_INDEX,
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        bam = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam",
        bai = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam.bai"
    params:
        index_prefix = DBCAN_BWA_PREFIX
    log:
        f"{LOG_DIR}/dbcan/bwa_mem_mapping/{{sample}}.log"
    threads:
        bwa_mapping_threads()
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.bam:q})"
        : > {log:q}

        if (( {threads} < 2 )); then
            echo "ERROR: bwa_mem_mapping requires at least 2 threads; received {threads}." >> {log:q}
            exit 1
        fi

        # Reserve approximately one quarter of the CPU budget for SAMtools.
        sort_threads=$(( {threads} / 4 ))
        if (( sort_threads < 1 )); then
            sort_threads=1
        fi
        bwa_threads=$(( {threads} - sort_threads ))

        # Reserve one CPU within the SAMtools budget for its main thread.
        sort_extra=$(( sort_threads - 1 ))

        echo "[$(date)] Total threads: {threads}; BWA: $bwa_threads; SAMtools budget: $sort_threads; SAMtools workers: $sort_extra" \
            >> {log:q}

        work=$(mktemp -d "$(dirname {output.bam:q})/.sort.XXXXXX")
        trap 'rm -rf -- "$work"' EXIT

        bwa mem \
            -t "$bwa_threads" \
            {params.index_prefix:q} \
            {input.R1:q} \
            {input.R2:q} \
            2>> {log:q} \
        | samtools sort \
            -@ "$sort_extra" -m 1G -T "$work/sort" \
            -o {output.bam:q} \
            - \
            >> {log:q} 2>&1

        samtools quickcheck -v {output.bam:q} >> {log:q} 2>&1
        samtools index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
        """


rule dbcan_depth:
    input:
        gff = f"{PRODIGAL_DIR}/coassembly.gff",
        bam = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam",
        bai = f"{DBCAN_SHARED}/samples/{{sample}}/mapping/{{sample}}.bam.bai"
    output:
        depth_file = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    log:
        f"{LOG_DIR}/dbcan/dbcan_depth/{{sample}}.log"
    params:
        overlap_base_ratio = config.get("dbcan_depth", {}).get(
            "overlap_base_ratio", 0.2
        ),
        mapping_quality = config.get("dbcan_depth", {}).get(
            "mapping_quality", 30
        ),
        identity = config.get("dbcan_depth", {}).get("identity", 0.98)
    threads:
        config.get("dbcan_depth", {}).get("threads", 4)
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.depth_file:q})"
        : > {log:q}

        dbcan_utils cal_coverage \
            -g {input.gff:q} \
            -i {input.bam:q} \
            -o {output.depth_file:q} \
            -t {threads} \
            --overlap_base_ratio {params.overlap_base_ratio} \
            --mapping_quality {params.mapping_quality} \
            --identity {params.identity} \
            >> {log:q} 2>&1
        """


# -------------------------------------------------------------------
# Abundance calculation
# -------------------------------------------------------------------

rule get_abundances_rpm:
    input:
        annotation = DBCAN_QUANTIFICATION,
        depth_file = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/{{sample}}.depth.txt"
    output:
        cazy_fam_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/fam_abund.out",
        cazy_subfam_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/subfam_abund.out",
        EC_number = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/EC_abund.out",
        substrate_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/fam_substrate_abund.out",
        cgc_ab = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/CGC_abund.out",
        substrate_ho = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/CGC_substrate_PUL_homology.out",
        cgc_substrate_voting = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/CGC_substrate_majority_voting.out"
    params:
        overview = f"{DBCAN_QUANTIFICATION}/overview.tsv",
        cgc_table = f"{DBCAN_QUANTIFICATION}/cgc_standard_out.tsv",
        substrate = f"{DBCAN_QUANTIFICATION}/substrate_prediction.tsv",
        dbcan_dir = DBCAN_QUANTIFICATION,
        ab_dir = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund",
        marker_file = f"{DBCAN_SHARED}/samples/{{sample}}/{{sample}}_abund/did_not_run_get_abundances_rpm.txt"
    log:
        f"{LOG_DIR}/dbcan/get_abundances_rpm/{{sample}}.log"
    conda:
        "../envs/dbcan.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" {params.ab_dir:q}
        : > {log:q}
        rm -f -- {params.marker_file:q} {output:q}

        # Return success only when a tabular file contains at least one
        # non-empty row after its header.
        has_data_rows() {{
            local input_file="$1"
            [[ -s "$input_file" ]] && \
                awk 'NR > 1 && NF {{ found = 1; exit }} END {{ exit !found }}' \
                    "$input_file"
        }}

        require_output() {{
            local output_file="$1"
            if [[ ! -e "$output_file" ]]; then
                echo "ERROR: Expected output was not created: $output_file" \
                    >> {log:q}
                exit 1
            fi
        }}

        record_skip() {{
            printf '%s\n' "$1" >> {params.marker_file:q}
        }}

        (
            cd {params.ab_dir:q}

            if has_data_rows {params.overview:q}; then
                dbcan_utils fam_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.cazy_fam_ab:q}
                require_output {output.cazy_subfam_ab:q}
                require_output {output.EC_number:q}

                dbcan_utils fam_substrate_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.substrate_ab:q}
            else
                : > {output.cazy_fam_ab:q}
                : > {output.cazy_subfam_ab:q}
                : > {output.EC_number:q}
                : > {output.substrate_ab:q}
                record_skip \
                    "Family-level abundance calculations were skipped because overview.tsv contained no annotation rows."
            fi

            if has_data_rows {params.cgc_table:q}; then
                dbcan_utils CGC_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.cgc_ab:q}
            else
                : > {output.cgc_ab:q}
                record_skip \
                    "CGC abundance calculation was skipped because cgc_standard_out.tsv contained no CGC rows."
            fi

            if has_data_rows {params.substrate:q}; then
                dbcan_utils CGC_substrate_abund \
                    -bt {input.depth_file:q} \
                    -i {params.dbcan_dir:q} \
                    -a RPM \
                    >> {log:q} 2>&1

                require_output {output.substrate_ho:q}
                require_output {output.cgc_substrate_voting:q}
            else
                : > {output.substrate_ho:q}
                : > {output.cgc_substrate_voting:q}
                record_skip \
                    "CGC substrate abundance calculations were skipped because substrate_prediction.tsv contained no prediction rows."
            fi
        )

        """
