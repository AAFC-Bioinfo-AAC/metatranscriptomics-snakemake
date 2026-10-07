"""Optional shared-reference CAZyme annotation and RNA fragment count tables."""
import hashlib

CAZYME_SETTINGS = config.get("cazyme", {})
if not isinstance(CAZYME_SETTINGS, dict):
    raise ValueError("cazyme must be a YAML mapping.")
CAZYME_ENABLED = CAZYME_SETTINGS.get("enabled", False)
if not isinstance(CAZYME_ENABLED, bool):
    raise ValueError("cazyme.enabled must be a YAML boolean (true or false).")
CAZYME_MIN_TOOLS = CAZYME_SETTINGS.get("min_tools", 2)
if isinstance(CAZYME_MIN_TOOLS, bool) or not isinstance(CAZYME_MIN_TOOLS, int) or CAZYME_MIN_TOOLS not in (2, 3):
    raise ValueError("cazyme.min_tools must be 2 or 3.")
if not isinstance(CAZYME_SETTINGS.get("database_release", "unspecified"), str) or not CAZYME_SETTINGS.get("database_release", "unspecified").strip():
    raise ValueError("cazyme.database_release must be a nonempty string.")
CAZYME_DIR = relpath("dbcan_output_dir")
CAZYME_TARGETS = [f"{CAZYME_DIR}/annotation", f"{CAZYME_DIR}/cazyme_annotations.tsv",
                  f"{CAZYME_DIR}/cazyme_gene_counts.tsv", f"{CAZYME_DIR}/cazyme_family_counts.tsv",
                  f"{CAZYME_DIR}/count_summary.json"]


def cazyme_database_inputs(wildcards):
    # Resolve lazily: unrelated targets can run without a CAZyme database.
    if not config.get("dbcan_DB_path"):
        raise ValueError("Set dbcan_DB_path to a prepared dbCAN database directory for CAZyme targets.")
    directory = relpath("dbcan_DB_path")
    return [os.path.join(directory, name) for name in
            ("CAZy.dmnd", "dbCAN.hmm", "dbCAN-sub.hmm", "fam-substrate-mapping.tsv")]


def cazyme_script_checksum(name):
    # Shell-invoked helper changes must also trigger reruns under the SLURM profile.
    with open(os.path.join(PIPELINE_DIR, "workflow", "scripts", name), "rb") as stream:
        return hashlib.sha256(stream.read()).hexdigest()


rule cazyme_all:
    input:
        CAZYME_TARGETS


rule prepare_cazyme_proteins:
    input:
        proteins = f"{PRODIGAL_DIR}/coassembly.faa",
        gff = f"{PRODIGAL_DIR}/coassembly.gff"
    output:
        proteins = f"{CAZYME_DIR}/reference_proteins.faa",
        mapping = f"{CAZYME_DIR}/protein_gene_ids.tsv"
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
          --output {output.proteins:q} --mapping {output.mapping:q} > {log:q} 2>&1
        """


rule cazyme_annotation:
    input:
        proteins = f"{CAZYME_DIR}/reference_proteins.faa",
        databases = cazyme_database_inputs
    output:
        annotation = directory(f"{CAZYME_DIR}/annotation")
    params:
        script = f"{PIPELINE_DIR}/workflow/scripts/run_cazyme_annotation.py",
        script_checksum = cazyme_script_checksum("run_cazyme_annotation.py"),
        db_dir = lambda wc, input: os.path.dirname(input.databases[0]),
        database_release = CAZYME_SETTINGS.get("database_release", "unspecified")
    threads: CAZYME_SETTINGS.get("threads", 8)
    resources:
        mem_mb = 16000
    conda:
        "../envs/dbcan.yaml"
    log:
        f"{LOG_DIR}/dbcan/cazyme_annotation.log"
    shell:
        """
        python {params.script:q} --proteins {input.proteins:q} --db-dir {params.db_dir:q} \
          --output {output.annotation:q} --log {log:q} --threads {threads} \
          --database-release {params.database_release:q}
        """


rule cazyme_rna_counts:
    input:
        annotation = f"{CAZYME_DIR}/annotation",
        mapping = f"{CAZYME_DIR}/protein_gene_ids.tsv",
        counts = expand(f"{FEATURECOUNTS_DIR}/{{sample}}_counts.txt", sample=SAMPLE_NAMES)
    output:
        annotations = f"{CAZYME_DIR}/cazyme_annotations.tsv",
        genes = f"{CAZYME_DIR}/cazyme_gene_counts.tsv",
        families = f"{CAZYME_DIR}/cazyme_family_counts.tsv",
        summary = f"{CAZYME_DIR}/count_summary.json"
    params:
        script = f"{PIPELINE_DIR}/workflow/scripts/summarize_cazyme_counts.py",
        script_checksum = cazyme_script_checksum("summarize_cazyme_counts.py"),
        samples = SAMPLE_NAMES,
        min_tools = CAZYME_MIN_TOOLS
    threads: 1
    resources:
        mem_mb = 8000
    conda:
        "../envs/cazyme_tables.yaml"
    log:
        f"{LOG_DIR}/dbcan/cazyme_rna_counts.log"
    shell:
        """
        python {params.script:q} --overview {input.annotation:q}/overview.tsv \
          --mapping {input.mapping:q} --samples {params.samples:q} --counts {input.counts:q} \
          --min-tools {params.min_tools} --annotations {output.annotations:q} \
          --gene-counts {output.genes:q} --family-counts {output.families:q} \
          --summary {output.summary:q} > {log:q} 2>&1
        """
