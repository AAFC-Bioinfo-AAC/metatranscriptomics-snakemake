'''
    Filename: sample_assembly.smk
    Author: Katherine James-Gzyl
    Date created: 2025/08/13
    Snakemake version: 9.6.0
'''
rule rna_spades:
    input:
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        fasta = f"{ASSEMBLIES_DIR}/{{sample}}.fasta",
    log:
        f"{LOG_DIR}/spades/{{sample}}.log"
    conda:
        "../envs/RNAspades.yaml"
    threads: config["rna_spades"].get("threads", 48)
    params:
        # SPAdes -m expects MB (megabytes)!
        memory = config["rna_spades"].get("memory", 64000)
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log})"
        mkdir -p "$(dirname {output.fasta})"

        # Per-job tmp base and rnaspades work dir
        tmpbase="${{TMPDIR:-/tmp}}"
        jobtmp="$tmpbase/spades_{wildcards.sample}_$RANDOM"
        mkdir -p "$jobtmp" || {{ echo "Failed to create job tmpdir $jobtmp" >> "{log}"; exit 1; }}
        echo "Using job tmpdir base: $jobtmp" >> "{log}"

        outdir=$(mktemp -d -p "$jobtmp" "rnaspades_{wildcards.sample}_XXXXXX") || {{ echo "Failed to create outdir under $jobtmp" >> "{log}"; exit 1; }}
        spades_tmp="$outdir/tmp"
        mkdir -p "$spades_tmp"

        cleanup() {{
            # Remove only if the temp directories look correct
            if [[ -n "$outdir" && -d "$outdir" && "$outdir" == "$jobtmp"/rnaspades_* ]]; then
            rm -rf -- "$outdir"
            fi
            if [[ -n "$jobtmp" && -d "$jobtmp" && "$jobtmp" == $tmpbase/spades_* ]]; then
            rm -rf -- "$jobtmp"
            fi
        }}
        trap cleanup EXIT

        # Run SPAdes but do not abort on nonzero exit; handle it and ensure fasta exists
        set +e
        spades.py --rna \
            -t {threads} \
            -m {params.memory} \
            --tmp-dir "$spades_tmp" \
            -1 "{input.R1}" \
            -2 "{input.R2}" \
            -o "$outdir" >> "{log}" 2>&1
        status=$?
        set -e

        if [ -s "$outdir/transcripts.fasta" ]; then
            mv -f "$outdir/transcripts.fasta" "{output.fasta}"
            echo "SPAdes finished (exit $status); transcripts written." >> "{log}"
        else
            echo "SPAdes failed or produced no transcripts (exit $status); creating empty placeholder {output.fasta}" >> "{log}"
            : > "{output.fasta}"
        fi
        """

rule rnaquast_busco:
    input:
        fasta = f"{ASSEMBLIES_DIR}/{{sample}}.fasta",
        busco_lineage = lambda wc: BUSCO_LINEAGES[wc.lineage]
    output:
        report_dir = directory(f"{RNAQUAST_DIR}/{{sample}}_{{lineage}}")
    log:
        f"{LOG_DIR}/rnaquast/{{sample}}_{{lineage}}.log"
    conda:
        "../envs/rnaquast.yaml"
    threads: config["rnaquast_busco"].get("threads", 4)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"

        if [ -s "{input.fasta}" ]; then
            mkdir -p "{output.report_dir}"
            rnaQUAST.py --transcripts "{input.fasta}" \
                --output "{output.report_dir}" \
                --threads {threads} \
                --busco "{input.busco_lineage}" \
                >> "{log}" 2>&1
        else
            echo "Skipping rnaQUAST for {wildcards.sample}: assembly is empty" >> "{log}"
            mkdir -p "{output.report_dir}"
        fi
        """