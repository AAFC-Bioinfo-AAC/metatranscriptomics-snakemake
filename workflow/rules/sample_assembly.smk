'''
    Filename: sample_assembly.smk
    Author: Katherine James-Gzyl
    Date created: 2025/08/13
    Snakemake version: 9.6.0
'''
## Need to work on making RNA samples an optional rule. If the assembly is bad then keep going with pipeline
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
        memory = config["rna_spades"].get("memory", "64000")
    shell:
        r"""
        set -euo pipefail

        export TMPDIR="${{TMPDIR:-/tmp}}/spades_{wildcards.sample}_$RANDOM"
        mkdir -p "$TMPDIR" || (echo "Failed to create TMPDIR $TMPDIR" && exit 1)
        echo "Using TMPDIR: $TMPDIR" >> {log}

        
        outdir=$(mktemp -d "$TMPDIR/rnaspades_{wildcards.sample}_XXXXXX")
        spades_tmp="$outdir/tmp"
        mkdir -p "$spades_tmp"

        spades.py --rna \
            -t {threads} \
            -m {params.memory} \
            --tmp-dir "$spades_tmp" \
            -1 {input.R1} \
            -2 {input.R2} \
            -o "$outdir" >> {log} 2>&1 || true

         # Always produce an output
        if [ ! -s "$outdir/transcripts.fasta" ]; then
            echo "Failed to assemble for {wildcards.sample}: transcripts.fasta not produced or empty." >> {log} 2>&1
            echo ">dummy_sequence" > {output.fasta}
        else
            cp "$outdir/transcripts.fasta" {output.fasta}
        fi

        echo "Temporary directory to be removed: $outdir" >> {log}
        rm -rf "$outdir"
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
        "envs/rnaquast.yaml"
    threads: config["rnaquast_busco"].get("threads", 4)
    shell:
        r"""
        set -euo pipefail

        rnaQUAST.py --transcripts {input.fasta} \
            --output {output.report_dir} \
            --threads {threads} \
            --busco {input.busco_lineage} \
            &>> {log}
        """