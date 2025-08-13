'''
    Filename: sortmerna.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/16
    Snakemake version: 9.6.0
'''
rule sortmerna_pe:
    input:
        ref=f"{RRNA_DB}",
        r1=f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        r2=f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        r1_out=f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        r2_out=f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz",
        stats=f"{RRNA_DEP_DIR}/{{sample}}_sortmerna_pe.stats"
    log:
        f"{LOG_DIR}/sortmerna/{{sample}}reads_pe.log"
    threads: config["sortmerna_pe"].get("threads", 48)
    conda:
        "../envs/sortmerna.yaml"
    shell:
        r"""
        set -euo pipefail

        export TMPDIR="${{TMPDIR:-/tmp}}"
        WORKDIR="$TMPDIR/smrn_{wildcards.sample}_$RANDOM"
        mkdir -p "$WORKDIR" || {{ echo "Failed to create TMP workdir $WORKDIR" >> {log}; exit 1; }}

        echo "Temporary workdir: $WORKDIR" >> {log}

        # Actually run SortMeRNA
        sortmerna \
            --ref {input.ref} \
            --reads {input.r1} \
            --reads {input.r2} \
            --aligned "$WORKDIR/aligned_reads" \
            --other "$WORKDIR/other_reads" \
            --fastx \
            --idx-dir idx \
            --paired_in \
            --out2 \
            --workdir "$WORKDIR" \
            --threads {threads} \
            >> {log} 2>&1

        # Move/copy output files to expected workflow outputs.
        cp "$WORKDIR/other_reads_fwd.fq.gz" {output.r1_out}
        cp "$WORKDIR/other_reads_rev.fq.gz" {output.r2_out}

        # Save the main SortMeRNA alignment/summary log as the .stats file
        if [ -f "$WORKDIR/aligned_reads.log" ]; then
            cp "$WORKDIR/aligned_reads.log" {output.stats}
            # Add a footer
            echo -e "\nSortMeRNA run finished at $(date)" >> {output.stats}
        else
            echo "SortMeRNA stats not found; see main log for errors." > {output.stats}
        fi

        echo "Temporary directory to be removed: $WORKDIR" >> {log}
        rm -rf "$WORKDIR"
        """