'''
    Filename: sortmerna.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/16
    Updated: 2026/10/09
    Snakemake version: 9.20.0
'''


rule sortmerna_pe:
    input:
        ref = f"{RRNA_DB}",
        r1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        r2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    output:
        r1_out = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        r2_out = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz",
        stats = f"{RRNA_DEP_DIR}/{{sample}}_sortmerna_pe.stats"
    log:
        f"{LOG_DIR}/sortmerna/{{sample}}reads_pe.log"
    threads: config.get("sortmerna_pe", {}).get("threads", 24)
    conda:
        "../envs/sortmerna.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.r1_out:q})"
        : > {log:q}

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"
        WORKDIR=$(mktemp -d "$tmpbase/smrn.XXXXXX")
        trap 'rm -rf -- "$WORKDIR"' EXIT

        echo "Temporary workdir: $WORKDIR" >> {log:q}

        # Each job owns its index, alignment database and intermediate files.
        # --paired_in removes both mates if either matches the rRNA reference.
        sortmerna \
            --ref {input.ref:q} \
            --reads {input.r1:q} \
            --reads {input.r2:q} \
            --aligned "$WORKDIR/aligned_reads" \
            --other "$WORKDIR/other_reads" \
            --fastx \
            --zip-out 1 \
            --idx-dir "$WORKDIR/idx" \
            --paired_in \
            --out2 \
            --workdir "$WORKDIR" \
            --threads {threads} \
            >> {log:q} 2>&1

        if [[ ! -f "$WORKDIR/other_reads_fwd.fq.gz" || \
              ! -f "$WORKDIR/other_reads_rev.fq.gz" ]]; then
            echo "ERROR: SortMeRNA did not produce both retained paired-read files." >> {log:q}
            exit 1
        fi

        # Test gzip integrity before publishing outputs. Valid empty files pass.
        pigz -p {threads} -t \
            "$WORKDIR/other_reads_fwd.fq.gz" \
            "$WORKDIR/other_reads_rev.fq.gz" >> {log:q} 2>&1

        if [[ ! -s "$WORKDIR/aligned_reads.log" ]]; then
            echo "ERROR: SortMeRNA did not produce a nonempty summary log." >> {log:q}
            exit 1
        fi

        cp -- "$WORKDIR/other_reads_fwd.fq.gz" {output.r1_out:q}
        cp -- "$WORKDIR/other_reads_rev.fq.gz" {output.r2_out:q}
        cp -- "$WORKDIR/aligned_reads.log" {output.stats:q}

        echo "SortMeRNA completed; paired reads and summary published." >> {log:q}
        """
