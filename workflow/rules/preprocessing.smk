'''
    Filename: preprocessing.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/24
    Updated: 2026/10/08
    Snakemake version: 9.20.0
'''

import os


def _host_index_files(wildcards):
    """Select a complete, prebuilt Bowtie2 host/PhiX index."""
    prefix = os.path.abspath(BOWTIE_INDEX)
    parts = ("1", "2", "3", "4", "rev.1", "rev.2")

    for suffix in ("bt2", "bt2l"):
        files = [f"{prefix}.{part}.{suffix}" for part in parts]
        if all(
            os.path.isfile(path) and os.path.getsize(path) > 0
            for path in files
        ):
            return files

    raise ValueError(
        f"No complete, nonempty .bt2 or .bt2l host/PhiX index found: {prefix}"
    )


for setting in ("cut_tail", "cut_front", "detect_adapter_for_pe"):
    if not isinstance(config.get("fastp", {}).get(setting, True), bool):
        raise ValueError(
            f"fastp.{setting} must be true or false, without quotes."
        )


rule fastp_pe:
    input:
        fastq1 = lambda wc: SAMPLES[wc.sample]["fastq_1"],
        fastq2 = lambda wc: SAMPLES[wc.sample]["fastq_2"]
    output:
        r1 = temp(f"{TRIMMED_DIR}/{{sample}}_r1.fastq.gz"),
        r2 = temp(f"{TRIMMED_DIR}/{{sample}}_r2.fastq.gz"),
        u1 = temp(f"{TRIMMED_DIR}/{{sample}}_u1.fastq.gz"),
        u2 = temp(f"{TRIMMED_DIR}/{{sample}}_u2.fastq.gz"),
        html = f"{LOG_DIR}/fastp/{{sample}}.fastp.html",
        json = f"{LOG_DIR}/fastp/{{sample}}.fastp.json"
    log:
        f"{LOG_DIR}/fastp/{{sample}}.fastp.log"
    params:
        cut_tail = "--cut_tail" if config.get("fastp", {}).get("cut_tail", True) else "",
        cut_front = "--cut_front" if config.get("fastp", {}).get("cut_front", True) else "",
        detect_adapter = "--detect_adapter_for_pe" if config.get("fastp", {}).get("detect_adapter_for_pe", True) else "",
        cut_mean_quality = config.get("fastp", {}).get("cut_mean_quality", 20),
        cut_window_size = config.get("fastp", {}).get("cut_window_size", 4),
        qualified_quality_phred = config.get("fastp", {}).get("qualified_quality_phred", 15),
        length_required = config.get("fastp", {}).get("length_required", 100)
    threads: config.get("fastp", {}).get("threads", 2)
    conda:
        "../envs/fastp.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"

        fastp \
            --in1 {input.fastq1:q} \
            --in2 {input.fastq2:q} \
            --out1 {output.r1:q} \
            --out2 {output.r2:q} \
            --unpaired1 {output.u1:q} \
            --unpaired2 {output.u2:q} \
            {params.cut_tail} {params.cut_front} {params.detect_adapter} \
            --cut_mean_quality {params.cut_mean_quality:q} \
            --cut_window_size {params.cut_window_size:q} \
            --qualified_quality_phred {params.qualified_quality_phred:q} \
            --length_required {params.length_required:q} \
            --json {output.json:q} \
            --html {output.html:q} \
            --thread {threads} \
            > {log:q} 2>&1
        """


rule bowtie2_align:
    input:
        r1 = f"{TRIMMED_DIR}/{{sample}}_r1.fastq.gz",
        r2 = f"{TRIMMED_DIR}/{{sample}}_r2.fastq.gz",
        idx = _host_index_files
    output:
        bam = temp(f"{TRIMMED_DIR}/bam/{{sample}}.bam")
    params:
        r1_path = lambda wc, input: os.path.abspath(input.r1),
        r2_path = lambda wc, input: os.path.abspath(input.r2),
        index_suffix = lambda wc, input: input.idx[0].rsplit(".", 1)[1],
        rg_id = lambda wc: wc.sample,
        rg_sm = lambda wc: f"SM:{wc.sample}"
    log:
        f"{LOG_DIR}/bowtie2/{{sample}}.log"
    threads: config.get("bowtie2_align", {}).get("threads", 12)
    conda:
        "../envs/bowtie2.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})"
        : > {log:q}

        if (( {threads} < 2 )); then
            echo "ERROR: bowtie2_align requires at least 2 threads; received {threads}." >> {log:q}
            exit 1
        fi

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"
        tmpjob=$(mktemp -d "$tmpbase/host_alignment.XXXXXX")
        tmpjob=$(cd "$tmpjob" && pwd)
        trap 'rm -rf -- "$tmpjob"' EXIT

        # Expose only the selected index format, avoiding stale index files.
        index_files=( {input.idx:q} )
        index_parts=( 1 2 3 4 rev.1 rev.2 )

        for i in "${{!index_files[@]}}"; do
            ln -s "${{index_files[$i]}}" \
                "$tmpjob/host.${{index_parts[$i]}}.{params.index_suffix}"
        done

        # Simple relative aliases also avoid Bowtie2's comma-list parsing
        # when the original read paths contain commas.
        ln -s {params.r1_path:q} "$tmpjob/R1.fastq.gz"
        ln -s {params.r2_path:q} "$tmpjob/R2.fastq.gz"

        # Reserve SAMtools sort's main thread alongside its worker threads.
        bt2_threads=$(( {threads} / 2 ))
        sort_extra=$(( {threads} - bt2_threads - 1 ))

        (
            cd "$tmpjob"

            bowtie2 \
                -x host \
                -1 R1.fastq.gz \
                -2 R2.fastq.gz \
                --threads "$bt2_threads" \
                --rg-id {params.rg_id:q} \
                --rg {params.rg_sm:q}
        ) 2>> {log:q} \
        | samtools sort \
            -@ "$sort_extra" \
            -T "$tmpjob/sort" \
            -O BAM \
            -o {output.bam:q} \
            - 2>> {log:q}

        samtools quickcheck -v {output.bam:q} >> {log:q} 2>&1
        """


rule extract_unmapped_fastq:
    input:
        bam = f"{TRIMMED_DIR}/bam/{{sample}}.bam"
    output:
        r1 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        r2 = f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    log:
        f"{LOG_DIR}/bedtools/{{sample}}.log"
    threads: config.get("extract_unmapped_fastq", {}).get("threads", 8)
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1:q})" "$(dirname {log:q})"
        : > {log:q}

        if (( {threads} < 5 )); then
            echo "ERROR: extract_unmapped_fastq requires at least 5 threads; received {threads}." >> {log:q}
            exit 1
        fi

        # Reserve one CPU each for view, sort, BEDTools and two compressors.
        sort_extra=$(( {threads} - 5 ))

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"
        tmpjob=$(mktemp -d "$tmpbase/bam2fq.XXXXXX")

        pid1=""
        pid2=""

        cleanup() {{
            for pid in "$pid1" "$pid2"; do
                if [[ -n "$pid" ]]; then
                    kill "$pid" 2>/dev/null || true
                    wait "$pid" 2>/dev/null || true
                fi
            done

            rm -rf -- "$tmpjob"
        }}
        trap cleanup EXIT

        mkfifo "$tmpjob/r1.fifo" "$tmpjob/r2.fifo"

        pigz -p 1 --fast \
            < "$tmpjob/r1.fifo" \
            > "$tmpjob/r1.fastq.gz" \
            2>> {log:q} &
        pid1=$!

        pigz -p 1 --fast \
            < "$tmpjob/r2.fifo" \
            > "$tmpjob/r2.fastq.gz" \
            2>> {log:q} &
        pid2=$!

        # Both mates must be unmapped.
        # Exclude secondary and supplementary records.
        # BEDTools requires name-sorted input for paired fastq output.
        samtools view \
            -u -f 12 -F 2304 \
            {input.bam:q} 2>> {log:q} \
        | samtools sort \
            -n -@ "$sort_extra" \
            -T "$tmpjob/sort" \
            -l 0 -O BAM \
            - 2>> {log:q} \
        | bedtools bamtofastq \
            -i - \
            -fq "$tmpjob/r1.fifo" \
            -fq2 "$tmpjob/r2.fifo" \
            2>> {log:q}

        compression_status=0

        if wait "$pid1"; then
            pid1=""
        else
            compression_status=$?
            pid1=""
        fi

        if wait "$pid2"; then
            pid2=""
        else
            compression_status=$?
            pid2=""
        fi

        if (( compression_status != 0 )); then
            echo "ERROR: fastq compression failed (exit $compression_status)." >> {log:q}
            exit "$compression_status"
        fi

        # Publish the files only after both gzip streams have completed.
        mv -- "$tmpjob/r1.fastq.gz" {output.r1:q}
        mv -- "$tmpjob/r2.fastq.gz" {output.r2:q}
        """
