'''
    Filename: preprocessing.smk
    Author: Katherine James-Gzyl
    Date created: 2025/07/16
    Snakemake version: 9.6.0
'''
rule fastp_pe:
    input:
        sample = lambda wc: [find_read_file(wc.sample, 1), find_read_file(wc.sample, 2)]
    output:
        trimmed = [f"{TRIMMED_DIR}/{{sample}}_r1.fastq.gz", f"{TRIMMED_DIR}/{{sample}}_r2.fastq.gz"],
        unpaired1 = f"{TRIMMED_DIR}/{{sample}}_u1.fastq.gz",
        unpaired2 = f"{TRIMMED_DIR}/{{sample}}_u2.fastq.gz",
        html = f"{TRIMMED_DIR}/{{sample}}.fastp.html",
        json = f"{TRIMMED_DIR}/{{sample}}.fastp.json"
    log:
        f"{LOG_DIR}/fastp/{{sample}}.fastp.log"
    threads: 4
    wrapper:
        "file://workflow/local-wrappers/bio/fastp"

rule bowtie2_align:
    input:
        sample=[f"{TRIMMED_DIR}/{{sample}}_r1.fastq.gz", f"{TRIMMED_DIR}/{{sample}}_r2.fastq.gz"],
        idx = BOWTIE_INDEX_FILES
    output:
        bam = temp(f"{TRIMMED_DIR}/bam/{{sample}}.bam")
    log:
        f"{LOG_DIR}/bowtie2/{{sample}}.log"
    params:
        bt2_threads = 12,
        view_threads = 4,
        sort_threads = 8,
        extra= lambda wc: f"-R '@RG\\tID:{wc.sample}\\tSM:{wc.sample}'"
    threads: 24
    conda:
        "../envs/bowtie2.yaml"
    shell:
        r"""
        set -euo pipefail

        echo 'Bowtie2 version:' > {log}
        bowtie2 --version >> {log}
        echo 'Samtools version:' >> {log}
        samtools --version | head -n 1 >> {log}
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" >> {log}
        # Run alignment, piping eqverything to the log
        bowtie2 -x {BOWTIE_INDEX} -1 {input.sample[0]} -2 {input.sample[1]} --threads {params.bt2_threads} {params.extra} 2>> {log} \
        | samtools view -u -@ {params.view_threads} 2>> {log} \
        | samtools sort -@ {params.sort_threads} -o {output.bam} 2>> {log}
        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: ${{hours}}h ${{mins}}m ${{secs}}s (total ${{runtime}} seconds)" >> {log}
        """  
rule extract_unmapped_fastq:
    input:
        bam=f"{TRIMMED_DIR}/bam/{{sample}}.bam"
    output:
        r1=f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        r2=f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    log:
        f"{LOG_DIR}/bedtools/{{sample}}.log"
    threads: 60
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
        # Split threads between samtools and pigz; bedtools is single-threaded
        SAMTOOLS_THREADS=$(( {threads} * 4 / 5 ))
        PIGZ_THREADS=$(( {threads} / 5 ))
        [ $SAMTOOLS_THREADS -lt 1 ] && SAMTOOLS_THREADS=1
        [ $PIGZ_THREADS -lt 1 ] && PIGZ_THREADS=1

        echo 'Bedtools version:' > {log}
        bedtools --version >> {log}
        echo 'Samtools version:' >> {log}
        samtools --version | head -n 1 >> {log}
        echo "SAMTOOLS_THREADS: $SAMTOOLS_THREADS" >> {log}
        echo "PIGZ_THREADS: $PIGZ_THREADS" >> {log}
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" >> {log}

        TMPDIR="${{TMPDIR:-/tmp}}"
        JOB_ID="${{SLURM_JOB_ID:-manual}}"
        TMPJOB=$(mktemp -d "${{TMPDIR}}/bam2fq_${{JOB_ID}}_XXXXXX" 2>> {log} || echo "/tmp/bam2fq_${{JOB_ID}}_manualtmp")

        samtools view -u -f 12 -F 256 --threads $SAMTOOLS_THREADS {input.bam} 2>> {log} \
        | samtools sort -n --threads $SAMTOOLS_THREADS -T $TMPJOB/{wildcards.sample}_sort_tmp -O BAM - 2>> {log} \
        | bedtools bamtofastq -i - 2>> {log} \
            -fq >(pigz -p $PIGZ_THREADS --fast > {output.r1}) \
            -fq2 >(pigz -p $PIGZ_THREADS --fast > {output.r2})

        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: ${{hours}}h ${{mins}}m ${{secs}}s (total ${{runtime}} seconds)" >> {log}
        echo "Temporary directory to be removed: $TMPJOB" >> {log}
        rm -rf "$TMPJOB"
        """