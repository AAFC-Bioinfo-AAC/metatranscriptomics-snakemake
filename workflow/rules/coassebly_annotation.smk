
'''
    Filename: coassembly_annotation.smk
    Author: Katherine James-Gzyl
    Date created: 2025/08/15
    Snakemake version: 9.6.0
'''
rule megahit_coassembly:
    input: 
        r1 = expand(f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz", sample=SAMPLES),
        r2 = expand(f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz", sample=SAMPLES)
    output:
        f"{MEGAHIT_DIR}/final.contigs.fa"
    threads: config.get("megahit_coassembly", {}).get("threads", 8)
    conda:
        "../envs/megahit.yaml"
    log:
        f"{LOG_DIR}/coassembly/megahit.log"
    shell:
        r"""
        set -euo pipefail

        r1_list=$(echo {input.r1} | tr ' ' ',')
        r2_list=$(echo {input.r2} | tr ' ' ',')
        
        export TMPDIR="${{TMPDIR:-/tmp}}/megahit_coassembly_${{RANDOM}}_$(date +%s)"
        mkdir -p "$TMPDIR" || (echo "Failed to create TMPDIR $TMPDIR" >> {log}; exit 1)
        echo "Using TMPDIR: $TMPDIR" >> {log}

        megahit \
          -1 "$r1_list" \
          -2 "$r2_list" \
          -t {threads} \
          -o {MEGAHIT_DIR} --force\
          --out-prefix final \
          --tmp-dir "$TMPDIR" \
          > {log} 2>&1

        echo "Removing tmpdir $TMPDIR" >> {log}
        rm -rf "$TMPDIR"
        """ 
rule index_coassembly:
    input:
        coassembly = f"{MEGAHIT_DIR}/final.contigs.fa"
    output:
        index = [
            f"{COASSEMBLY_INDEX}/coassembly.1.bt2",
            f"{COASSEMBLY_INDEX}/coassembly.2.bt2",
            f"{COASSEMBLY_INDEX}/coassembly.3.bt2",
            f"{COASSEMBLY_INDEX}/coassembly.4.bt2",
            f"{COASSEMBLY_INDEX}/coassembly.rev.1.bt2",
            f"{COASSEMBLY_INDEX}/coassembly.rev.2.bt2"
        ]
    log:
        f"{LOG_DIR}/coassembly/coassembly_index.log"
    threads: config.get("index_coassembly", {}).get("threads", 8)
    conda:
        "../envs/bowtie2.yaml"  
    params:
        prefix = f"{COASSEMBLY_INDEX}/coassembly"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {COASSEMBLY_INDEX} 

        bowtie2-build --threads {threads} {input.coassembly} {params.prefix} &>> {log}
        """
rule bowtie2_map_transcripts: 
    input:
        index = f"{COASSEMBLY_INDEX}/coassembly.1.bt2",
        r1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        r2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam"
    log:
        f"{LOG_DIR}/sorted_bam/{{sample}}_sorted.log"
    threads: config.get("bowtie2_map_transcripts", {}).get("threads", 16)
    conda:
        "../envs/bowtie2.yaml"
    params:
         index_prefix = f"{COASSEMBLY_INDEX}/coassembly"
    shell:
        r"""
        set -euo pipefail

        #divide threads
        t_bowtie2=$(( ({threads} + 1) / 2 ))
        t_sort=$(( {threads} - t_bowtie2 ))

        bowtie2 -x {params.index_prefix} -1 {input.r1} -2 {input.r2} --local -p {t_bowtie2} 2>> {log}\
        | samtools view -bS - 2>> {log} \
        | samtools sort -@ {t_sort} -o {output.bam} 2>> {log}
        """
rule assembly_stats_depth:
    input:
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam"
    output:
        stats = f"{ASSEMBLY_MAPPING}/{{sample}}.flagstat.txt",
        depth = f"{ASSEMBLY_MAPPING}/{{sample}}.coverage.txt.gz",
        idxstats = f"{ASSEMBLY_MAPPING}/{{sample}}.idxstats.txt.gz"
    threads: config.get("assembly_stats_depth", {}).get("threads", 2)
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} &&
        samtools depth {input.bam} | pigz -p {threads} > {output.depth} &&
        samtools idxstats {input.bam} | pigz -p {threads} > {output.idxstats}
        """
rule prodigal_genes:
    input:
        coassembly = f"{MEGAHIT_DIR}/final.contigs.fa"
    output:
        proteins = f"{PRODIGAL_DIR}/coassembly.faa",
        nucs = f"{PRODIGAL_DIR}/coassembly.fna",
        gff = f"{PRODIGAL_DIR}/coassembly.gff",
        saf = f"{PRODIGAL_DIR}/coassembly.saf"
    threads: 1
    log:
        f"{LOG_DIR}/prodigal/coassembly_prodigal.log"
    conda:
        "../envs/prodigal.yaml"
    shell:
        r"""
        set -euo pipefail

        # Ensure output directory exists
        mkdir -p {PRODIGAL_DIR}

        prodigal -i {input.coassembly} \
            -a {output.proteins} \
            -d {output.nucs} \
            -o {output.gff} \
            -p meta \
            -f gff 2>&1 | grep -v "^Finding genes in sequence" | grep -v "^Request:" >> {log}

        awk '$3=="CDS" {{OFS="\t"; split($9,a,";"); for(i in a){{if(a[i]~/^ID=/)id=substr(a[i],4)}}; print id, $1, $4, $5, $7}}' {output.gff} | \
        awk 'BEGIN{{print "GeneID\tChr\tStart\tEnd\tStrand"}} 1' > {output.saf}
        """
rule featurecounts:
    input:
        saf = f"{PRODIGAL_DIR}/coassembly.saf",
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam"
    output:
        counts = f"{FEATURECOUNTS_DIR}/{{sample}}_counts.txt"
    threads: config.get("featurecounts", {}).get("threads", 4)
    log:
        f"{LOG_DIR}/featurecounts/{{sample}}_featurecounts.log"
    conda:
        "../envs/featurecounts.yaml"
    shell:
        r"""
        set -euo pipefail

        # Ensure output directory exists
        mkdir -p {FEATURECOUNTS_DIR}

        featureCounts -a {input.saf} -F SAF -p -T {threads} -o {output.counts} {input.bam} &>> {log}
        """