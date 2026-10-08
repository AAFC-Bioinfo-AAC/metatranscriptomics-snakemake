'''
    Filename: coassembly_annotation.smk
    Original author: Katherine James-Gzyl
    Date created: 2025/08/15
    Updated: 2026/10/08
    Workflow controller version: Snakemake 9.20.0
'''

import os


# Use a supplied reference assembly when configured; otherwise assemble RNA.
reference = config.get("reference_assembly")

if reference is not None and not isinstance(reference, str):
    raise ValueError("reference_assembly must be a path or null.")

reference = reference.strip() if reference else None

if reference:
    reference = os.path.expandvars(os.path.expanduser(reference))
    COASSEMBLY_FASTA = os.path.abspath(
        os.path.join(PROJECT_ROOT, reference)
    )
else:
    COASSEMBLY_FASTA = f"{MEGAHIT_DIR}/final.contigs.fa"


# Explicitly select the Bowtie2 index format.
large_index = config.get("index_coassembly", {}).get(
    "large_index", False
)

if not isinstance(large_index, bool):
    raise ValueError(
        "index_coassembly.large_index must be true or false."
    )

index_suffix = "bt2l" if large_index else "bt2"

COASSEMBLY_INDEX_FILES = [
    f"{COASSEMBLY_INDEX}/coassembly.{part}.{index_suffix}"
    for part in ("1", "2", "3", "4", "rev.1", "rev.2")
]


# Match this setting to the library preparation protocol.
strandedness = config.get("featurecounts", {}).get(
    "strandedness", 0
)

if (
    isinstance(strandedness, bool)
    or not isinstance(strandedness, int)
    or strandedness not in (0, 1, 2)
):
    raise ValueError(
        "featurecounts.strandedness must be 0, 1 or 2."
    )


rule megahit_coassembly:
    input:
        r1 = expand(
            f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
            sample=SAMPLES
        ),
        r2 = expand(
            f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz",
            sample=SAMPLES
        )
    output:
        contigs = f"{MEGAHIT_DIR}/final.contigs.fa"
    threads:
        config.get("megahit_coassembly", {}).get("threads", 8)
    resources:
        mem_mb = 256000
    params:
        r1_files = lambda wildcards, input: [
            os.path.abspath(path) for path in input.r1
        ],
        r2_files = lambda wildcards, input: [
            os.path.abspath(path) for path in input.r2
        ],
        memory_bytes = lambda wildcards, resources: int(
            resources.mem_mb * 1024 * 1024 * 0.9
        )
    conda:
        "../envs/megahit.yaml"
    log:
        f"{LOG_DIR}/coassembly/megahit.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.contigs:q})" \
                 "$(dirname {log:q})"
        : > {log:q}

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"

        run_dir="$(mktemp -d "$tmpbase/megahit_coassembly.XXXXXX")"
        run_dir="$(cd "$run_dir" && pwd)"
        trap 'rm -rf -- "$run_dir"' EXIT

        # MEGAHIT 1.2.9 uses unquoted paths in some internal commands.
        if [[ "$run_dir" =~ [^a-zA-Z0-9_./-] ]]; then
            echo "ERROR: Use a TMPDIR path containing only letters, digits, underscores, dots, slashes or hyphens for MEGAHIT 1.2.9." >> {log:q}
            exit 1
        fi

        # Local aliases preserve source paths containing spaces or commas.
        r1_files=( {params.r1_files:q} )
        r2_files=( {params.r2_files:q} )

        if (( ${{#r1_files[@]}} == 0 || ${{#r1_files[@]}} != ${{#r2_files[@]}} )); then
            echo "ERROR: Coassembly requires matching nonempty lists of read pairs." >> {log:q}
            exit 1
        fi

        r1_list=""
        r2_list=""

        for (( i=0; i<${{#r1_files[@]}}; i++ )); do
            ln -s "${{r1_files[$i]}}" "$run_dir/r1_$i.fastq.gz"
            ln -s "${{r2_files[$i]}}" "$run_dir/r2_$i.fastq.gz"

            r1_list+="${{r1_list:+,}}$run_dir/r1_$i.fastq.gz"
            r2_list+="${{r2_list:+,}}$run_dir/r2_$i.fastq.gz"
        done

        # The directory supplied to --tmp-dir must already exist.
        mkdir -p "$run_dir/tmp"

        echo "MEGAHIT working directory: $run_dir" >> {log:q}

        megahit \
            -1 "$r1_list" \
            -2 "$r2_list" \
            -t {threads} \
            -m {params.memory_bytes} \
            -o "$run_dir/assembly" \
            --out-prefix final \
            --tmp-dir "$run_dir/tmp" \
            >> {log:q} 2>&1

        if [[ ! -s "$run_dir/assembly/final.contigs.fa" ]]; then
            echo "ERROR: MEGAHIT produced no nonempty contig file." >> {log:q}
            exit 1
        fi

        cp -- "$run_dir/assembly/final.contigs.fa" \
            {output.contigs:q}
        """


rule index_coassembly:
    input:
        coassembly = COASSEMBLY_FASTA
    output:
        index = COASSEMBLY_INDEX_FILES
    params:
        prefix = f"{COASSEMBLY_INDEX}/coassembly",
        reference = lambda wildcards, input: os.path.abspath(
            input.coassembly
        ),
        index_option = "--large-index" if large_index else ""
    threads:
        config.get("index_coassembly", {}).get("threads", 8)
    conda:
        "../envs/bowtie2.yaml"
    log:
        f"{LOG_DIR}/coassembly/coassembly_index.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {params.prefix:q})" \
                 "$(dirname {log:q})"
        : > {log:q}

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"

        run_dir="$(mktemp -d "$tmpbase/coassembly_index.XXXXXX")"
        run_dir="$(cd "$run_dir" && pwd)"
        trap 'rm -rf -- "$run_dir"' EXIT

        if [[ "$run_dir" == *,* ]]; then
            echo "ERROR: Bowtie2 indexing requires a TMPDIR path without commas." >> {log:q}
            exit 1
        fi

        ln -s {params.reference:q} "$run_dir/reference.fa"

        bowtie2-build \
            {params.index_option} \
            --threads {threads} \
            "$run_dir/reference.fa" \
            "$run_dir/coassembly" \
            >> {log:q} 2>&1

        for index_file in {output.index:q}; do
            index_name="${{index_file##*/}}"

            if [[ ! -s "$run_dir/$index_name" ]]; then
                echo "ERROR: Missing or empty Bowtie2 index: $index_file. If a large index was generated, set index_coassembly.large_index: true." >> {log:q}
                exit 1
            fi
        done

        for index_file in {output.index:q}; do
            cp -- "$run_dir/${{index_file##*/}}" "$index_file"
        done
        """


rule bowtie2_map_transcripts:
    input:
        index = COASSEMBLY_INDEX_FILES,
        r1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        r2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam",
        bai = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam.bai"
    params:
        index_files = lambda wildcards, input: [
            os.path.abspath(path) for path in input.index
        ]
    threads:
        config.get("bowtie2_map_transcripts", {}).get("threads", 16)
    conda:
        "../envs/bowtie2.yaml"
    log:
        f"{LOG_DIR}/sorted_bam/{{sample}}_sorted.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.bam:q})" \
                 "$(dirname {log:q})"
        : > {log:q}

        if (( {threads} < 2 )); then
            echo "ERROR: bowtie2_map_transcripts requires at least 2 threads." >> {log:q}
            exit 1
        fi

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"

        run_dir="$(mktemp -d "$tmpbase/coassembly_map.XXXXXX")"
        run_dir="$(cd "$run_dir" && pwd)"
        trap 'rm -rf -- "$run_dir"' EXIT

        # Expose only the selected index format to Bowtie2.
        # This avoids using stale indexes of the other format.
        for index_file in {params.index_files:q}; do
            ln -s "$index_file" "$run_dir/${{index_file##*/}}"
        done

        t_bowtie2=$(( {threads} / 2 ))
        sort_extra=$(( {threads} - t_bowtie2 - 1 ))

        # SAMtools sort accepts Bowtie2's SAM stream directly.
        bowtie2 \
            -x "$run_dir/coassembly" \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            --local \
            -p "$t_bowtie2" \
            2>> {log:q} \
        | samtools sort \
            -@ "$sort_extra" \
            -T "$run_dir/sort" \
            -O BAM \
            -o {output.bam:q} \
            - \
            2>> {log:q}

        samtools quickcheck {output.bam:q} >> {log:q} 2>&1

        samtools index \
            -@ $(( {threads} - 1 )) \
            {output.bam:q} \
            {output.bai:q} \
            >> {log:q} 2>&1
        """


rule assembly_stats_depth:
    input:
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam",
        bai = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam.bai"
    output:
        stats = f"{ASSEMBLY_MAPPING}/{{sample}}.flagstat.txt",
        depth = f"{ASSEMBLY_MAPPING}/{{sample}}.coverage.txt.gz",
        idxstats = f"{ASSEMBLY_MAPPING}/{{sample}}.idxstats.txt.gz"
    threads:
        config.get("assembly_stats_depth", {}).get("threads", 2)
    conda:
        "../envs/bedtools.yaml"
    log:
        f"{LOG_DIR}/coassembly/{{sample}}_stats_depth.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.stats:q})" \
                 "$(dirname {log:q})"
        : > {log:q}

        if (( {threads} < 2 )); then
            echo "ERROR: assembly_stats_depth requires at least 2 threads." >> {log:q}
            exit 1
        fi

        t_pigz=$(( {threads} - 1 ))

        samtools flagstat {input.bam:q} \
            > {output.stats:q} 2>> {log:q}

        # Preserve default depth behavior: omit zero-depth positions and
        # count both mates where their aligned sequences overlap.
        samtools depth {input.bam:q} 2>> {log:q} \
            | pigz -p "$t_pigz" \
            > {output.depth:q} 2>> {log:q}

        samtools idxstats {input.bam:q} 2>> {log:q} \
            | pigz -p "$t_pigz" \
            > {output.idxstats:q} 2>> {log:q}
        """


rule prodigal_genes:
    input:
        coassembly = COASSEMBLY_FASTA
    output:
        proteins = f"{PRODIGAL_DIR}/coassembly.faa",
        nucs = f"{PRODIGAL_DIR}/coassembly.fna",
        gff = f"{PRODIGAL_DIR}/coassembly.gff",
        saf = f"{PRODIGAL_DIR}/coassembly.saf"
    threads: 1
    conda:
        "../envs/prodigal.yaml"
    log:
        f"{LOG_DIR}/prodigal/coassembly_prodigal.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.gff:q})" \
                 "$(dirname {log:q})"

        prodigal \
            -i {input.coassembly:q} \
            -a {output.proteins:q} \
            -d {output.nucs:q} \
            -o {output.gff:q} \
            -p meta \
            -f gff \
            > {log:q} 2>&1

        # Preserve GFF ID attributes for featureCounts and CAZyme joins.
        awk -F '\t' '
            BEGIN {{
                OFS="\t"
                print "GeneID", "Chr", "Start", "End", "Strand"
            }}
            /^#/ {{ next }}
            $3=="CDS" {{
                id=""
                id_count=0
                split($9, attributes, ";")

                for (i in attributes) {{
                    if (attributes[i] ~ /^ID=/) {{
                        id=substr(attributes[i], 4)
                        id_count++
                    }}
                }}

                if (NF!=9 || id_count!=1 || id=="" || (id in seen) ||
                    $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ ||
                    $4<1 || $5<$4 || ($7!="+" && $7!="-")) {{
                    print "ERROR: Invalid or duplicate CDS annotation at GFF line " NR > "/dev/stderr"
                    bad=1
                    exit 1
                }}

                seen[id]=1
                count++
                print id, $1, $4, $5, $7
            }}
            END {{
                if (!bad && count==0) {{
                    print "ERROR: No CDS features were predicted; gene counting cannot proceed." > "/dev/stderr"
                    exit 1
                }}
            }}
        ' {output.gff:q} > {output.saf:q} 2>> {log:q}

        for result_file in {output:q}; do
            if [[ ! -s "$result_file" ]]; then
                echo "ERROR: Missing or empty Prodigal output: $result_file" >> {log:q}
                exit 1
            fi
        done
        """


rule featurecounts:
    input:
        saf = f"{PRODIGAL_DIR}/coassembly.saf",
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.coassembly.sorted.bam"
    output:
        counts = f"{FEATURECOUNTS_DIR}/{{sample}}_counts.txt",
        summary = f"{FEATURECOUNTS_DIR}/{{sample}}_counts.txt.summary"
    params:
        strandedness = strandedness
    threads:
        config.get("featurecounts", {}).get("threads", 4)
    conda:
        "../envs/featurecounts.yaml"
    log:
        f"{LOG_DIR}/featurecounts/{{sample}}_featurecounts.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.counts:q})" \
                 "$(dirname {log:q})"

        featureCounts \
            -a {input.saf:q} \
            -F SAF \
            -p \
            --countReadPairs \
            -s {params.strandedness} \
            -T {threads} \
            -o {output.counts:q} \
            {input.bam:q} \
            > {log:q} 2>&1

        for result_file in {output:q}; do
            if [[ ! -s "$result_file" ]]; then
                echo "ERROR: Missing or empty featureCounts output: $result_file" >> {log:q}
                exit 1
            fi
        done
        """
