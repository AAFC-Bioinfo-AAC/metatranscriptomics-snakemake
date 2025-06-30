configfile: "config/config.yaml"
## REMEMBER TO REMOVE TEST SUBSET
# For test rules:
SAMPLES_SUBSET = config["sample_subset"]

from dotenv import load_dotenv
import os
import sys
import glob 

TMPDIR="/gpfs/fs7/aafc/scratch/kjg000/tmpdir"
print(f"DEBUG: Using TMPDIR = {TMPDIR}")

## Load .env 
load_dotenv('.env')
#TMPDIR = os.getenv("TMPDIR", "/gpfs/fs7/aafc/scratch/kjg000/tmpdir")
#print(f"DEBUG: Using TMPDIR = {TMPDIR}")
RGI_CARD = os.getenv("RGI_CARD") or config.get("card_latest") #fallback if user do not have common card database
if not RGI_CARD:
    raise ValueError("You must set the environment variable RGI_CARD to your CARD DB directory")

## Early exit if RGI CARD database cannot be loaded.
if not RGI_CARD or not os.path.exists(RGI_CARD):
    sys.exit(f"\nERROR: RGI_CARD not set or path not found: {RGI_CARD}\nEdit .env or config.yaml.\n")

## Variable set up and sample loading
with open(config["samplelist"]) as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

READS_DIR = config["reads_dir"]
TRIMMED_DIR = config["reads_trimmed"]
HOST_DEP_DIR = config["reads_host_dep"]
RRNA_DEP_DIR = config["reads_rRNA_dep"]
ASSEMBLIES_DIR =config["assemblies_dir"]
KRAKEN_OUTPUT_DIR = config["taxonomy_short_reads_dir"]
CARD_RGI_OUTPUT_DIR = config["amr_screaning_dir"]
RNAQUAST_DIR = config["rnaquast_dir"]
LOG_DIR = config["log_files"]
BOWTIE_INDEX = config["bowtie2_index"]
BOWTIE_INDEX_FILES = [
    f"{BOWTIE_INDEX}.1.bt2",
    f"{BOWTIE_INDEX}.2.bt2",
    f"{BOWTIE_INDEX}.3.bt2",
    f"{BOWTIE_INDEX}.4.bt2",
    f"{BOWTIE_INDEX}.rev.1.bt2",
    f"{BOWTIE_INDEX}.rev.2.bt2"
]
RRNA_DB = config["sortmerna_DB"]
TAXONOMY_DB = config["gtbd_DB"]
BUSCO_LINEAGES = config["busco_lineages"]
LINEAGES = BUSCO_LINEAGES.keys()
ASSEMBLY_INDEX = config["bowtie2_assembly_index"]
ASSEMBLY_MAPPING = config["mapping_back_dir"]

## Read samples 
def find_read_file(sample, readnum):
    patterns = [
        f"{READS_DIR}/{sample}_R{readnum}.fastq.gz",
        f"{READS_DIR}/{sample}_r{readnum}.fastq.gz"
    ]
    for pattern in patterns:
        found = glob.glob(pattern)
        if found:
            return found[0]
    raise ValueError(f"No FASTQ found for sample {sample}, read {readnum}.")

# Optional: check specifically for card_reference.fasta
def has_fasta(path):
    return (
        os.path.exists(os.path.join(path, "card_reference.fasta")) or
        any(name.endswith('.fasta') for name in os.listdir(path))
    )

if not (os.path.exists(f"{RGI_CARD}/card.json")
        and has_fasta(RGI_CARD)):
    sys.exit(
        f"\nERROR: CARD DB not properly prepared at {RGI_CARD} (missing card.json or .fasta).\n"
        "Please follow manual setup instructions (see pipeline README).\n"
    )

# Set to the last step in the pipeline
rule all:
    input:
        f"{RNAQUAST_DIR}/LLC82Sep06GR_archaea"

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
    threads: 10
    wrapper:
        "file://snakemake-wrappers/bio/fastp"

rule bowtie2_align:
    input:
        sample=[f"{TRIMMED_DIR}/{{sample}}_r1.fastq.gz", f"{TRIMMED_DIR}/{{sample}}_r2.fastq.gz"],
        idx = BOWTIE_INDEX_FILES
    output:
        bam = temp("results/{sample}.bam")
    log:
        f"{LOG_DIR}/bowtie2/{{sample}}.log"
    params:
        bt2_threads = 44,
        view_threads = 4,
        sort_threads = 12,
        extra= lambda wc: f"-R '@RG\\tID:{wc.sample}\\tSM:{wc.sample}'"
    threads: 60
    conda:
        "envs/bowtie2.yaml"
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
        bam="results/{sample}.bam"
    output:
        r1=f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R1.fastq.gz",
        r2=f"{HOST_DEP_DIR}/{{sample}}_trimmed_clean_R2.fastq.gz"
    log:
        f"{LOG_DIR}/bedtools/{{sample}}.log"
    threads: 60
    conda:
        "envs/bedtools.yaml"
    shell:
        r"""
        echo 'Bedtools version:' > {log}
        bedtools --version >> {log}
        echo 'Samtools version:' >> {log}
        samtools --version | head -n 1 >> {log}
        echo 'Wall time:' >> {log}
        start_time=$(date '+%s')
        # Extract unmapped reads (-f 12), not secondary alignments (-F 256), sort by read name
        samtools view -u -f 12 -F 256 --threads {threads} {input.bam} 2>> {log} \
        | samtools sort -n --threads {threads} -O BAM - 2>> {log} \
        | bedtools bamtofastq -i - 2>> {log} \
            -fq >(pigz --fast > {output.r1}) \
            -fq2 >(pigz --fast > {output.r2})
        end_time=$(date '+%s')
        runtime=$((end_time - start_time))
        echo "Elapsed seconds: $runtime" >> {log}
        """
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
    threads: 60
    conda:
        "envs/sortmerna.yaml"
    shell:
        r"""
        set -euo pipefail

        # Wall time tracking
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" > {log}

        export TMPDIR={TMPDIR}
        WORKDIR="$TMPDIR/smrn_{wildcards.sample}_$RANDOM"
        mkdir -p "$WORKDIR"

        echo "Temporary workdir: $WORKDIR" >> {log}
        echo "Using SortMeRNA version:" >> {log}
        sortmerna --version >> {log}

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

        # Wall time logging
        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: $hours"h" $mins"m" $secs"s" (total $runtime seconds)" >> {log}
        echo "Temporary directory to be removed: $WORKDIR" >> {log}
        rm -rf "$WORKDIR"
        """

rule rna_spades:
    input:
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        fasta = f"{ASSEMBLIES_DIR}/{{sample}}.fasta"
    log:
        f"{LOG_DIR}/spades/{{sample}}.log"
    conda:
        "envs/RNAspades.yaml"
    threads: 60
    params:
        memory = "64000"
    shell:
        r"""
        set -euo pipefail

        export TMPDIR={TMPDIR}
        mkdir -p "$TMPDIR"

        echo 'SPAdes version:' > {log}
        spades.py --version >> {log}
        echo 'Wall time:' >> {log}
        start_time=$(date '+%s')

        outdir=$(mktemp -d "$TMPDIR/rnaspades_{wildcards.sample}_XXXXXX")
        spades_tmp="$outdir/tmp"
        mkdir -p "$spades_tmp"

        spades.py --rna \
            -t {threads} \
            -m {params.memory} \
            --tmp-dir "$spades_tmp" \
            -1 {input.R1} \
            -2 {input.R2} \
            -o "$outdir" \
            &>> {log}

        # Try to copy the rnaSPAdes output
        if [ -f "$outdir/transcripts.fasta" ]; then
            cp "$outdir/transcripts.fasta" {output.fasta}
        else
            echo "ERROR: Expected transcripts.fasta not found in $outdir!" >&2
            ls -lh "$outdir"
            exit 1
        fi

        end_time=$(date '+%s')
        runtime=$((end_time - start_time))
        echo "Elapsed seconds: $runtime" >> {log}
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
    threads: 4 
    shell:
        r"""
        set -euo pipefail

        export TMPDIR={TMPDIR}/rnaquast_{wildcards.sample}_{wildcards.lineage}_$RANDOM
        mkdir -p "$TMPDIR"
        echo "TMPDIR is $TMPDIR" >> {log}

        echo 'Wall time:' >> {log}
        start_time=$(date '+%s')

        rnaQUAST.py --transcripts {input.fasta} \
            --output {output.report_dir} \
            --threads {threads} \
            --busco {input.busco_lineage} \
            &>> {log}
            
        end_time=$(date '+%s')
        runtime=$((end_time - start_time))
        echo "Elapsed seconds: $runtime" >> {log}
        echo "Temporary files can be found at $TMPDIR" >> {log}
        # rm -rf "$TMPDIR"  # commented-out or removed, so tmpdir is kept.
        """
rule bowtie_build_transcript:
    input:
        fasta = f"{ASSEMBLIES_DIR}/{{sample}}.fasta"
    output:
        index = [
            f"{ASSEMBLY_INDEX}/{{sample}}.index.1.bt2",
            f"{ASSEMBLY_INDEX}/{{sample}}.index.2.bt2",
            f"{ASSEMBLY_INDEX}/{{sample}}.index.3.bt2",
            f"{ASSEMBLY_INDEX}/{{sample}}.index.4.bt2",
            f"{ASSEMBLY_INDEX}/{{sample}}.index.rev.1.bt2",
            f"{ASSEMBLY_INDEX}/{{sample}}.index.rev.2.bt2"
        ]
    log:
        f"{LOG_DIR}/mapping_back/{{sample}}_index.log"
    threads: 8  
    conda:
        "envs/bowtie2.yaml"  
    params:
        prefix = f"{ASSEMBLY_INDEX}/{{sample}}.index"
    shell:
        r"""
        set -euo pipefail

        echo 'Bowtie2 version:' > {log}
        bowtie2 --version >> {log}
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" >> {log}

        bowtie2-build --threads {threads} {input.fasta} {params.prefix} &> {log}
        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: ${{hours}}h ${{mins}}m ${{secs}}s (total ${{runtime}} seconds)" >> {log}
        """

rule bowtie2_map_transcripts: 
    input:
        index = f"{ASSEMBLY_INDEX}/{{sample}}.index.1.bt2", 
        r1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        r2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.sorted.bam"
    log:
        f"{LOG_DIR}/mapping_back/{{sample}}_sorted.log"
    threads: 40
    conda:
        "envs/bowtie2.yaml"
    params:
         index_prefix = f"{ASSEMBLY_INDEX}/{{sample}}.index"
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

        bowtie2 -x {params.index_prefix} -1 {input.r1} -2 {input.r2} --local -p {threads} 2>> {log}\
        | samtools view -bS - 2>> {log} \
        | samtools sort -@ {threads} -o {output.bam} 2>> {log}
        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: ${{hours}}h ${{mins}}m ${{secs}}s (total ${{runtime}} seconds)" >> {log}
        """
rule assembly_stats_depth:
    input:
        bam = f"{ASSEMBLY_MAPPING}/{{sample}}.sorted.bam"
    output:
        stats = f"{ASSEMBLY_MAPPING}/{{sample}}.flagstat.txt",
        depth = f"{ASSEMBLY_MAPPING}/{{sample}}.coverage.txt.gz",
        idxstats = f"{ASSEMBLY_MAPPING}/{{sample}}.idxstats.txt.gz"
    threads: 2
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} &&
        samtools depth {input.bam} | pigz -p {threads} > {output.depth} &&
        samtools idxstats {input.bam} | pigz -p {threads} > {output.idxstats}
        """

rule kraken2: 
    input:
        ref = f"{TAXONOMY_DB}",
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
        
    output:
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt",
        kraken = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.kraken" 
    log:
        f"{LOG_DIR}/kraken2/{{sample}}.log"
    conda:
        "envs/kraken2.yaml"
    threads: 2
    params:
        conf_threshold = "0.5"
    shell:
        r"""
        set -euo pipefail

        # Wall time tracking
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" > {log}

        echo 'Kraken2 version:' >> {log}
        kraken2 --version >> {log}
        # echo 'Database is GTDB 220' >> {log} #uncomment if there is a text file with the version
        # cat {input.ref}/VERSION.txt >> {log} #uncomment if there is a text file with the version 

        kraken2 --use-names \
            --threads {threads} \
            --db {input.ref} \
            --confidence {params.conf_threshold} \
            --report-zero-counts \
            --paired {input.R1} {input.R2} \
            --report {output.report} \
            --output {output.kraken} \
            &>> {log}
            
        # Wall time logging
        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: $hours"h" $mins"m" $secs"s" (total $runtime seconds)" >> {log}
        """

rule bracken:
    input:
        ref = f"{TAXONOMY_DB}",
        report = f"{KRAKEN_OUTPUT_DIR}/{{sample}}.report.txt"
    output:
        species = f"{KRAKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt",
        genus = f"{KRAKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt",
        phylum = f"{KRAKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt"
    log:f"{LOG_DIR}/bracken/{{sample}}.log"
    conda:
        "envs/kraken2.yaml"
    threads: 10
    params:
        readlen = 150
    shell:
        r"""
        set -euo pipefail

        # Wall time tracking
        start_time=$(date +%s)
        start_hr=$(date)
        echo "Started at: $start_hr" > {log}

        echo 'bracken version:' >> {log}
        bracken --version >> {log}

        mkdir -p $(dirname {output.species}) $(dirname {output.genus}) $(dirname {output.phylum})

        bracken -r {params.readlen} -t {threads} -d {input.ref} -i {input.report} -l S -o {output.species} &>> {log}

        bracken -r {params.readlen} -t {threads} -d {input.ref} -i {input.report} -l G -o {output.genus} &>> {log}

        bracken -r {params.readlen} -t {threads} -d {input.ref} -i {input.report} -l P -o {output.phylum} &>> {log}

        # Wall time logging
        end_time=$(date +%s)
        end_hr=$(date)
        runtime=$((end_time - start_time))
        hours=$((runtime / 3600))
        mins=$(((runtime % 3600) / 60))
        secs=$((runtime % 60))
        echo "Finished at: $end_hr" >> {log}
        echo "Wall time: $hours"h" $mins"m" $secs"s" (total $runtime seconds)" >> {log}
        """
rule combine_bracken_outputs:
    input:
        species = expand(f"{KRAKEN_OUTPUT_DIR}/species/{{sample}}_bracken.species.report.txt", sample=SAMPLES_SUBSET),
        genus = expand(f"{KRAKEN_OUTPUT_DIR}/genus/{{sample}}_bracken.genus.report.txt", sample=SAMPLES_SUBSET),
        phylum = expand(f"{KRAKEN_OUTPUT_DIR}/phylum/{{sample}}_bracken.phylum.report.txt", sample=SAMPLES_SUBSET)
    output:
        species = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_species.txt",
        genus = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_genus.txt",
        phylum = f"{KRAKEN_OUTPUT_DIR}/merged_abundance_phylum.txt"
    log: f"{LOG_DIR}/bracken/combine_bracken_outputs.log"
    conda:
        "envs/kraken2.yaml"
    shell:
        r"""
        set -euo pipefail

        echo "Combining species bracken outputs..." > {log}
        combine_bracken_outputs.py --files {input.species} --output {output.species} &>> {log}

        echo "Combining genus bracken outputs..." >> {log}
        combine_bracken_outputs.py --files {input.genus} --output {output.genus} &>> {log}

        echo "Combining phylum bracken outputs..." >> {log}
        combine_bracken_outputs.py --files {input.phylum} --output {output.phylum} &>> {log}
        """
rule bracken_extract:
    input:
        species_table = f"{KRAKEN_OUTPUT_DIR}/testing/merged_abundance_species.txt",
        genus_table = f"{KRAKEN_OUTPUT_DIR}/testing/merged_abundance_genus.txt",
        phylum_table = f"{KRAKEN_OUTPUT_DIR}/testing/merged_abundance_phylum.txt"
    output:
        species_raw = f"{KRAKEN_OUTPUT_DIR}/Bracken_species_raw_abundance.csv",
        species_rel = f"{KRAKEN_OUTPUT_DIR}/Bracken_species_relative_abundance.csv", 
        genus_raw = f"{KRAKEN_OUTPUT_DIR}/Bracken_genus_raw_abundance.csv",
        genus_rel = f"{KRAKEN_OUTPUT_DIR}/Bracken_genus_relative_abundance.csv",
        phylum_raw = f"{KRAKEN_OUTPUT_DIR}/Bracken_phylum_raw_abundance.csv",
        phylum_rel = f"{KRAKEN_OUTPUT_DIR}/Bracken_phylum_relative_abundance.csv"
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        scripts/extract_bracken_columns.py --input {input.species_table} --level species --out-raw {output.species_raw} --out-rel {output.species_rel}
        scripts/extract_bracken_columns.py --input {input.genus_table} --level genus --out-raw {output.genus_raw} --out-rel {output.genus_rel}
        scripts/extract_bracken_columns.py --input {input.phylum_table} --level phylum --out-raw {output.phylum_raw} --out-rel {output.phylum_rel}

        """

rule rgi_reload_database:
    input:
        card_json = f"{RGI_CARD}/card.json",
        card_fasta = f"{RGI_CARD}/card_reference.fasta"
    output:
        f"{LOG_DIR}/rgi_reload_db.done"
    log:
        f"{LOG_DIR}/rgi/rgi_reload_db.log"
    conda:
        "envs/rgi.yaml"
    shell:
        r"""
        set -e

        LOCALDB="{workflow.basedir}/localDB"

        # Standard "is file, not dir" check
        if [ -f "$LOCALDB" ]; then
            echo "ERROR: localDB exists as a file, but should be a directory." >&2
            exit 1
        fi
        mkdir -p "$LOCALDB"

        # Debug output
        echo "Current files in $LOCALDB before loading:" >> {log}
        ls -lh "$LOCALDB" >> {log}

        # Only reload if DB is not present (or is empty)
        if [ -f "$LOCALDB/card.json" ] && [ -s "$LOCALDB/card.json" ] && \
           [ -f "$LOCALDB/card_reference.fasta" ] && [ -s "$LOCALDB/card_reference.fasta" ]; then
            echo "CARD local database already present in $LOCALDB, skipping reload." >> {log}
        else
            echo "Reloading CARD into $LOCALDB..." >> {log}
            rgi load --card_json {input.card_json} --card_annotation {input.cexitard_fasta} --local --local_database "$LOCALDB" >> {log} 2>&1
        fi

        touch {output}
        """

rule rgi_bwt:
    input:
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz",
        reload_marker = f"{LOG_DIR}/rgi_reload_db.done"
    output:
        json = f"{CARD_RGI_OUTPUT_DIR}/{{sample}}_paired.allele_mapping_data.json",
        allele = f"{CARD_RGI_OUTPUT_DIR}/{{sample}}_paired.allele_mapping_data.txt"
    params:
        outprefix = lambda wc: f"{CARD_RGI_OUTPUT_DIR}/{wc.sample}_paired"
    log:
        f"{LOG_DIR}/rgi/bwt_{{sample}}.log"
    threads: 40
    conda:
        "envs/rgi.yaml"
    shell:
        r"""
        set -o pipefail
        export TMPDIR={TMPDIR}/bwt_{wildcards.sample}_$RANDOM
        mkdir -p "$TMPDIR"
        echo "Using TMPDIR: $TMPDIR" >> {log}
        export RGI_DATA_PATH={RGI_CARD}
        echo "RGI_DATA_PATH: {RGI_CARD}" >> {log}
        start_time=$(date +%s)
        echo "Started at: $(date)" >> {log}
        echo 'RGI version:' >> {log}
        rgi main --version >> {log}
        echo 'bwt without wildcard' >> {log}
        mkdir -p {CARD_RGI_OUTPUT_DIR}
        (rgi bwt \
            -1 {input.R1} \
            -2 {input.R2} \
            -a kma \
            -n {threads} \
            -o {params.outprefix} \
            --local \
            --clean \
            ) 2>&1 | egrep -av '\[W::sam_parse1\]|mapped query cannot have zero coordinate' | awk 'NF' >> {log}

        rgistatus=${{PIPESTATUS[0]}}
        if [ $rgistatus -ne 0 ]; then
        echo "ERROR: rgi bwt failed with exit code $rgistatus" >> {log}
        exit $rgistatus
        fi

        end_time=$(date +%s)
        runtime=$((end_time - start_time))

        echo "Finished at: $(date)" >> {log}
        echo "Wall time: $((runtime/3600))h $(((runtime%3600)/60))m $((runtime%60))s (total $runtime s)" >> {log}

        echo "Cleaning up temp dir: $TMPDIR" >> {log}
        rm -rf "$TMPDIR"
        echo "Cleanup completed for: $TMPDIR" >> {log}
        """


