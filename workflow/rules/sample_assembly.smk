'''
    Filename: sample_assembly.smk
    Author: Katherine James-Gzyl
    Date created: 2025/08/13
    Updated: 2026/10/09
    Snakemake version: 9.20.0
'''

import os


if "memory" in config.get("rna_spades", {}):
    raise ValueError(
        "Replace rna_spades.memory with memory_gb; SPAdes -m uses GB."
    )

RNA_SPADES_MEMORY_GB = config.get(
    "rna_spades", {}
).get("memory_gb", 60)

if (
    isinstance(RNA_SPADES_MEMORY_GB, bool)
    or not isinstance(RNA_SPADES_MEMORY_GB, int)
    or RNA_SPADES_MEMORY_GB < 1
):
    raise ValueError(
        "rna_spades.memory_gb must be a positive integer."
    )


rule rna_spades:
    input:
        R1 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R1.fastq.gz",
        R2 = f"{RRNA_DEP_DIR}/{{sample}}_rRNAdep_R2.fastq.gz"
    output:
        fasta = f"{ASSEMBLIES_DIR}/{{sample}}.fasta"
    log:
        f"{LOG_DIR}/spades/{{sample}}.log"
    conda:
        "../envs/RNAspades.yaml"
    threads: config.get("rna_spades", {}).get("threads", 24)
    resources:
        mem_mb = 64000
    params:
        memory_gb = RNA_SPADES_MEMORY_GB,
        R1_path = lambda wc, input: os.path.abspath(input.R1),
        R2_path = lambda wc, input: os.path.abspath(input.R2)
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.fasta:q})"
        : > {log:q}

        if (( {params.memory_gb} * 1024 > {resources.mem_mb} )); then
            echo "ERROR: rna_spades.memory_gb exceeds the allocated mem_mb; increase the job's memory allocation or reduce memory_gb." >> {log:q}
            exit 1
        fi

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"
        run_dir=$(mktemp -d "$tmpbase/rnaspades.XXXXXX")
        run_dir=$(cd "$run_dir" && pwd)
        trap 'rm -rf -- "$run_dir"' EXIT

        # Simple working paths avoid problems in tool-generated commands.
        if [[ "$run_dir" =~ [^a-zA-Z0-9_./-] ]]; then
            echo "ERROR: Use a TMPDIR path containing only letters, digits, underscores, dots, slashes or hyphens." >> {log:q}
            exit 1
        fi

        ln -s {params.R1_path:q} "$run_dir/R1.fastq.gz"
        ln -s {params.R2_path:q} "$run_dir/R2.fastq.gz"
        mkdir -p "$run_dir/tmp"

        echo "rnaSPAdes working directory: $run_dir" >> {log:q}

        # A nonzero exit must fail the rule, even if a partial fasta exists.
        spades.py --rna \
            -t {threads} \
            -m {params.memory_gb} \
            --tmp-dir "$run_dir/tmp" \
            -1 "$run_dir/R1.fastq.gz" \
            -2 "$run_dir/R2.fastq.gz" \
            -o "$run_dir/assembly" >> {log:q} 2>&1

        if [[ ! -s "$run_dir/assembly/transcripts.fasta" ]]; then
            echo "ERROR: rnaSPAdes completed without a nonempty transcripts.fasta." >> {log:q}
            exit 1
        fi

        mv -- "$run_dir/assembly/transcripts.fasta" {output.fasta:q}

        echo "rnaSPAdes completed; transcripts written." >> {log:q}
        """


rule rnaquast_busco:
    input:
        fasta = f"{ASSEMBLIES_DIR}/{{sample}}.fasta",
        busco_lineage = lambda wc: BUSCO_LINEAGES[wc.lineage],
        lineage_config = lambda wc: os.path.join(
            BUSCO_LINEAGES[wc.lineage], "dataset.cfg"
        )
    output:
        report_dir = directory(
            f"{RNAQUAST_DIR}/{{sample}}_{{lineage}}"
        )
    log:
        f"{LOG_DIR}/rnaquast/{{sample}}_{{lineage}}.log"
    conda:
        "../envs/rnaquast.yaml"
    threads: config.get("rnaquast_busco", {}).get("threads", 4)
    params:
        fasta_path = lambda wc, input: os.path.abspath(input.fasta),
        lineage_path = lambda wc, input: os.path.abspath(
            input.busco_lineage
        ),
        lineage_name = lambda wc, input: os.path.basename(
            os.path.normpath(input.busco_lineage)
        ),
        label = lambda wc: wc.sample
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {log:q})" "$(dirname {output.report_dir:q})"
        : > {log:q}

        if [[ ! -s {input.fasta:q} ]]; then
            echo "ERROR: Assembly is empty; regenerate the assembly before running rnaQUAST." >> {log:q}
            exit 1
        fi

        label={params.label:q}
        if [[ "$label" =~ [^a-zA-Z0-9_.-] ]]; then
            echo "ERROR: Sample IDs used by rnaQUAST must contain only letters, digits, underscores, dots or hyphens." >> {log:q}
            exit 1
        fi

        lineage_name={params.lineage_name:q}
        if [[ "$lineage_name" =~ [^a-zA-Z0-9_.-] ]]; then
            echo "ERROR: The BUSCO dataset directory name must contain only letters, digits, underscores, dots or hyphens." >> {log:q}
            exit 1
        fi

        tmpbase="${{TMPDIR:-/tmp}}"
        mkdir -p "$tmpbase"
        run_dir=$(mktemp -d "$tmpbase/rnaquast.XXXXXX")
        run_dir=$(cd "$run_dir" && pwd)
        trap 'rm -rf -- "$run_dir"' EXIT

        # rnaQUAST builds the BUSCO command without quoting its paths.
        if [[ "$run_dir" =~ [^a-zA-Z0-9_./-] ]]; then
            echo "ERROR: Use a TMPDIR path containing only letters, digits, underscores, dots, slashes or hyphens." >> {log:q}
            exit 1
        fi

        mkdir -p "$run_dir/lineages"

        ln -s {params.fasta_path:q} "$run_dir/transcripts.fasta"
        ln -s {params.lineage_path:q} "$run_dir/lineages/$lineage_name"

        printf 'Assembly: %s\nBUSCO lineage: %s\n' \
            {params.fasta_path:q} {params.lineage_path:q} >> {log:q}

        # --debug retains native BUSCO summaries that rnaQUAST otherwise deletes.
        rnaQUAST.py \
            --transcripts "$run_dir/transcripts.fasta" \
            --output_dir "$run_dir/results" \
            --threads {threads} \
            --labels {params.label:q} \
            --prokaryote \
            --busco "$run_dir/lineages/$lineage_name" \
            --debug >> {log:q} 2>&1

        # BUSCO failures can be non-fatal to rnaQUAST itself.
        if grep -Eiq 'ERROR![[:space:]]+busco[[:space:]]+failed|non-fatal ERROR:[[:space:]]*busco[[:space:]]+failed' {log:q}; then
            echo "ERROR: rnaQUAST reported a BUSCO failure." >> {log:q}
            exit 1
        fi

        python - "$run_dir/results" "$lineage_name" >> {log:q} 2>&1 <<'PY'
import re
import sys
from pathlib import Path

results = Path(sys.argv[1])
lineage = sys.argv[2]
report = results / "short_report.txt"

if not report.is_file() or report.stat().st_size == 0:
    raise SystemExit("ERROR: rnaQUAST did not produce short_report.txt.")

summaries = list(results.glob("tmp/*_BUSCO/short_summary.*.txt"))
if len(summaries) != 1:
    raise SystemExit("ERROR: Expected one retained BUSCO summary.")

text = summaries[0].read_text()
dataset = re.search(r"^# The lineage dataset is:\s+(\S+)", text, re.M)
if not dataset or dataset.group(1) != lineage:
    raise SystemExit(
        "ERROR: BUSCO summary lineage does not match the configured dataset."
    )

patterns = (
    r"Complete BUSCOs \(C\)",
    r"Complete and single-copy BUSCOs \(S\)",
    r"Complete and duplicated BUSCOs \(D\)",
    r"Fragmented BUSCOs \(F\)",
    r"Missing BUSCOs \(M\)",
    r"Total BUSCO groups searched",
)

counts = []
for pattern in patterns:
    match = re.search(r"^\s*(\d+)\s+" + pattern, text, re.M)
    if not match:
        raise SystemExit("ERROR: BUSCO summary is incomplete.")
    counts.append(int(match.group(1)))

complete, single, duplicated, fragmented, missing, total = counts
if (
    total < 1
    or complete != single + duplicated
    or complete + fragmented + missing != total
):
    raise SystemExit("ERROR: BUSCO summary counts are inconsistent.")

print("Validated rnaQUAST report and BUSCO summary:", summaries[0])
PY

        # Publish the directory only after required results are checked.
        mkdir -p {output.report_dir:q}
        cp -a "$run_dir/results/." {output.report_dir:q}
        """
