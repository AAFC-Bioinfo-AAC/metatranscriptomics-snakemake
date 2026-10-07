"""Join dbCAN recommendations to raw featureCounts paired-fragment counts."""
import argparse
import csv
import json
from pathlib import Path
import re


def read_table(path, required):
    with open(path, newline="") as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing required columns {sorted(required)}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}: malformed table row")
            yield row


def families(recommendation):
    result = set()
    for token in re.split(r"[|+;,]", recommendation):
        token = token.strip()
        match = re.fullmatch(r"((?:GH|GT|PL|CE|AA|CBM)\d+)(?:_\d+)*(?:\.hmm)?(?:\(\d+-\d+\))?", token)
        if not match:
            raise ValueError(f"Unrecognized dbCAN recommended family: {token!r}")
        # Collapse subfamilies and repeated domains to one parent-family membership.
        result.add(match.group(1))
    return sorted(result)


def read_counts(path, reference, selected):
    expected_genes = set(reference)
    counts, seen = {}, set()
    with open(path, newline="") as stream:
        reader = csv.reader((line for line in stream if not line.startswith("#")), delimiter="\t")
        header = next(reader, [])
        if header[:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"] or len(header) != 7:
            raise ValueError(f"{path}: expected a single-sample featureCounts table with seven columns")
        for row in reader:
            if len(row) != 7:
                raise ValueError(f"{path}: malformed count row")
            gene, value = row[0], row[6]
            if gene in seen:
                raise ValueError(f"{path}: duplicate count GeneID {gene}")
            if not re.fullmatch(r"\d+", value):
                raise ValueError(f"{path}: raw fragment count must be a nonnegative integer for {gene}")
            if gene in reference:
                annotation = reference[gene]
                coordinates = [annotation[key] for key in ("Chr", "Start", "End", "Strand")]
                length = int(annotation["End"]) - int(annotation["Start"]) + 1
                if row[1:5] != coordinates or row[5] != str(length):
                    raise ValueError(f"{path}: count coordinates/length differ from the protein/GFF reference for {gene}")
            seen.add(gene)
            if gene in selected:
                counts[gene] = int(value)
    if seen != expected_genes:
        raise ValueError(f"{path}: count gene IDs differ from the protein/GFF reference "
                         f"({len(expected_genes - seen)} missing, {len(seen - expected_genes)} extra)")
    return counts


def write_table(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def summarize(overview, mapping, samples, counts, min_tools, annotations, gene_counts, family_counts, summary):
    if min_tools not in (2, 3):
        raise ValueError("min_tools must be 2 or 3 (dbCAN recommendations require at least two methods)")
    if len(samples) != len(counts) or not samples or len(set(samples)) != len(samples):
        raise ValueError("Provide one count table per unique sample ID")
    if {"GeneID", "Family"}.intersection(samples):
        raise ValueError("Sample IDs GeneID and Family conflict with count-matrix headers")
    genes, proteins = {}, set()
    for row in read_table(mapping, {"GeneID", "ProteinID", "Chr", "Start", "End", "Strand"}):
        gene = row["GeneID"]
        if not gene or gene in genes or not row["ProteinID"] or row["ProteinID"] in proteins:
            raise ValueError(f"{mapping}: empty or duplicate gene/protein ID")
        genes[gene] = row
        proteins.add(row["ProteinID"])
    if not genes:
        raise ValueError(f"{mapping}: empty reference gene map")
    selected, seen = {}, set()
    for row in read_table(overview, {"Gene ID", "#ofTools", "Recommend Results"}):
        gene = row["Gene ID"]
        if gene in seen or gene not in genes:
            raise ValueError(f"{overview}: duplicate or unknown Gene ID {gene}")
        seen.add(gene)
        tools = row["#ofTools"]
        if tools not in ("0", "1", "2", "3"):
            raise ValueError(f"{overview}: invalid #ofTools {tools!r}")
        if int(tools) >= min_tools:
            selected[gene] = {"families": families(row["Recommend Results"]), "tools": tools,
                              "recommendation": row["Recommend Results"], "ec": row.get("EC#", "-")}
    gene_ids = sorted(selected)
    family_ids = sorted({family for gene in selected.values() for family in gene["families"]})
    matrix = {gene: [] for gene in gene_ids}
    totals = {family: [] for family in family_ids}
    for path in counts:
        sample_counts = read_counts(path, genes, selected)
        family_totals = dict.fromkeys(family_ids, 0)
        for gene in gene_ids:
            value = sample_counts[gene]
            matrix[gene].append(value)
            for family in selected[gene]["families"]:
                family_totals[family] += value
        for family in family_ids:
            totals[family].append(family_totals[family])
    write_table(annotations, ["GeneID", "ProteinID", "Chr", "Start", "End", "Strand", "Tools",
                              "Families", "RecommendedResults", "EC"],
                ([gene, genes[gene]["ProteinID"], genes[gene]["Chr"], genes[gene]["Start"],
                  genes[gene]["End"], genes[gene]["Strand"], selected[gene]["tools"],
                  ";".join(selected[gene]["families"]), selected[gene]["recommendation"],
                  selected[gene]["ec"]] for gene in gene_ids))
    write_table(gene_counts, ["GeneID", *samples], ([gene, *matrix[gene]] for gene in gene_ids))
    write_table(family_counts, ["Family", *samples], ([family, *totals[family]] for family in family_ids))
    Path(summary).parent.mkdir(parents=True, exist_ok=True)
    Path(summary).write_text(json.dumps({"samples": samples, "reference_genes": len(genes),
        "annotated_genes_before_filter": len(seen), "accepted_cazyme_genes": len(selected),
        "families": len(family_ids), "min_tools": min_tools,
        "count_unit": "raw assigned paired fragments",
        "family_assignment": "full count to each unique parent family; family totals overlap",
        "gene_counts": str(gene_counts), "family_counts": str(family_counts)}, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("overview", "mapping", "annotations", "gene-counts", "family-counts", "summary"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--counts", nargs="+", required=True)
    parser.add_argument("--min-tools", type=int, default=2)
    summarize(**vars(parser.parse_args()))
