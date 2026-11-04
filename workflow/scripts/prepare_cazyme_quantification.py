"""Apply RNA CAZyme tool support before dbcan_utils abundance estimation.

The gene-depth table still covers ALL predicted genes, preserving the RPM
denominator. CGC tables retain run_dbCAN's complete gene neighbourhoods.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import tempfile


def filtered_table(source, destination, gene_column, allowed, min_tools=None):
    kept = 0
    with open(source, newline="") as stream, open(destination, "w", newline="") as target:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or gene_column not in reader.fieldnames:
            raise ValueError(f"{source}: missing {gene_column}")
        if min_tools is not None and "#ofTools" not in reader.fieldnames:
            raise ValueError(f"{source}: missing #ofTools")
        writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        seen = set()
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"{source}: malformed row")
            gene = row[gene_column]
            if gene not in allowed:
                # Unknown overview genes are reference errors. Subfamily hits
                # from rejected genes are expected and simply filtered out.
                if min_tools is not None:
                    raise ValueError(f"{source}: unknown reference gene {gene}")
                continue
            if min_tools is not None:
                if gene in seen or row["#ofTools"] not in ("0", "1", "2", "3"):
                    raise ValueError(f"{source}: duplicate gene or invalid #ofTools")
                seen.add(gene)
                if int(row["#ofTools"]) < min_tools:
                    continue
            writer.writerow(row)
            kept += 1
    return kept


def prepare(annotation, mapping, output, min_tools):
    if min_tools not in (2, 3):
        raise ValueError("min_tools must be 2 or 3")
    annotation, output = Path(annotation), Path(output)
    with open(mapping, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if "GeneID" not in (reader.fieldnames or []):
            raise ValueError(f"{mapping}: missing GeneID")
        ids = [row["GeneID"] for row in reader]
    if not ids or any(not i for i in ids) or len(set(ids)) != len(ids):
        raise ValueError(f"{mapping}: empty or duplicate reference IDs")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".quantification-", dir=output.parent) as temporary:
        work = Path(temporary)
        accepted = filtered_table(annotation / "overview.tsv", work / "overview.tsv",
                                  "Gene ID", set(ids), min_tools)
        with (work / "overview.tsv").open() as stream:
            accepted_ids = {row["Gene ID"] for row in csv.DictReader(stream, delimiter="\t")}
        filtered_table(annotation / "dbCANsub_hmm_results.tsv", work / "dbCANsub_hmm_results.tsv",
                       "Target Name", accepted_ids)
        for name in ("cgc_standard_out.tsv", "substrate_prediction.tsv"):
            shutil.copy2(annotation / name, work / name)
        (work / "filter_summary.json").write_text(json.dumps({
            "reference_genes": len(ids), "accepted_cazyme_genes": accepted,
            "minimum_tools": min_tools,
            "rpm_denominator": "sum of counts assigned to all predicted genes in each sample",
            "cgc_counts": "complete run_dbcan gene neighbourhoods, without this family-level filter"
        }, indent=2) + "\n")
        if output.exists():
            raise ValueError(f"Refusing to replace an existing directory: {output}")
        shutil.move(str(work), str(output))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("annotation", "mapping", "output"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--min-tools", type=int, choices=(2, 3), default=2)
    prepare(**vars(p.parse_args()))
