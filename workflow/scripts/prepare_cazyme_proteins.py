"""Use Prodigal's GFF IDs as FASTA IDs so dbCAN and featureCounts can join."""
import argparse
import csv
import re
from pathlib import Path


def fasta_records(path):
    header, sequence = None, []
    with open(path) as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            elif header is None:
                raise ValueError(f"{path}: sequence before FASTA header")
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def prepare(proteins, gff, output, mapping):
    genes = {}
    with open(gff) as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gff}: expected nine GFF columns")
            if fields[2] != "CDS":
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            gene = attrs.get("ID", "")
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", gene) or len(gene) > 50:
                raise ValueError(f"{gff}: missing or unsafe CDS ID {gene!r}")
            if gene in genes:
                raise ValueError(f"{gff}: duplicate CDS ID {gene}")
            genes[gene] = [fields[0], fields[3], fields[4], fields[6]]
    if not genes:
        raise ValueError(f"{gff}: no predicted CDSs")
    for path in (output, mapping):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    seen, protein_ids = set(), set()
    with open(output, "w") as fasta, open(mapping, "w", newline="") as table:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(["GeneID", "ProteinID", "Chr", "Start", "End", "Strand"])
        for header, sequence in fasta_records(proteins):
            protein = header.split()[0]
            match = re.search(r"(?:^|[;\s])ID=([^;\s]+)", header)
            if not match:
                raise ValueError(f"{proteins}: missing Prodigal ID= attribute for {protein}")
            gene = match.group(1)
            if gene not in genes:
                raise ValueError(f"{proteins}: gene {gene} absent from GFF")
            if gene in seen or protein in protein_ids:
                raise ValueError(f"{proteins}: duplicate protein or gene ID {protein}/{gene}")
            if not sequence or not sequence.rstrip("*"):
                raise ValueError(f"{proteins}: empty protein {protein}")
            seen.add(gene)
            protein_ids.add(protein)
            # A short, unique GFF ID also avoids dbCAN's protein-ID truncation.
            fasta.write(f">{gene}\n")
            for start in range(0, len(sequence), 80):
                fasta.write(sequence[start:start + 80] + "\n")
            writer.writerow([gene, protein, *genes[gene]])
    if seen != set(genes):
        raise ValueError(f"{proteins}: {len(set(genes) - seen)} GFF CDSs lack proteins")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("proteins", "gff", "output", "mapping"):
        parser.add_argument(f"--{name}", required=True)
    prepare(**vars(parser.parse_args()))
