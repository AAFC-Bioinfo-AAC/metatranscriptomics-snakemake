"""Merge per-sample dbcan_utils RPM tables without rescaling their values."""
import argparse
import csv
import json
import math
from pathlib import Path


def merge(samples, directories, output, summary):
    if not samples or len(samples) != len(directories) or len(set(samples)) != len(samples):
        raise ValueError("Provide one abundance directory per unique sample")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    tables = (("fam_abund.out", "cazyme_family_RPM.tsv", "Family"),
              ("subfam_abund.out", "cazyme_subfamily_RPM.tsv", "Subfamily"),
              ("fam_substrate_abund.out", "cazyme_substrate_RPM.tsv", "Substrate"))
    sizes = {}
    for source_name, output_name, label in tables:
        matrices = []
        for directory in directories:
            path = Path(directory) / source_name
            values = {}
            with path.open(newline="") as stream:
                reader = csv.reader(stream, delimiter="\t")
                header = next(reader, None)
                if header is not None and (len(header) < 2 or "Abundance" not in header[1]):
                    raise ValueError(f"{path}: unexpected RPM table header")
                for row in reader:
                    if len(row) < 2 or not row[0] or row[0] in values:
                        raise ValueError(f"{path}: empty, duplicate or malformed feature")
                    value = float(row[1])
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(f"{path}: invalid RPM value")
                    values[row[0]] = value
            matrices.append(values)
        features = sorted(set().union(*(set(m) for m in matrices)))
        with (output / output_name).open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow([label, *samples])
            for feature in features:
                writer.writerow([feature, *(m.get(feature, 0.0) for m in matrices)])
        sizes[label] = len(features)
    Path(summary).write_text(json.dumps({"samples": samples, "features": sizes,
        "normalization": "dbcan_utils RPM; no further scaling",
        "denominator": "sum of counts assigned to all predicted genes per sample",
        "count_unit": "filtered mapped reads, not paired fragments",
        "strandedness": "unstranded",
        "family_definition": "dbcan_utils assignments; families can overlap"
    }, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--directories", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--summary", required=True)
    merge(**vars(p.parse_args()))
