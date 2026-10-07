"""Run all three dbCAN methods; reject logged errors and incomplete results."""
import argparse
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

DATABASE_FILES = ("CAZy.dmnd", "dbCAN.hmm", "dbCAN-sub.hmm", "fam-substrate-mapping.tsv")
RESULT_COLUMNS = {
    "diamond.out": {"Gene ID", "CAZy ID"},
    "dbCAN_hmm_results.tsv": {"Target Name", "HMM Name", "Target From", "Target To", "i-Evalue"},
    "dbCANsub_hmm_results.tsv": {"Target Name", "Subfam Name", "Subfam EC", "Target From", "Target To", "i-Evalue"},
    "overview.tsv": {"Gene ID", "#ofTools", "Recommend Results", "dbCAN_hmm", "dbCAN_sub", "DIAMOND"},
}


def validate_results(directory, log):
    with open(log) as stream:
        for line in stream:
            if re.search(r" - (ERROR|CRITICAL) - ", line):
                raise ValueError(f"dbCAN logged an error; inspect {log}: {line.strip()}")
    observed = set()
    overview = set()
    for filename, columns in RESULT_COLUMNS.items():
        with open(Path(directory) / filename, newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if not columns.issubset(reader.fieldnames or []):
                raise ValueError(f"{filename}: missing expected dbCAN 5 columns")
            id_column = "Gene ID" if "Gene ID" in columns else "Target Name"
            for row in reader:
                if None in row or any(row.get(key) is None for key in columns):
                    raise ValueError(f"{filename}: malformed annotation row")
                gene = row[id_column]
                if not gene:
                    raise ValueError(f"{filename}: empty gene ID")
                if filename == "overview.tsv":
                    if gene in overview:
                        raise ValueError(f"{filename}: duplicate gene {gene}")
                    overview.add(gene)
                else:
                    observed.add(gene)
    if overview != observed:
        raise ValueError("dbCAN overview does not cover exactly the genes in its method results")


def annotate(proteins, db_dir, output, log, threads, database_release):
    db_dir = Path(db_dir).resolve()
    output, log = Path(output), Path(log)
    output.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name in DATABASE_FILES:
        path = db_dir / name
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Required dbCAN database file missing or empty: {path}")
        with path.open("rb") as stream:
            digest = hashlib.sha256()
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        manifest[name] = {"sha256": digest.hexdigest(), "bytes": path.stat().st_size}
    with tempfile.TemporaryDirectory(prefix=".dbcan-", dir=output.parent) as temporary:
        results = Path(temporary) / "annotation"
        command = ["run_dbcan", "CAZyme_annotation", "--input_raw_data", str(Path(proteins).resolve()),
                   "--mode", "protein", "--output_dir", str(results), "--db_dir", str(db_dir),
                   "--methods", "diamond,hmm,dbCANsub", "--threads", str(threads), "--log-level", "INFO"]
        with log.open("w") as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
        validate_results(results, log)
        # Keep the original method outputs as well as the merged overview.
        shutil.copy2(log, results / "run_dbcan.log")
        provenance = {"dbcan_version": version("dbcan"), "database_directory": str(db_dir),
                      "database_release": database_release, "database_files": manifest,
                      "command": command}
        (results / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        if output.exists():
            raise ValueError(f"Refusing to replace an existing annotation directory: {output}")
        shutil.move(str(results), str(output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("proteins", "db-dir", "output", "log", "database-release"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--threads", type=int, required=True)
    annotate(**vars(parser.parse_args()))
