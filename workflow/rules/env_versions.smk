'''
    Filename: env_versions.smk
    Author: Katherine James-Gzyl
    Date created: 2025/10/15
    Updated: 2026/10/08
    Snakemake version: 9.20.0
'''

import html
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


# Track edits to the workflow's environment definitions.
ENV_VERSION_DEFINITIONS = sorted(
    str(path)
    for pattern in ("*.yaml", "*.yml")
    for path in (Path(workflow.basedir) / "envs").glob(pattern)
)

# These are package names, which may differ from executable names.
# For example, featureCounts is supplied by the subread package.
KEY_BIOINFORMATICS_PACKAGES = (
    "bedtools", "bowtie2", "bracken", "busco", "dbcan", "diamond",
    "fastp", "hmmer", "kma", "kraken2", "megahit", "pigz",
    "prodigal", "pyhmmer", "rgi", "rnaquast", "samtools",
    "sortmerna", "spades", "subread",
)


rule software_report:
    input:
        env_definitions = ENV_VERSION_DEFINITIONS
    output:
        summary = f"{SOFTWARE_VERSIONS}/software_versions_summary.txt"
    params:
        conda_prefix = lambda wc: workflow.persistence.conda_env_path
    log:
        f"{LOG_DIR}/software_versions/software_report.log"
    threads: 1
    run:
        env_base = Path(params.conda_prefix).absolute()
        Path(output.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)

        with open(log[0], "w", encoding="utf-8") as logfile:
            if not env_base.is_dir():
                message = f"Conda environment cache not found: {env_base}"
                print(message, file=logfile)
                raise ValueError(
                    message + "; check the effective --conda-prefix."
                )

            environments = sorted(
                path for path in env_base.iterdir()
                if path.is_dir() and (path / "conda-meta").is_dir()
            )

            if not environments:
                message = (
                    f"No installed Conda environments found under {env_base}"
                )
                print(message, file=logfile)
                raise ValueError(
                    message + "; create the workflow environments first."
                )

            conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda")
            if not conda_exe:
                raise ValueError(
                    "Conda is unavailable; activate the Snakemake environment."
                )

            with open(output.summary, "w", encoding="utf-8") as summary:
                print("Conda environment package inventory", file=summary)
                print("===================================", file=summary)
                print(
                    f"Generated (UTC): "
                    f"{datetime.now(timezone.utc).isoformat(timespec='seconds')}",
                    file=summary,
                )
                print(f"Conda environment cache: {env_base}", file=summary)
                print(
                    "Scope: all installed environments directly beneath this cache.",
                    file=summary,
                )
                print(
                    "Includes cached environments from earlier runs; "
                    "installed packages do not prove tool execution.",
                    file=summary,
                )

                for environment in environments:
                    print(
                        f"Listing packages: {environment}",
                        file=logfile,
                        flush=True,
                    )

                    result = subprocess.run(
                        [
                            conda_exe,
                            "list",
                            "--prefix",
                            str(environment),
                            "--json",
                        ],
                        stdout=subprocess.PIPE,
                        stderr=logfile,
                        text=True,
                        encoding="utf-8",
                    )

                    if result.returncode != 0:
                        raise ValueError(
                            f"conda list failed for {environment} "
                            f"(exit {result.returncode}); see {log[0]}"
                        )

                    try:
                        packages = json.loads(result.stdout)
                    except json.JSONDecodeError as error:
                        raise ValueError(
                            f"Invalid conda list JSON for {environment}"
                        ) from error

                    if not isinstance(packages, list) or not packages:
                        raise ValueError(
                            f"No package records returned for {environment}"
                        )

                    rows = []
                    for package in packages:
                        if not isinstance(package, dict):
                            raise ValueError(
                                f"Invalid package record in {environment}"
                            )

                        name = package.get("name")
                        version = package.get("version")

                        if (
                            not isinstance(name, str)
                            or not name
                            or not isinstance(version, str)
                            or not version
                        ):
                            raise ValueError(
                                f"Package name/version missing in {environment}"
                            )

                        build = package.get(
                            "build_string",
                            package.get("build", "unknown"),
                        )
                        channel = package.get("channel", "unknown")
                        row = (name, version, str(build), str(channel))

                        if any(
                            any(char in field for char in "\t\r\n")
                            for field in row
                        ):
                            raise ValueError(
                                f"Invalid package field in {environment}"
                            )

                        rows.append(row)

                    print(
                        f"\n### Environment: {environment.name}",
                        file=summary,
                    )
                    print(f"Path: {environment}", file=summary)
                    print("Package\tVersion\tBuild\tChannel", file=summary)

                    for row in sorted(
                        rows,
                        key=lambda item: item[0].casefold(),
                    ):
                        print("\t".join(row), file=summary)


rule filter_key_bioinformatics_versions:
    input:
        summary = f"{SOFTWARE_VERSIONS}/software_versions_summary.txt"
    output:
        key_versions = f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.txt"
    params:
        key_packages = KEY_BIOINFORMATICS_PACKAGES
    threads: 1
    run:
        key_packages = {
            name.casefold() for name in params.key_packages
        }
        lines = Path(input.summary).read_text(
            encoding="utf-8"
        ).splitlines()

        selected = [
            "Key bioinformatics packages installed in cached Conda environments",
            "================================================================",
            f"Source: {input.summary}",
            "Installed packages do not prove tool execution.",
        ]

        section = []
        emitted = False
        found = False

        for line in lines:
            if line.startswith(
                ("Generated (UTC):", "Conda environment cache:")
            ):
                selected.append(line)

            elif line.startswith("### Environment:"):
                section = [line]
                emitted = False

            elif line.startswith("Path:"):
                section.append(line)

            elif (
                line.strip()
                and line.split()[0].casefold() in key_packages
            ):
                if not emitted:
                    selected.extend([
                        "",
                        *section,
                        "Package\tVersion\tBuild\tChannel",
                    ])
                    emitted = True

                selected.append(line)
                found = True

        if not found:
            selected.extend([
                "",
                "No matching bioinformatics packages were found.",
            ])

        Path(output.key_versions).parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        Path(output.key_versions).write_text(
            "\n".join(selected) + "\n",
            encoding="utf-8",
        )


rule filter_key_bioinformatics_html:
    input:
        key_versions = f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.txt"
    output:
        html = report(
            f"{SOFTWARE_VERSIONS}/key_bioinformatics_software.html",
            category="Software versions",
        )
    threads: 1
    run:
        escaped_text = html.escape(
            Path(input.key_versions).read_text(encoding="utf-8")
        )

        document = (
            '<!doctype html>\n'
            '<html lang="en">\n'
            '<head>\n'
            '<meta charset="utf-8">\n'
            '<title>Bioinformatics package versions</title>\n'
            '</head>\n'
            '<body>\n'
            '<h1>Bioinformatics package versions</h1>\n'
            f'<pre>{escaped_text}</pre>\n'
            '</body>\n'
            '</html>\n'
        )

        Path(output.html).parent.mkdir(parents=True, exist_ok=True)
        Path(output.html).write_text(document, encoding="utf-8")
