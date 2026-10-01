from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import cpgi_nr_prod as prod
import pancgi_contract as contract
import pancgi_core as core
import pancgi_hal_runtime as hal_runtime
import pancgi_mapping as mapping


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def validate_label(label: str) -> None:
    if not label or any(char in label for char in "/\\\t\r\n"):
        raise ValueError(f"Invalid catalog label: {label!r}")


def prepare_bed4(row: dict[str, str], path: Path, contig_index: dict[tuple[str, str], str]) -> tuple[set[str], str]:
    records = []
    identifiers = set()
    for record in prod.iter_cpgi_bed(row["bed"], row["kind"]):
        fid = str(record["fid"])
        parsed = core.parse_fid(fid)
        if str(parsed["sample"]) != row["graph_sample"] or str(parsed["hap"]) != row["graph_hap"]:
            raise ValueError(f"{row['label']}: BED identifier does not match graph_sample/graph_hap: {fid}")
        if fid in identifiers:
            raise ValueError(f"{row['label']}: duplicate CpGI identifier: {fid}")
        identifiers.add(fid)
        contig_id = str(parsed["contig"])
        key = (row["graph_sample"], contig_id)
        if key not in contig_index:
            raise ValueError(f"{row['label']}: internal contig has no confirmed HAL mapping: {key}")
        hal_sequence = contig_index[key]
        records.append(f"{hal_sequence}\t{record['start0']}\t{record['end0']}\t{fid}\n")
    atomic_text(path, "".join(records))
    return identifiers, sha256_file(path)


def inspect_image(docker_bin: str, image: str) -> str:
    return hal_runtime.inspect_image(docker_bin, image)


def parse_psl(path: Path, identifiers: set[str]) -> tuple[int, Counter, list[str]]:
    counts = Counter()
    unknown = []
    records = 0
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 22:
                raise ValueError(f"{path}: line {line_number} has {len(fields)} columns; expected 22")
            fid = fields[0]
            if fid not in identifiers:
                unknown.append(fid)
            counts[fid] += 1
            records += 1
    if unknown:
        raise ValueError(f"{path}: PSL contains identifiers absent from source BED4: {unknown[:5]}")
    return records, counts, unknown


def signature(
    row: dict[str, str],
    bed_sha256: str,
    hal: Path,
    target: str,
    image: str,
    image_id: str,
    executable: str,
) -> dict[str, object]:
    native_identity = None
    if os.environ.get("PANCGI_HAL_RUNTIME", "docker") == "native":
        program = Path(hal_runtime.container_command("", image, [], executable, [])[0]).resolve()
        native_identity = {"path": str(program), "file_state": contract.file_state(program)}
    return {
        "label": row["label"],
        "source_genome": row["hal_genome"],
        "target_genome": target,
        "bed4_sha256": bed_sha256,
        "hal_path": str(hal.resolve()),
        "hal_file_state": contract.file_state(hal),
        "docker_image": image,
        "docker_image_id": image_id,
        "hal_liftover": executable,
        "hal_liftover_identity": native_identity,
        "parameters": ["--bedType", "4", "--outPSLWithName"],
    }


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def run_one(
    row: dict[str, str],
    hal: Path,
    target: str,
    root: Path,
    docker_bin: str,
    image: str,
    image_id: str,
    executable: str,
    overwrite: bool,
    contig_index: dict[tuple[str, str], str],
) -> dict[str, object]:
    label = row["label"]
    bed_path = root / "bed4" / f"{label}.bed4"
    psl_path = root / "psl" / f"{label}.to_primary.psl"
    metadata_path = root / "metadata" / f"{label}.json"
    log_path = root / "logs" / f"{label}.log"
    identifiers, bed_sha256 = prepare_bed4(row, bed_path, contig_index)
    expected = signature(row, bed_sha256, hal, target, image, image_id, executable)
    status = "generated"
    if psl_path.is_file() or metadata_path.is_file():
        if psl_path.is_file() and metadata_path.is_file() and load_json(metadata_path) == expected:
            records, counts, _ = parse_psl(psl_path, identifiers)
            status = "validated_existing"
        elif not overwrite:
            raise RuntimeError(f"{label}: existing HAL output does not match current inputs")
    if status == "generated":
        if not identifiers:
            atomic_text(psl_path, "")
            atomic_text(log_path, "Empty CGI input; no HAL projection required.\n")
            records, counts = 0, Counter()
        else:
            psl_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            temporary = psl_path.with_name(psl_path.name + f".partial.{os.getpid()}")
            if temporary.exists():
                temporary.unlink()
            mounts = [
                (hal.parent.resolve(), "ro"),
                (bed_path.parent.resolve(), "ro"),
                (psl_path.parent.resolve(), "rw"),
            ]
            command = hal_runtime.container_command(docker_bin, image, mounts, executable, [
                "--bedType",
                "4",
                "--outPSLWithName",
                str(hal.resolve()),
                row["hal_genome"],
                str(bed_path.resolve()),
                target,
                str(temporary.resolve()),
            ])
            result = subprocess.run(command, check=False, capture_output=True, text=True)
            atomic_text(log_path, result.stdout + result.stderr)
            if result.returncode != 0:
                if temporary.exists():
                    temporary.unlink()
                raise RuntimeError(f"{label}: halLiftover failed with exit code {result.returncode}; see {log_path}")
            if not temporary.exists():
                raise RuntimeError(f"{label}: halLiftover did not create its PSL output")
            records, counts, _ = parse_psl(temporary, identifiers)
            os.replace(temporary, psl_path)
        atomic_text(metadata_path, json.dumps(expected, indent=2, sort_keys=True) + "\n")
    mapped = set(counts)
    unmapped = sorted(identifiers - mapped)
    multimapped = sorted(fid for fid, count in counts.items() if count > 1)
    return {
        "label": label,
        "kind": row["kind"],
        "hal_genome": row["hal_genome"],
        "target_genome": target,
        "input_intervals": len(identifiers),
        "psl_records": records,
        "mapped_ids": len(mapped),
        "unmapped_ids": len(unmapped),
        "multimap_ids": len(multimapped),
        "status": status,
        "bed4": str(bed_path),
        "psl": str(psl_path),
        "unmapped": unmapped,
        "multimapped": [(fid, counts[fid]) for fid in multimapped],
    }


def write_reports(root: Path, rows: list[dict[str, object]]) -> None:
    summary_path = root / "hal_liftover_summary.tsv"
    fields = ["label", "kind", "hal_genome", "target_genome", "input_intervals", "psl_records", "mapped_ids", "unmapped_ids", "multimap_ids", "status", "bed4", "psl"]
    lines = ["\t".join(fields) + "\n"]
    for row in rows:
        lines.append("\t".join(str(row[field]) for field in fields) + "\n")
    atomic_text(summary_path, "".join(lines))
    unmapped_path = root / "hal_liftover_unmapped_ids.tsv.gz"
    unmapped_tmp = unmapped_path.with_name(unmapped_path.name + f".partial.{os.getpid()}")
    with gzip.open(unmapped_tmp, "wt", encoding="utf-8", newline="") as handle:
        handle.write("label\tfid\n")
        for row in rows:
            for fid in row["unmapped"]:
                handle.write(f"{row['label']}\t{fid}\n")
    os.replace(unmapped_tmp, unmapped_path)
    multimap_path = root / "hal_liftover_multimap_ids.tsv.gz"
    multimap_tmp = multimap_path.with_name(multimap_path.name + f".partial.{os.getpid()}")
    with gzip.open(multimap_tmp, "wt", encoding="utf-8", newline="") as handle:
        handle.write("label\tfid\tn_psl_records\n")
        for row in rows:
            for fid, count in row["multimapped"]:
                handle.write(f"{row['label']}\t{fid}\t{count}\n")
    os.replace(multimap_tmp, multimap_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--hal", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--docker-bin", default="docker")
    parser.add_argument("--docker-image", required=True)
    parser.add_argument("--hal-stats", default="/opt/hal/bin/halStats")
    parser.add_argument("--hal-liftover", default="/opt/hal/bin/halLiftover")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    hal = Path(args.hal).resolve()
    root = Path(args.out_dir).resolve()
    if not hal.is_file() or hal.stat().st_size == 0:
        raise FileNotFoundError(hal)
    if args.threads < 1:
        raise ValueError("threads must be positive")
    rows = list(prod.iter_internal_catalog(args.catalog))
    contig_rows = mapping.read_tsv(args.contigs, mapping.CONTIG_COLUMNS)
    contig_index = {}
    for contig_row in contig_rows:
        if contig_row["mapping_status"] != "confirmed":
            raise ValueError(f"HAL liftover requires mapping_status=confirmed: line {contig_row['_line_number']}")
        key = (contig_row["genome_id"], contig_row["contig_id"])
        if key in contig_index:
            raise ValueError(f"Duplicate HAL contig mapping: {key}")
        contig_index[key] = contig_row["hal_sequence"]
    primary_rows = [row for row in rows if row["kind"] == "reference_primary"]
    if len(primary_rows) != 1:
        raise ValueError("Catalog requires exactly one primary reference")
    target_genome = primary_rows[0]["hal_genome"]
    image_id = inspect_image(args.docker_bin, args.docker_image)
    genomes = set(hal_runtime.hal_genomes(args.hal, args.docker_bin, args.docker_image, args.hal_stats))
    expected = {row["hal_genome"] for row in rows}
    absent = sorted(expected - genomes)
    if absent:
        raise ValueError(f"Catalog HAL genomes absent from genome list: {absent}")
    if target_genome not in genomes:
        raise ValueError(f"Primary target genome absent from genome list: {target_genome}")
    tasks = []
    for row in rows:
        validate_label(row["label"])
        if row["kind"] == "reference_primary":
            if row["hal_genome"] != target_genome:
                raise ValueError(f"{row['label']}: primary reference row does not match target genome")
            continue
        tasks.append(row)
    root.mkdir(parents=True, exist_ok=True)
    completed = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(
                run_one,
                row,
                hal,
                target_genome,
                root,
                args.docker_bin,
                args.docker_image,
                image_id,
                args.hal_liftover,
                args.overwrite,
                contig_index,
            ): row["label"]
            for row in tasks
        }
        for future in as_completed(futures):
            completed.append(future.result())
    order = {row["label"]: index for index, row in enumerate(rows)}
    completed.sort(key=lambda row: order[row["label"]])
    write_reports(root, completed)
    print(json.dumps({"status": "pass", "tasks": len(completed), "out_dir": str(root)}, sort_keys=True))


if __name__ == "__main__":
    main()
