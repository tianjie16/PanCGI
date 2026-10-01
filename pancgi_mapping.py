from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sqlite3
import time
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

import pancgi_hal_runtime as hal_runtime


GENOME_COLUMNS = [
    "genome_id",
    "role",
    "individual_id",
    "haplotype_id",
    "hal_genome",
    "cpgi_bed",
    "sv_file",
    "sv_format",
    "sv_schema",
    "sv_evidence_status",
]

CONTIG_COLUMNS = [
    "genome_id",
    "contig_id",
    "gfa_path_id",
    "hal_sequence",
    "cpgi_contig",
    "sv_contig",
    "mapping_status",
    "mapping_method",
    "mapping_note",
]

GFA_INVENTORY_COLUMNS = [
    "gfa_path_id",
    "record_type",
    "raw_path_name",
    "w_sample_id",
    "w_haplotype_index",
    "sequence_id",
    "start0",
    "end0",
    "path_length_bp",
    "fragment_count",
    "source_lines",
    "overlap_status",
]

HAL_INVENTORY_COLUMNS = [
    "hal_genome",
    "hal_sequence",
    "length_bp",
    "top_segments",
    "bottom_segments",
]

HAL_GENOME_INVENTORY_COLUMNS = [
    "hal_genome",
    "n_sequences",
    "total_sequence_bp",
    "sequence_status",
]

VALID_ROLES = {"primary_reference", "comparison_reference", "sample"}
VALID_SV_EVIDENCE = {"available", "confirmed_none", "unavailable"}
INTERNAL_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
WALK_TOKEN_RE = re.compile(r"([<>])([^<>,\s]+)")
P_TOKEN_RE = re.compile(r"([^,;\s]+)([+-])")
ACCEPTED_PATH_STATUSES = {"not_applicable", "zero", "zero_verified_from_links", "zero_single_segment", "requires_link_verification"}


def open_text(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def stable_path_id(record_type: str, values: Sequence[str]) -> str:
    payload = "\0".join([record_type, *map(str, values)]).encode("utf-8")
    return "GP_" + hashlib.sha256(payload).hexdigest()[:24]


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    with temporary.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    os.replace(temporary, path)


def read_tsv(path: str, required: Sequence[str]) -> List[Dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if fields != list(required):
            raise ValueError(f"{path}: expected exact columns {list(required)}, got {fields}")
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'{path}:{line_number}: incorrect row width')
            if any(value != value.strip() for value in row.values()):
                raise ValueError(f'{path}:{line_number}: unexpected surrounding whitespace')
            normalized = dict(row)
            normalized["_line_number"] = str(line_number)
            rows.append(normalized)
    if not rows:
        raise ValueError(f"{path} contains no data rows")
    return rows


def segment_length(fields: Sequence[str], gfa_major: int) -> int:
    if gfa_major == 2:
        if len(fields) < 4:
            raise ValueError("Malformed GFA2 S record")
        declared = int(fields[2])
        sequence = fields[3]
        if sequence != "*" and len(sequence) != declared:
            raise ValueError(f"GFA2 segment {fields[1]} length does not match sequence")
        return declared
    if len(fields) < 3:
        raise ValueError("Malformed GFA1 S record")
    if fields[2] != "*":
        return len(fields[2])
    for tag in fields[3:]:
        if tag.startswith("LN:i:"):
            return int(tag[5:])
    raise ValueError(f"GFA1 segment {fields[1]} has no sequence or LN:i tag")


def canonical_numeric(value: str) -> bool:
    return value.isdigit() and str(int(value)) == value


def source_fingerprint(path: str) -> Dict[str, object]:
    source = Path(path).resolve()
    stat = source.stat()
    return dict(path=str(source), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                ctime_ns=stat.st_ctime_ns, device=stat.st_dev, inode=stat.st_ino)


def check_source_fingerprint(path: str, expected: Dict[str, object], *, cross_node: bool = False) -> Dict[str, object]:
    try:
        actual = source_fingerprint(path)
    except OSError as exc:
        raise ValueError("GFA source changed; rerun inventory-gfa") from exc
    valid = isinstance(expected, dict) and set(expected) == set(actual)
    valid = valid and all(type(expected[key]) is type(value) for key, value in actual.items())
    compared = set(actual) - ({"device"} if cross_node else set())
    if not valid or any(actual[key] != expected[key] for key in compared):
        raise ValueError("GFA source changed or inventory fingerprint mismatched; rerun inventory-gfa")
    return actual


def inventory_metadata_path(inventory: str) -> Path:
    return Path(str(Path(inventory).resolve()) + ".json")


def write_json(path: Path, value: Dict[str, object]) -> None:
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


class ProgressLog:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = path.open("at", encoding="utf-8", buffering=1)
        self.last_update = time.monotonic()

    def emit(self, phase: str, event: str, **values) -> None:
        self.handle.write(json.dumps(dict(time=time.time(), phase=phase, event=event, **values), sort_keys=True) + "\n")
        self.last_update = time.monotonic()

    def update(self, phase: str, line_number: int, **values) -> None:
        if line_number % 100000 == 0 or time.monotonic() - self.last_update >= 10:
            self.emit(phase, "progress", source_lines=line_number, **values)

    def close(self) -> None:
        self.handle.close()


def p_segment_tokens(path_field: str) -> Iterator[Tuple[str, str]]:
    found = False
    cursor = 0
    for match in P_TOKEN_RE.finditer(path_field):
        separator = path_field[cursor:match.start()]
        if ";" in separator:
            raise ValueError("GFA P jump separators are unsupported for coordinate unfolding")
        if separator.strip(" ,\t\r\n"):
            raise ValueError("Malformed GFA P path")
        yield match.group(1), match.group(2)
        found = True
        cursor = match.end()
    if ";" in path_field[cursor:]:
        raise ValueError("GFA P jump separators are unsupported for coordinate unfolding")
    if path_field[cursor:].strip(" ,\t\r\n") or not found:
        raise ValueError("Malformed or empty GFA P path")


def canonical_link_key(from_name: str, from_orientation: str, to_name: str, to_orientation: str) -> Tuple[str, str, str, str]:
    if from_orientation not in {"+", "-"} or to_orientation not in {"+", "-"}:
        raise ValueError("Invalid GFA link orientation")
    reverse = (
        to_name,
        "+" if to_orientation == "-" else "-",
        from_name,
        "+" if from_orientation == "-" else "-",
    )
    direct = (from_name, from_orientation, to_name, to_orientation)
    return min(direct, reverse)


class LinkStore:
    def __init__(self, path: Path):
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA journal_mode=OFF")
        self.connection.execute("PRAGMA synchronous=OFF")
        self.connection.execute(
            "CREATE TABLE link (from_name TEXT, from_orientation TEXT, to_name TEXT, to_orientation TEXT, overlap TEXT, match_count INTEGER, PRIMARY KEY (from_name, from_orientation, to_name, to_orientation))"
        )
        self.pending: List[Tuple[str, str, str, str, str]] = []

    def add(self, key: Tuple[str, str, str, str], overlap: str) -> None:
        self.pending.append((*key, overlap))
        if len(self.pending) >= 100000:
            self.flush()

    def flush(self) -> None:
        if not self.pending:
            return
        self.connection.executemany(
            "INSERT INTO link VALUES (?, ?, ?, ?, ?, 1) ON CONFLICT(from_name, from_orientation, to_name, to_orientation) DO UPDATE SET overlap = CASE WHEN overlap = excluded.overlap THEN overlap ELSE '__CONFLICT__' END, match_count = match_count + 1",
            self.pending,
        )
        self.connection.commit()
        self.pending = []

    @lru_cache(maxsize=1000000)
    def status(self, key: Tuple[str, str, str, str]) -> str:
        self.flush()
        row = self.connection.execute(
            "SELECT overlap, match_count FROM link WHERE from_name = ? AND from_orientation = ? AND to_name = ? AND to_orientation = ?",
            key,
        ).fetchone()
        if row is None:
            return "missing_link"
        if int(row[1]) != 1:
            return "multiple_matching_links"
        if row[0] != "0M":
            return "link_overlap_unsupported"
        return "zero"

    def close(self) -> None:
        self.flush()
        self.status.cache_clear()
        self.connection.close()


def overlap_status(value: str, segment_count: int) -> str:
    if value == "*":
        return "zero_single_segment" if segment_count == 1 else "requires_link_verification"
    entries = value.split(",") if value else []
    if len(entries) != max(0, segment_count - 1):
        return "invalid_count"
    if all(entry == "0M" for entry in entries):
        return "zero"
    return "nonzero_unsupported"


def inventory_gfa(gfa: str, output: str, work_dir: str = "") -> Dict[str, object]:
    gfa_path = Path(gfa).resolve()
    if not gfa_path.is_file() or gfa_path.stat().st_size == 0:
        raise FileNotFoundError(gfa_path)
    fingerprint = source_fingerprint(str(gfa_path))
    output_path = Path(output).resolve()
    root = Path(work_dir).resolve() if work_dir else output_path.parent
    root.mkdir(parents=True, exist_ok=True)
    progress = ProgressLog(root / "progress.jsonl")
    gfa_major = 1
    stats = dict(n_segments=0, max_numeric=-1, all_numeric=True)
    grouped: Dict[str, Dict[str, object]] = {}
    line_number = 0
    progress.emit("inventory", "start", gfa=str(gfa_path))
    try:
        with open_text(str(gfa_path), "rt") as handle:
            for line_number, raw in enumerate(handle, start=1):
                progress.update("inventory", line_number, n_segments=stats["n_segments"], n_paths=len(grouped))
                if raw.startswith("S\t"):
                    fields = raw.split("\t", 2)
                    name = fields[1].rstrip("\r\n")
                    stats["n_segments"] += 1
                    if canonical_numeric(name):
                        stats["max_numeric"] = max(stats["max_numeric"], int(name))
                    else:
                        stats["all_numeric"] = False
                    continue
                if not raw.startswith(("H\t", "W\t", "P\t", "O\t")):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if fields[0] == "H":
                    for tag in fields[1:]:
                        if tag.startswith("VN:Z:"):
                            major = int(tag[5:].split(".", 1)[0])
                            if major not in {1, 2} or ((stats["n_segments"] or grouped) and major != gfa_major):
                                raise ValueError(f"Unsupported or inconsistent GFA major version at line {line_number}")
                            gfa_major = major
                elif fields[0] == "W":
                    if len(fields) < 7:
                        raise ValueError(f"Malformed GFA W record at line {line_number}")
                    _, sample, hap, sequence_id, start0, end0, walk = fields[:7]
                    path_id = stable_path_id("W", [sample, hap, sequence_id])
                    known = start0 != "*" and end0 != "*"
                    start, end = (int(start0), int(end0)) if known else ("", "")
                    if known and (start < 0 or end <= start):
                        raise ValueError(f"Invalid GFA W coordinates at line {line_number}")
                    row = grouped.setdefault(path_id, {
                        "gfa_path_id": path_id, "record_type": "W", "raw_path_name": "",
                        "w_sample_id": sample, "w_haplotype_index": hap, "sequence_id": sequence_id,
                        "start0": start, "end0": end, "path_length_bp": 0,
                        "fragment_count": 0, "source_lines": [],
                        "overlap_status": "not_applicable",
                    })
                    if not known or row["overlap_status"] == "w_coordinates_unknown_unsupported":
                        row.update(overlap_status="w_coordinates_unknown_unsupported", start0="", end0="", path_length_bp="")
                    else:
                        if row["fragment_count"] and row["end0"] != start:
                            raise ValueError(f"Non-contiguous W fragments for {path_id} at line {line_number}")
                        row["end0"] = end
                        row["path_length_bp"] += end - start
                    row["fragment_count"] += 1
                    row["source_lines"].append(line_number)
                elif (fields[0] == "P" and gfa_major == 1) or (fields[0] == "O" and gfa_major == 2):
                    record_type = fields[0]
                    if len(fields) < (4 if record_type == "P" else 3):
                        raise ValueError(f"Malformed GFA {record_type} record at line {line_number}")
                    path_name = fields[1]
                    path_id = stable_path_id(record_type, [path_name])
                    if path_id in grouped:
                        raise ValueError(f"Duplicate GFA {record_type} path {path_name!r}")
                    status = "gfa2_ordered_group_unsupported"
                    if record_type == "P":
                        if ";" in fields[2]:
                            raise ValueError("GFA P jump separators are unsupported for coordinate unfolding")
                        if fields[3] == "*":
                            status = "zero_single_segment" if P_TOKEN_RE.fullmatch(fields[2]) else "requires_link_verification"
                        else:
                            status = "zero" if all(v == "0M" for v in fields[3].split(",")) else "nonzero_unsupported"
                    grouped[path_id] = {
                        "gfa_path_id": path_id, "record_type": record_type, "raw_path_name": path_name,
                        "w_sample_id": "", "w_haplotype_index": "", "sequence_id": path_name,
                        "start0": 0, "end0": "", "path_length_bp": "",
                        "fragment_count": 1, "source_lines": [line_number],
                        "overlap_status": status,
                    }
        check_source_fingerprint(str(gfa_path), fingerprint)
        if not grouped:
            raise ValueError("GFA contains no explicit full-sequence W, P, or O paths; strict PanCGI genotyping cannot use an S/L-only or reference-only rGFA")
        stats["dense_numeric"] = stats["all_numeric"] and stats["max_numeric"] <= max(1024, stats["n_segments"] * 2)
        rows = []
        for path_id in sorted(grouped):
            row = grouped[path_id]
            row["source_lines"] = ",".join(map(str, row["source_lines"]))
            rows.append(row)
        write_tsv(output_path, rows, GFA_INVENTORY_COLUMNS)
        summary = {
            "schema_version": 2, "gfa": str(gfa_path), "gfa_major_version": gfa_major,
            "source_fingerprint": fingerprint, "node_stats": stats, "n_segments": stats["n_segments"],
            "n_paths": len(rows), "n_w_paths": sum(row["record_type"] == "W" for row in rows),
            "n_p_paths": sum(row["record_type"] == "P" for row in rows),
            "n_o_paths": sum(row["record_type"] == "O" for row in rows), "output": str(output_path),
            "metadata_file": str(inventory_metadata_path(output)), "source_read_passes": 1,
            "segment_index_builds": 0, "link_index_builds": 0, "progress_file": str(root / "progress.jsonl"),
        }
        write_json(inventory_metadata_path(output), summary)
        progress.emit("inventory", "complete", source_lines=line_number, n_segments=stats["n_segments"], n_paths=len(rows), source_read_passes=1)
        return summary
    except Exception as exc:
        progress.emit("inventory", "error", message=str(exc), source_lines=line_number)
        raise
    finally:
        progress.close()


def inventory_hal(args: argparse.Namespace) -> Dict[str, object]:
    genomes = hal_runtime.hal_genomes(args.hal, args.docker_bin, args.docker_image, args.hal_stats)
    rows: List[Dict[str, object]] = []
    genome_rows: List[Dict[str, object]] = []
    for genome in genomes:
        sequence_rows = hal_runtime.hal_sequence_stats(
            args.hal,
            genome,
            args.docker_bin,
            args.docker_image,
            args.hal_stats,
            allow_empty=True,
        )
        rows.extend(sequence_rows)
        genome_rows.append({
            "hal_genome": genome,
            "n_sequences": len(sequence_rows),
            "total_sequence_bp": sum(int(row["length_bp"]) for row in sequence_rows),
            "sequence_status": "available" if sequence_rows else "no_sequences",
        })
    write_tsv(Path(args.output).resolve(), rows, HAL_INVENTORY_COLUMNS)
    genome_output = Path(args.genome_output).resolve() if args.genome_output else Path(args.output).resolve().with_name("hal_genomes.tsv")
    write_tsv(genome_output, genome_rows, HAL_GENOME_INVENTORY_COLUMNS)
    return {"hal": str(Path(args.hal).resolve()), "n_genomes": len(genomes), "n_genomes_without_sequences": sum(not int(row["n_sequences"]) for row in genome_rows), "n_sequences": len(rows), "output": str(Path(args.output).resolve()), "genome_output": str(genome_output)}


def make_public_templates(args):
    from pancgi_contract import GENOMES, PATHS
    hal_rows = read_tsv(args.hal_inventory, HAL_INVENTORY_COLUMNS)
    gfa_rows = read_tsv(args.gfa_inventory, GFA_INVENTORY_COLUMNS)
    genomes = sorted({r['hal_genome'] for r in hal_rows})
    output = Path(args.out_dir).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    write_tsv(output / 'genomes.tsv', [dict(hal_genome=g, role='', cpgi_bed='', sv_tsv='') for g in genomes], GENOMES)
    write_tsv(output / 'paths.tsv', [dict(hal_genome=r['hal_genome'], hal_sequence=r['hal_sequence']) for r in hal_rows], PATHS)
    write_tsv(output / 'gfa_paths.tsv', [{k: r[k] for k in GFA_INVENTORY_COLUMNS} for r in gfa_rows], GFA_INVENTORY_COLUMNS)
    write_tsv(output / 'hal_sequences.tsv', hal_rows, HAL_INVENTORY_COLUMNS)
    sequence_counts = Counter(row['hal_sequence'] for row in hal_rows)
    w_paths = defaultdict(list)
    p_paths = defaultdict(list)
    for row in gfa_rows:
        if row['record_type'] == 'W':
            w_paths[(row['w_sample_id'], row['sequence_id'])].append(row)
        elif row['record_type'] == 'P':
            p_paths[row['raw_path_name']].append(row)
    columns = ['hal_genome', 'hal_sequence', 'hal_length_bp', 'candidate_status', 'candidate_count',
               'candidate_gfa_path_ids', 'candidate_gfa_path_id', 'candidate_record_type',
               'candidate_raw_path_name', 'candidate_w_sample_id', 'candidate_w_haplotype_index',
               'candidate_sequence_id', 'candidate_path_length_bp', 'candidate_overlap_status',
               'sequence_identity_checked']
    review = []
    for row in hal_rows:
        candidates = list(w_paths[(row['hal_genome'], row['hal_sequence'])])
        if sequence_counts[row['hal_sequence']] == 1:
            candidates.extend(p_paths[row['hal_sequence']])
        item = dict(hal_genome=row['hal_genome'], hal_sequence=row['hal_sequence'],
                    hal_length_bp=row['length_bp'], candidate_count=len(candidates),
                    candidate_gfa_path_ids=json.dumps(sorted(candidate['gfa_path_id'] for candidate in candidates)),
                    candidate_status='no_exact_literal_match', sequence_identity_checked='false')
        if len(candidates) == 1:
            item['candidate_status'] = 'exact_literal_candidate'
            for column in ('gfa_path_id', 'record_type', 'raw_path_name', 'w_sample_id',
                           'w_haplotype_index', 'sequence_id', 'path_length_bp', 'overlap_status'):
                item['candidate_' + column] = candidates[0][column]
        elif candidates:
            item['candidate_status'] = 'multiple_exact_literal_matches'
        review.append(item)
    write_tsv(output / 'mapping_review.tsv', review, columns)
    return dict(genome_n=len(genomes), sequence_n=len(hal_rows), gfa_path_n=len(gfa_rows),
                matching='explicit_user_assignment', mapping_review=str(output / 'mapping_review.tsv'),
                exact_literal_candidate_n=sum(row['candidate_status'] == 'exact_literal_candidate' for row in review))


def main() -> None:
    parser = argparse.ArgumentParser(description="PanCGI pangenome inventory and confirmed mapping interface")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("inventory-gfa")
    p.add_argument("--gfa", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--work-dir", default="")
    p.set_defaults(func=lambda a: inventory_gfa(a.gfa, a.output, a.work_dir))
    p = sub.add_parser("inventory-hal")
    p.add_argument("--hal", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--genome-output", default="")
    p.add_argument("--docker-bin", default="docker")
    p.add_argument("--docker-image", default="")
    p.add_argument("--hal-runtime", choices=['native','docker'], default=os.environ.get('PANCGI_HAL_RUNTIME', 'docker'))
    p.add_argument("--hal-stats", default="halStats")
    p.set_defaults(func=inventory_hal)
    p = sub.add_parser("make-mapping-template")
    p.add_argument("--gfa-inventory", required=True)
    p.add_argument("--hal-inventory", required=True)
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=make_public_templates)
    args = parser.parse_args()
    if hasattr(args, 'hal_runtime'):
        os.environ['PANCGI_HAL_RUNTIME'] = args.hal_runtime
    print(json.dumps(args.func(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
