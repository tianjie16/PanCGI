from __future__ import annotations

import argparse
import csv
import gzip
import mmap
import hashlib
import io
import json
import os
import re
import sqlite3
from multiprocessing.util import Finalize
from array import array
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Tuple

import pancgi_mapping as mapping


WALK_TOKEN_RE = re.compile(r"([<>])([^<>,\s]+)")
P_TOKEN_RE = re.compile(r"([^,;\s]+)([+-])")


def open_text(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def safe_component(value: str, field_name: str) -> str:
    value = str(value)
    if not value or value in {".", ".."} or "/" in value or "\\" in value or "\0" in value:
        raise ValueError(f"Unsafe or empty {field_name}: {value!r}")
    return value


def canonical_numeric(value: str) -> bool:
    return mapping.canonical_numeric(value)


def load_inventory_metadata(gfa: str, inventory: str) -> Dict[str, object]:
    path = mapping.inventory_metadata_path(inventory)
    try:
        metadata = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Missing or invalid GFA inventory metadata {path}; rerun inventory-gfa") from exc
    if not isinstance(metadata, dict) or metadata.get("schema_version") != 2 or metadata.get("gfa_major_version") not in {1, 2}:
        raise ValueError("Invalid GFA inventory metadata; rerun inventory-gfa")
    metadata["run_source_fingerprint"] = mapping.check_source_fingerprint(
        gfa, metadata.get("source_fingerprint"), cross_node=True)
    stats = metadata.get("node_stats")
    if not isinstance(stats, dict) or any(type(stats.get(key)) is not int for key in ("n_segments", "max_numeric")) or any(type(stats.get(key)) is not bool for key in ("all_numeric", "dense_numeric")):
        raise ValueError("Invalid GFA inventory node statistics; rerun inventory-gfa")
    if stats["n_segments"] <= 0:
        raise ValueError("GFA contains no S records")
    dense = stats["all_numeric"] and stats["max_numeric"] <= max(1024, stats["n_segments"] * 2)
    if stats["max_numeric"] < -1 or (stats["all_numeric"] and stats["max_numeric"] < 0) or stats["dense_numeric"] != dense:
        raise ValueError("Invalid GFA inventory node statistics; rerun inventory-gfa")
    return metadata


class SegmentIndex:
    def __init__(self, work_dir: Path, stats: Dict[str, object]):
        self.mode = "dense_numeric" if stats["dense_numeric"] else "sqlite_opaque"
        self.lengths = array("I")
        self.connection = None
        self.cursor = None
        self.n_segments = 0
        self.max_numeric = -1
        self.all_numeric = True
        self.pending = []
        self.map_handle = None
        self.map_path = work_dir / "segment_id_map.tsv.gz"
        work_dir.mkdir(parents=True, exist_ok=True)
        if self.mode == "dense_numeric":
            self.lengths = array("I", [0]) * (int(stats["max_numeric"]) + 1)
        else:
            database = work_dir / "segment_index.sqlite"
            if database.exists():
                database.unlink()
            self.connection = sqlite3.connect(database)
            self.connection.execute("PRAGMA journal_mode=OFF")
            self.connection.execute("PRAGMA synchronous=OFF")
            self.connection.execute("CREATE TABLE segment (name TEXT PRIMARY KEY, internal_id INTEGER UNIQUE, length_bp INTEGER NOT NULL)")
            self.map_handle = gzip.open(self.map_path, "wt", encoding="utf-8", newline="")
            self.writer = csv.writer(self.map_handle, delimiter="\t", lineterminator="\n")
            self.writer.writerow(["internal_node_id", "gfa_segment_id", "length_bp"])

    def add(self, fields: List[str], gfa_major: int, line_number: int) -> None:
        try:
            length = mapping.segment_length(fields, gfa_major)
            if fields[3 if gfa_major == 2 else 2] == '*':
                raise ValueError('Sequence-bearing GFA segments are required')
            name = fields[1]
            if not name or length <= 0:
                raise ValueError(f"Invalid GFA segment: name={name!r}, length={length}")
            numeric = canonical_numeric(name)
            if numeric:
                self.max_numeric = max(self.max_numeric, int(name))
            else:
                self.all_numeric = False
            if self.mode == "dense_numeric":
                if not numeric or int(name) >= len(self.lengths):
                    raise ValueError("GFA inventory node statistics mismatched; rerun inventory-gfa")
                node_id = int(name)
                if self.lengths[node_id] != 0:
                    raise ValueError(f"Duplicate GFA segment ID: {name}")
                self.lengths[node_id] = length
            else:
                internal_id = self.n_segments
                self.pending.append((name, internal_id, length))
                self.writer.writerow([internal_id, name, length])
                if len(self.pending) >= 100000:
                    self.flush()
            self.n_segments += 1
        except (ValueError, OverflowError) as exc:
            raise ValueError(f"GFA segment error at line {line_number}: {exc}") from exc

    def flush(self) -> None:
        if not self.pending:
            return
        try:
            self.connection.executemany("INSERT INTO segment VALUES (?, ?, ?)", self.pending)
            self.connection.commit()
        except sqlite3.IntegrityError as exc:
            raise ValueError("Duplicate GFA segment ID") from exc
        self.pending.clear()

    def finish(self, stats: Dict[str, object]) -> None:
        self.flush()
        for key in ("n_segments", "max_numeric", "all_numeric"):
            if getattr(self, key) != stats[key]:
                raise ValueError("GFA inventory node statistics mismatched; rerun inventory-gfa")
        if self.map_handle is not None:
            self.map_handle.close()
            self.map_handle = None
        if self.connection is not None:
            self.connection.execute("BEGIN")
            self.cursor = self.connection.cursor()

    @lru_cache(maxsize=1000000)
    def resolve(self, name: str) -> Tuple[int, int]:
        if self.mode == "dense_numeric":
            if not canonical_numeric(name):
                raise ValueError(f"Graph path references non-numeric segment in dense-numeric graph: {name!r}")
            node_id = int(name)
            if node_id >= len(self.lengths) or self.lengths[node_id] == 0:
                raise ValueError(f"Graph path references undefined segment {name!r}")
            return node_id, int(self.lengths[node_id])
        row = self.cursor.execute("SELECT internal_id, length_bp FROM segment WHERE name = ?", (name,)).fetchone()
        if row is None:
            raise ValueError(f"Graph path references undefined segment {name!r}")
        return int(row[0]), int(row[1])

    def close(self) -> None:
        self.resolve.cache_clear()
        if self.map_handle is not None:
            self.map_handle.close()
        if self.connection is not None:
            self.connection.close()


class ReadOnlySegmentIndex(SegmentIndex):
    def __init__(self, work_dir, stats):
        self.mode = 'dense_numeric' if stats['dense_numeric'] else 'sqlite_opaque'
        self.connection = None
        self.cursor = None
        self.map_handle = None
        self.lengths = ()
        self.length_file = None
        self.length_map = None
        if self.mode == 'dense_numeric':
            self.length_file = (work_dir / 'segment_lengths.bin').open('rb')
            self.length_map = mmap.mmap(self.length_file.fileno(), 0, access=mmap.ACCESS_READ)
            self.lengths = memoryview(self.length_map).cast('I')
            if len(self.lengths) != stats['max_numeric'] + 1:
                self.close()
                raise ValueError('Shared segment length index size mismatch')
        else:
            self.connection = sqlite3.connect((work_dir / 'segment_index.sqlite').resolve().as_uri() + '?mode=ro', uri=True)
            self.connection.execute('PRAGMA query_only=ON')
            self.connection.execute('BEGIN')
            self.cursor = self.connection.cursor()

    def close(self):
        super().close()
        if isinstance(self.lengths, memoryview):
            self.lengths.release()
            self.lengths = ()
        if self.length_map is not None:
            self.length_map.close()
            self.length_map = None
        if self.length_file is not None:
            self.length_file.close()
            self.length_file = None


class ReadOnlyLinkStore(mapping.LinkStore):
    def __init__(self, path):
        self.connection = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)
        self.connection.execute('PRAGMA query_only=ON')
        self.connection.execute('BEGIN')
        self.pending = []


_WORKER_INDEX = None
_WORKER_LINKS = None


def close_worker_indexes():
    global _WORKER_INDEX, _WORKER_LINKS
    for resource in (_WORKER_INDEX, _WORKER_LINKS):
        if resource is not None:
            resource.close()
    _WORKER_INDEX = _WORKER_LINKS = None


def init_unfold_worker(work_dir, stats, has_links):
    global _WORKER_INDEX, _WORKER_LINKS
    work_dir = Path(work_dir)
    _WORKER_INDEX = ReadOnlySegmentIndex(work_dir, stats)
    _WORKER_LINKS = ReadOnlyLinkStore(work_dir / 'links.sqlite') if has_links else None
    Finalize(None, close_worker_indexes, exitpriority=10)


def iter_w_tokens(walk: str) -> Iterator[Tuple[str, str]]:
    cursor = 0
    found = 0
    for match in WALK_TOKEN_RE.finditer(walk):
        if walk[cursor:match.start()].strip(" ,\t\r\n"):
            raise ValueError("Malformed GFA W walk")
        yield match.group(1), match.group(2)
        cursor = match.end()
        found += 1
    if walk[cursor:].strip(" ,\t\r\n") or found == 0:
        raise ValueError("Malformed or empty GFA W walk")


def iter_p_tokens(path_field: str) -> Iterator[Tuple[str, str]]:
    cursor = 0
    found = 0
    for match in P_TOKEN_RE.finditer(path_field):
        separator = path_field[cursor:match.start()]
        if ";" in separator:
            raise ValueError("GFA P jump separators are unsupported for coordinate unfolding")
        if separator.strip(" ,\t\r\n"):
            raise ValueError("Malformed GFA P/O path")
        yield (">" if match.group(2) == "+" else "<"), match.group(1)
        cursor = match.end()
        found += 1
    if ";" in path_field[cursor:]:
        raise ValueError("GFA P jump separators are unsupported for coordinate unfolding")
    if path_field[cursor:].strip(" ,\t\r\n") or found == 0:
        raise ValueError("Malformed or empty GFA P/O path")


@dataclass
class OutputState:
    genome_id: str
    contig_id: str
    gfa_path_id: str
    final_path: Path
    partial_path: Path
    source_lines: List[int] = field(default_factory=list)
    node_rows: int = 0
    span_bp: int = 0
    fragments: int = 0
    end0: Optional[int] = None
    digest: object = field(default_factory=hashlib.sha256)


def write_rows(state: OutputState, rows: Iterable[Tuple[str, int, int, str]], compression: str) -> None:
    first = state.fragments == 0
    if compression == "gzip":
        raw_handle = open(state.partial_path, "wb" if first else "ab")
        gz_handle = gzip.GzipFile(fileobj=raw_handle, mode="wb" if first else "ab", compresslevel=6, mtime=0)
        handle = io.TextIOWrapper(gz_handle, encoding="ascii", newline="")
    else:
        handle = open(state.partial_path, "wt" if first else "at", encoding="ascii", newline="")
    try:
        for contig, start0, end0, token in rows:
            line = f"{contig}\t{start0}\t{end0}\t{token}\n"
            handle.write(line)
            state.digest.update(line.encode("ascii"))
            state.node_rows += 1
            state.span_bp += end0 - start0
    finally:
        handle.close()
    state.fragments += 1


def rows_for_tokens(contig_id: str, start0: int, tokens: Iterable[Tuple[str, str]], index: SegmentIndex, links: Optional[mapping.LinkStore] = None) -> Iterable[Tuple[str, int, int, str]]:
    current = int(start0)
    previous = None
    for orientation, segment_name in tokens:
        internal_id, length = index.resolve(segment_name)
        if length <= 0:
            raise ValueError(f"Non-positive GFA segment length: {segment_name!r}")
        token = (segment_name, "+" if orientation == ">" else "-")
        if links is not None and previous is not None:
            status = links.status(mapping.canonical_link_key(*previous, *token))
            if status != "zero":
                raise ValueError(f"GFA P unspecified overlap: {status} for {previous!r} -> {token!r}")
        previous = token
        yield contig_id, current, current + length, f"{orientation}{internal_id}"
        current += length


def load_confirmed_mapping(contigs: str, inventory: str) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    contig_rows = mapping.read_tsv(contigs, mapping.CONTIG_COLUMNS)
    inventory_rows = mapping.read_tsv(inventory, mapping.GFA_INVENTORY_COLUMNS)
    inventory_by_id = {row["gfa_path_id"]: row for row in inventory_rows}
    selected = {}
    for row in contig_rows:
        if row["mapping_status"] != "confirmed":
            raise ValueError(
                f"Formal graph unfolding requires mapping_status=confirmed for every row; "
                f"line {row['_line_number']} is {row['mapping_status']!r}"
            )
        path_id = row["gfa_path_id"]
        if path_id not in inventory_by_id:
            raise ValueError(f"Unknown gfa_path_id in contig mapping: {path_id}")
        if inventory_by_id[path_id]["overlap_status"] not in mapping.ACCEPTED_PATH_STATUSES:
            raise ValueError(
                f"Selected GFA path cannot be unfolded without unsupported coordinate assumptions: "
                f"{path_id} ({inventory_by_id[path_id]['overlap_status']})"
            )
        if path_id in selected:
            raise ValueError(f"gfa_path_id is mapped more than once: {path_id}")
        selected[path_id] = row
    return selected, inventory_by_id


def build_index_and_spool(gfa: str, metadata: Dict[str, object], selected: Dict[str, Dict[str, str]], work_dir: Path, index: SegmentIndex, links: Optional[mapping.LinkStore], progress: mapping.ProgressLog) -> Dict[str, int]:
    spool_dir = work_dir / "selected_paths"
    spool_dir.mkdir()
    counters = dict(source_read_passes=0, source_lines=0, selected_records_spooled=0, links_indexed=0)
    progress.emit("index", "start", segment_index_builds=1, link_index_builds=int(links is not None))
    mapping.check_source_fingerprint(gfa, metadata["run_source_fingerprint"])
    with open_text(gfa, "rt") as handle:
        counters["source_read_passes"] += 1
        for line_number, raw in enumerate(handle, start=1):
            counters["source_lines"] = line_number
            if raw.startswith("S\t"):
                index.add(raw.rstrip("\r\n").split("\t"), metadata["gfa_major_version"], line_number)
            elif raw.startswith(("W\t", "P\t")):
                fields = raw.rstrip("\r\n").split("\t", 6)
                if fields[0] == "W":
                    if len(fields) < 7:
                        raise ValueError(f"Malformed GFA W record at line {line_number}")
                    path_id = mapping.stable_path_id("W", fields[1:4])
                else:
                    if len(fields) < 4:
                        raise ValueError(f"Malformed GFA P record at line {line_number}")
                    path_id = mapping.stable_path_id("P", [fields[1]])
                if path_id in selected:
                    with (spool_dir / path_id).open("at", encoding="utf-8", newline="") as spool:
                        spool.write(f"{line_number}\t")
                        spool.write(raw)
                        if not raw.endswith("\n"):
                            spool.write("\n")
                    counters["selected_records_spooled"] += 1
            elif links is not None and raw.startswith("L\t"):
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 6:
                    raise ValueError(f"Malformed GFA L record at line {line_number}")
                try:
                    links.add(mapping.canonical_link_key(fields[1], fields[2], fields[3], fields[4]), fields[5])
                except ValueError as exc:
                    raise ValueError(f"GFA link error at line {line_number}: {exc}") from exc
                counters["links_indexed"] += 1
            progress.update("index", line_number, n_segments=index.n_segments, selected_records_spooled=counters["selected_records_spooled"], links_indexed=counters["links_indexed"])
    mapping.check_source_fingerprint(gfa, metadata["run_source_fingerprint"])
    index.finish(metadata["node_stats"])
    if links is not None:
        links.flush()
        links.connection.execute("BEGIN")
    progress.emit("index", "complete", n_segments=index.n_segments, **counters)
    return counters


def unfold_spooled_path(state: OutputState, spool_path: Path, index: SegmentIndex, links: Optional[mapping.LinkStore], compression: str, progress: mapping.ProgressLog) -> None:
    if not spool_path.is_file():
        raise ValueError(f"Confirmed mapping refers to path not observed during unfolding: {state.gfa_path_id}")
    progress.emit("unfold", "path_start", gfa_path_id=state.gfa_path_id)
    with spool_path.open("rt", encoding="utf-8") as spool:
        for raw in spool:
            source_line, record = raw.split("\t", 1)
            line_number = int(source_line)
            fields = record.rstrip("\r\n").split("\t")
            try:
                old_span, old_count = state.span_bp, state.node_rows
                if fields[0] == "W":
                    if fields[4] == "*" or fields[5] == "*":
                        raise ValueError("Selected GFA W path has unknown coordinates")
                    start0, end0 = int(fields[4]), int(fields[5])
                    if start0 < 0 or end0 <= start0:
                        raise ValueError("Invalid GFA W coordinates")
                    if state.end0 is not None and state.end0 != start0:
                        raise ValueError(f"Non-contiguous W fragments for {state.gfa_path_id}")
                    write_rows(state, rows_for_tokens(state.contig_id, start0, iter_w_tokens(fields[6]), index), compression)
                    computed_end = start0 + state.span_bp - old_span
                    if computed_end != end0:
                        raise ValueError("GFA W span mismatch")
                else:
                    if state.fragments:
                        raise ValueError(f"Duplicate GFA P path {fields[1]!r}")
                    unspecified = fields[3] == "*"
                    single = mapping.P_TOKEN_RE.fullmatch(fields[2]) is not None
                    if unspecified and not single and links is None:
                        raise ValueError("GFA P path requires link verification but no link index was built")
                    write_rows(state, rows_for_tokens(state.contig_id, 0, iter_p_tokens(fields[2]), index,
                        links if unspecified and not single else None), compression)
                    count = state.node_rows - old_count
                    status = mapping.overlap_status(fields[3], count)
                    if status not in {"zero", "zero_single_segment", "requires_link_verification"}:
                        raise ValueError(f"Selected GFA P path has unsupported overlaps: {status}")
                    if status == "requires_link_verification" and links is None:
                        raise ValueError("GFA P path requires link verification but no link index was built")
                    end0 = state.span_bp - old_span
                state.end0 = end0
                state.source_lines.append(line_number)
                progress.emit("unfold", "fragment_complete", gfa_path_id=state.gfa_path_id, source_line=line_number, node_rows=state.node_rows, span_bp=state.span_bp, fragments=state.fragments)
            except ValueError as exc:
                raise ValueError(f"GFA path {state.gfa_path_id} at line {line_number}: {exc}") from exc
    progress.emit("unfold", "path_complete", gfa_path_id=state.gfa_path_id, node_rows=state.node_rows, span_bp=state.span_bp)


def validated_path_record(state, row, expected_length):
    path_id = state.gfa_path_id
    if row['path_length_bp'] and state.span_bp != int(row['path_length_bp']):
        raise ValueError(f"Unfolded span mismatch for {path_id}: {state.span_bp} versus {row['path_length_bp']}")
    if expected_length is not None and (state.span_bp != expected_length or state.end0 != expected_length):
        raise ValueError(f'HAL span mismatch for {path_id}: unfolded span={state.span_bp}, end={state.end0}, expected={expected_length}')
    source_lines = ','.join(map(str, state.source_lines))
    if source_lines != row['source_lines'] or state.fragments != int(row['fragment_count']):
        raise ValueError(f'GFA path inventory metadata mismatched for {path_id}; rerun inventory-gfa')
    return dict(genome_id=state.genome_id, contig_id=state.contig_id, gfa_path_id=path_id,
        node_rows=state.node_rows, span_bp=state.span_bp, sha256_uncompressed=state.digest.hexdigest(),
        output_file=state.final_path.name, source_lines=source_lines)


def unfold_path_task(task):
    genome_id, contig_id, path_id, final_path, work_dir, row, expected_length, compression = task
    final_path, work_dir = Path(final_path), Path(work_dir)
    state = OutputState(genome_id, contig_id, path_id, final_path, Path(str(final_path) + '.partial'))
    progress = mapping.ProgressLog(work_dir / 'path_progress' / (path_id + '.jsonl'))
    try:
        unfold_spooled_path(state, work_dir / 'selected_paths' / path_id, _WORKER_INDEX, _WORKER_LINKS, compression, progress)
        result = validated_path_record(state, row, expected_length)
        progress.emit('validation', 'path_complete', gfa_path_id=path_id, expected_length_bp=expected_length)
        return result
    except Exception as exc:
        progress.emit('unfold', 'error', gfa_path_id=path_id, message=str(exc))
        raise
    finally:
        progress.close()


def unfold_gfa(gfa: str, out_dir: str, *, contigs: str, gfa_inventory: str, compression: str = "gzip", expected_lengths: Optional[Mapping[str, int]] = None, threads: int = 1) -> Dict[str, object]:
    if type(threads) is not int or threads < 1:
        raise ValueError('Worker count must be a positive integer')
    if compression not in {"gzip", "none"}:
        raise ValueError(f"Unsupported pathBED compression: {compression}")
    selected, inventory_by_id = load_confirmed_mapping(contigs, gfa_inventory)
    metadata = load_inventory_metadata(gfa, gfa_inventory)
    if metadata["gfa_major_version"] != 1:
        raise ValueError("GFA2 paths are unsupported for coordinate unfolding")
    if expected_lengths is not None:
        if not isinstance(expected_lengths, Mapping) or set(expected_lengths) != set(selected):
            raise ValueError("Expected lengths must cover exactly the selected gfa_path_id values")
        if any(type(length) is not int or length <= 0 for length in expected_lengths.values()):
            raise ValueError("Expected lengths must be positive integer HAL lengths")
    output_dir = Path(out_dir).resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Graph unfolding output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "index"
    progress = mapping.ProgressLog(output_dir / "progress.jsonl")
    index = None
    links = None
    try:
        progress.emit("unfold", "start", n_confirmed_paths=len(selected))
        states: Dict[str, OutputState] = {}
        suffix = ".bed.gz" if compression == "gzip" else ".bed"
        filenames = set()
        for path_id, row in selected.items():
            safe_component(path_id, "gfa_path_id")
            genome_id = safe_component(row["genome_id"], "genome_id")
            contig_id = safe_component(row["contig_id"], "contig_id")
            filename = f"{genome_id}__{contig_id}{suffix}"
            if filename in filenames:
                raise ValueError(f"Duplicate pathBED output filename: {filename}")
            filenames.add(filename)
            final_path = output_dir / filename
            states[path_id] = OutputState(genome_id, contig_id, path_id, final_path, Path(str(final_path) + ".partial"))
        index = SegmentIndex(work_dir, metadata["node_stats"])
        if any(inventory_by_id[path_id]["record_type"] == "P" and inventory_by_id[path_id]["overlap_status"] in {"requires_link_verification", "zero_verified_from_links"} for path_id in selected):
            links = mapping.LinkStore(work_dir / "links.sqlite")
        counters = build_index_and_spool(gfa, metadata, selected, work_dir, index, links, progress)
        inventory_rows = []
        effective_workers = min(threads, len(states))
        if effective_workers <= 1:
            for path_id in sorted(states):
                state = states[path_id]
                unfold_spooled_path(state, work_dir / "selected_paths" / path_id, index, links, compression, progress)
                expected = expected_lengths[path_id] if expected_lengths is not None else None
                inventory_rows.append(validated_path_record(state, inventory_by_id[path_id], expected))
                progress.emit("validation", "path_complete", gfa_path_id=path_id, expected_length_bp=expected,
                    paths_completed=len(inventory_rows), paths_total=len(states))
        else:
            from pancgi_parallel import completed_tasks
            if index.mode == 'dense_numeric':
                with (work_dir / 'segment_lengths.bin').open('xb') as handle:
                    index.lengths.tofile(handle)
                index.lengths = array('I')
            def tasks():
                for path_id in sorted(states):
                    state = states[path_id]
                    progress.emit('unfold', 'path_queued', gfa_path_id=path_id)
                    yield (state.genome_id, state.contig_id, path_id, str(state.final_path), str(work_dir),
                        inventory_by_id[path_id], expected_lengths[path_id] if expected_lengths is not None else None, compression)
            results = completed_tasks(unfold_path_task, tasks(), effective_workers, init_unfold_worker,
                (str(work_dir), metadata['node_stats'], links is not None))
            try:
                for _, result in results:
                    inventory_rows.append(result)
                    path_id = result['gfa_path_id']
                    progress.emit('unfold', 'path_complete', gfa_path_id=path_id,
                        node_rows=result['node_rows'], span_bp=result['span_bp'])
                    progress.emit('validation', 'path_complete', gfa_path_id=path_id,
                        expected_length_bp=expected_lengths[path_id] if expected_lengths is not None else None,
                        paths_completed=len(inventory_rows), paths_total=len(states))
            finally:
                results.close()
            inventory_rows.sort(key=lambda row: row['gfa_path_id'])
        if len(inventory_rows) != len(states) or {row['gfa_path_id'] for row in inventory_rows} != set(states):
            raise RuntimeError('Incomplete unfolded path set')
        mapping.check_source_fingerprint(gfa, metadata["run_source_fingerprint"])
        progress.emit("validation", "complete", n_paths=len(states))
        for state in states.values():
            os.replace(state.partial_path, state.final_path)
        mapping.write_tsv(output_dir / "pathbed_inventory.tsv", inventory_rows, ["genome_id", "contig_id", "gfa_path_id", "node_rows", "span_bp", "sha256_uncompressed", "output_file", "source_lines"])
        summary = {
            "implementation": "PanCGI generic confirmed-mapping graph unfolding",
            "gfa": str(Path(gfa).resolve()), "gfa_major_version": metadata["gfa_major_version"],
            "inventory_source_fingerprint": metadata["source_fingerprint"],
            "run_source_fingerprint": metadata["run_source_fingerprint"],
            "segment_index_mode": index.mode, "n_segments": index.n_segments,
            "n_confirmed_paths": len(states), "n_pathbed_rows": sum(row['node_rows'] for row in inventory_rows),
            "pathbed_span_bp": sum(row['span_bp'] for row in inventory_rows), "mapping_status_required": "confirmed",
            "workers_requested": threads, "workers_used": effective_workers,
            "pathbed_inventory": str(output_dir / "pathbed_inventory.tsv"),
            "segment_index_builds": 1, "link_index_builds": int(links is not None),
            "inventory_node_stats_reused": True, "inventory_source_read_passes": metadata["source_read_passes"],
            "selected_spool_read_passes": 1, "expected_lengths_checked": expected_lengths is not None,
            "progress_file": str(output_dir / "progress.jsonl"), **counters,
        }
        mapping.write_json(output_dir / "unfold_summary.json", summary)
        progress.emit("unfold", "complete", n_paths=len(states), segment_index_builds=1, **counters)
        return summary
    except Exception as exc:
        progress.emit("unfold", "error", message=str(exc))
        raise
    finally:
        if index is not None:
            index.close()
        if links is not None:
            links.close()
        progress.close()


def run_from_args(args: argparse.Namespace) -> None:
    expected_lengths = None
    if args.expected_lengths:
        expected_lengths = json.loads(Path(args.expected_lengths).read_text(encoding="utf-8"))
    print(json.dumps(unfold_gfa(args.gfa, args.out_dir, contigs=args.contigs, gfa_inventory=args.gfa_inventory, compression=args.compression, expected_lengths=expected_lengths, threads=args.threads), indent=2, sort_keys=True))


def add_cli_parser(subparsers) -> argparse.ArgumentParser:
    parser = subparsers.add_parser("unfold-graph", help="Expand confirmed GFA paths into pathBED files")
    parser.add_argument('--threads', type=int, default=1, help='Concurrent path workers sharing one node index')
    parser.add_argument("--gfa", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--gfa-inventory", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--compression", choices=["gzip", "none"], default="gzip")
    parser.add_argument("--expected-lengths", help="JSON object mapping selected gfa_path_id values to HAL lengths")
    parser.set_defaults(func=run_from_args)
    return parser


def main() -> None:
    parser = argparse.ArgumentParser(description="PanCGI confirmed-mapping graph path unfolding")
    parser.add_argument('--threads', type=int, default=1, help='Concurrent path workers sharing one node index')
    parser.add_argument("--gfa", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--gfa-inventory", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--compression", choices=["gzip", "none"], default="gzip")
    parser.add_argument("--expected-lengths", help="JSON object mapping selected gfa_path_id values to HAL lengths")
    args = parser.parse_args()
    run_from_args(args)


if __name__ == "__main__":
    main()
