from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, TextIO, Tuple


def open_text(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def inspect_image(docker_bin: str, image: str) -> str:
    mode = os.environ.get('PANCGI_HAL_RUNTIME', 'docker')
    if mode == 'native':
        if image:
            raise ValueError('Native HAL runtime cannot select a Docker image')
        return 'native'
    if mode != 'docker' or not image:
        raise ValueError('Select native or an explicit Docker HAL image')
    result = subprocess.run(
        [docker_bin, "image", "inspect", image, "--format", "{{.Id}}"],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        raise RuntimeError(f"Docker image is unavailable: {image}\n{result.stderr.strip()}")
    return result.stdout.strip()


def container_command(
    docker_bin: str,
    image: str,
    mounts: Sequence[Tuple[Path, str]],
    executable: str,
    arguments: Sequence[str],
) -> List[str]:
    mode = os.environ.get('PANCGI_HAL_RUNTIME', 'docker')
    if mode == 'native':
        if image:
            raise ValueError('Native HAL runtime cannot select a Docker image')
        program = shutil.which(executable)
        if program is None:
            raise FileNotFoundError(f'Explicit native HAL executable not found: {executable}')
        return [program, *map(str, arguments)]
    if mode != 'docker' or not image:
        raise ValueError('Select native or an explicit Docker HAL image')
    memory = os.environ.get('PANCGI_HAL_MEMORY', '1g')
    cpus = os.environ.get('PANCGI_HAL_CPUS', '1')
    if not cpus.isdecimal() or int(cpus) < 1:
        raise ValueError('PANCGI_HAL_CPUS must be a positive integer')
    if not re.fullmatch(r'[1-9][0-9]*[mMgG]', memory):
        raise ValueError('PANCGI_HAL_MEMORY must be positive MiB/GiB, for example 512m or 2g')
    command = [docker_bin, 'run', '--rm', '--pull', 'never', '--network', 'none',
               '--cpus', cpus, '--memory', memory, '--memory-swap', memory, '--pids-limit', '128',
               '--user', f'{os.getuid()}:{os.getgid()}', '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges']
    seen = set()
    for directory, mode in mounts:
        resolved = directory.resolve()
        key = (str(resolved), mode)
        if key in seen:
            continue
        seen.add(key)
        command.extend(["-v", f"{resolved}:{resolved}:{mode}"])
    command.extend(['--entrypoint', executable, image])
    command.extend(str(value) for value in arguments)
    return command


def hal_genomes(
    hal: str,
    docker_bin: str,
    image: str,
    hal_stats: str,
) -> List[str]:
    hal_path = Path(hal).resolve()
    command = container_command(
        docker_bin,
        image,
        [(hal_path.parent, "ro")],
        hal_stats,
        ["--genomes", str(hal_path)],
    )
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"halStats --genomes failed: {result.stderr.strip()}")
    genomes = sorted(set(result.stdout.split()))
    if not genomes:
        raise RuntimeError("halStats --genomes returned no genomes")
    return genomes


def hal_sequence_stats(
    hal: str,
    genome: str,
    docker_bin: str,
    image: str,
    hal_stats: str,
    allow_empty: bool = False,
) -> List[Dict[str, object]]:
    hal_path = Path(hal).resolve()
    command = container_command(
        docker_bin,
        image,
        [(hal_path.parent, "ro")],
        hal_stats,
        [str(hal_path), "--sequenceStats", genome],
    )
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"halStats --sequenceStats failed for {genome}: {result.stderr.strip()}"
        )
    rows: List[Dict[str, object]] = []
    reader = csv.reader(result.stdout.splitlines(), skipinitialspace=True)
    for line_number, fields in enumerate(reader, start=1):
        if not fields:
            continue
        normalized = [field.strip() for field in fields]
        if normalized[0] == "SequenceName":
            continue
        if len(normalized) != 4:
            raise RuntimeError(
                f"Unexpected halStats --sequenceStats row for {genome} at output line {line_number}: {normalized}"
            )
        rows.append(
            {
                "hal_genome": genome,
                "hal_sequence": normalized[0],
                "length_bp": int(normalized[1]),
                "top_segments": int(normalized[2]),
                "bottom_segments": int(normalized[3]),
            }
        )
    if not rows and not allow_empty:
        raise RuntimeError(f"halStats --sequenceStats returned no sequences for {genome}")
    return rows


def parse_full_chrom(value: str) -> Tuple[str, str, str]:
    parts = str(value).split("#", 2)
    if len(parts) != 3 or any(not part for part in parts):
        raise ValueError(
            "Internal CpGI BED column 1 must contain genome_id#0#contig_id"
        )
    return parts[0], parts[1], parts[2]


def read_cpgi_intervals(path: str, genome_id: str, contig_index: Dict[Tuple[str, str], str]) -> List[Tuple[str, int, int, str]]:
    records: List[Tuple[str, int, int, str]] = []
    identifiers = set()
    with open_text(path, "rt") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise ValueError(f"{path}: CpGI BED line {line_number} must have 10 columns")
            sample, hap, contig = parse_full_chrom(fields[0])
            if sample != genome_id or hap != "0":
                raise ValueError(f"{path}: internal CpGI key does not match genome_id={genome_id}: {fields[0]}")
            key = (sample, contig)
            if key not in contig_index:
                raise ValueError(f"{path}: internal contig has no confirmed HAL mapping: {key}")
            start0 = int(fields[1])
            end0 = int(fields[2])
            if start0 < 0 or end0 <= start0:
                raise ValueError(f"{path}: invalid interval at line {line_number}")
            fid = f"{sample}#{hap}#{contig}:{start0}-{end0}"
            if fid in identifiers:
                raise ValueError(f"{path}: duplicate CpGI identifier {fid}")
            identifiers.add(fid)
            records.append((contig_index[key], start0, end0, fid))
    return records


def extract_intervals_from_fasta_stream(
    stream: Iterable[str],
    records: Sequence[Tuple[str, int, int, str]],
) -> Tuple[Dict[str, str], Dict[str, int]]:
    by_contig: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for contig, start0, end0, fid in records:
        by_contig[str(contig)].append((int(start0), int(end0), str(fid)))
    for contig in by_contig:
        by_contig[contig].sort(key=lambda row: (row[0], row[1], row[2]))
    pieces: Dict[str, List[str]] = {fid: [] for _, _, _, fid in records}
    observed: Dict[str, int] = {fid: 0 for _, _, _, fid in records}
    current = ""
    offset = 0
    next_index = 0
    active: List[Tuple[int, int, str]] = []
    intervals: List[Tuple[int, int, str]] = []
    headers = set()

    def begin_contig(name: str) -> None:
        nonlocal current, offset, next_index, active, intervals
        current = name
        offset = 0
        next_index = 0
        active = []
        intervals = by_contig.get(name, [])

    for raw in stream:
        line = raw.rstrip("\r\n")
        if not line.strip():
            continue
        if line.lstrip().startswith(">"):
            name = line[1:]
            if (not line.startswith(">") or not name or name != name.strip()
                    or any(ord(char) < 32 or ord(char) == 127 for char in name)):
                raise ValueError(f"HAL FASTA output contains a malformed header: {line!r}")
            if name in headers:
                raise ValueError(f"HAL FASTA output contains a duplicate header: {name!r}")
            headers.add(name)
            begin_contig(name)
            continue
        if not current:
            raise ValueError("HAL FASTA output contains sequence before a header")
        chunk = line.strip().upper()
        chunk_start = offset
        chunk_end = offset + len(chunk)
        while next_index < len(intervals) and intervals[next_index][0] < chunk_end:
            active.append(intervals[next_index])
            next_index += 1
        retained: List[Tuple[int, int, str]] = []
        for start0, end0, fid in active:
            if end0 > chunk_start and start0 < chunk_end:
                left = max(start0, chunk_start) - chunk_start
                right = min(end0, chunk_end) - chunk_start
                sequence = chunk[left:right]
                pieces[fid].append(sequence)
                observed[fid] += len(sequence)
            if end0 > chunk_end:
                retained.append((start0, end0, fid))
        active = retained
        offset = chunk_end
    return {fid: "".join(parts) for fid, parts in pieces.items()}, observed


def write_fasta(path: Path, records: Sequence[Tuple[str, int, int, str]], sequences: Dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".partial.{os.getpid()}")
    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "wt"
    with opener(temporary, mode, encoding="ascii", newline="") as handle:
        for _, _, _, fid in records:
            sequence = sequences[fid]
            handle.write(f">{fid}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index:index + 60] + "\n")
    os.replace(temporary, path)


def extract_cpgi_fasta_from_hal(
    row: Dict[str, str],
    hal: str,
    docker_bin: str,
    image: str,
    hal2fasta: str,
    log_path: str,
    contig_index: Dict[Tuple[str, str], str],
) -> Dict[str, object]:
    hal_path = Path(hal).resolve()
    output = Path(row["cpgi_fa"]).resolve()
    log = Path(log_path).resolve()
    records = read_cpgi_intervals(row["bed"], row["graph_sample"], contig_index)
    if not records:
        write_fasta(output, records, {})
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text("Empty CGI input; no HAL FASTA extraction required.\n", encoding="utf-8")
        return {
            "label": row["label"],
            "records": 0,
            "missing": [],
            "output": str(output),
        }
    command = container_command(
        docker_bin,
        image,
        [(hal_path.parent, "ro")],
        hal2fasta,
        ["--upper", str(hal_path), row["hal_genome"]],
    )
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("wt", encoding="utf-8") as error_handle:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=error_handle,
            text=True,
            encoding="ascii",
            errors="strict",
            bufsize=1024 * 1024,
        )
        if process.stdout is None:
            process.kill()
            raise RuntimeError(f"{row['label']}: hal2fasta stdout is unavailable")
        try:
            sequences, observed = extract_intervals_from_fasta_stream(process.stdout, records)
        except Exception:
            process.kill()
            process.wait()
            raise
        return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"{row['label']}: hal2fasta failed with exit code {return_code}; see {log}")
    missing = []
    for _, start0, end0, fid in records:
        expected = end0 - start0
        actual = observed[fid]
        if actual != expected:
            missing.append((fid, expected, actual))
    if missing:
        return {
            "label": row["label"],
            "records": len(records),
            "missing": missing,
            "output": str(output),
        }
    write_fasta(output, records, sequences)
    return {
        "label": row["label"],
        "records": len(records),
        "missing": [],
        "output": str(output),
    }
