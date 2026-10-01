from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cpgi_nr_prod as prod
import pancgi_mapping as mapping


PATHBED_COLUMNS = [
    "genome_id",
    "contig_id",
    "gfa_path_id",
    "node_rows",
    "span_bp",
    "sha256_uncompressed",
    "output_file",
    "source_lines",
]


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_exact_tsv(path: Path, columns: Sequence[str]) -> List[Dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if list(reader.fieldnames or []) != list(columns):
            raise ValueError(f"Unexpected columns in {path}: {reader.fieldnames}")
        rows = [{key: str(value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path} contains no data rows")
    return rows


def validate_catalog(catalog: str, genomes: Sequence[Dict[str, str]], pathbed_dir: str, require_cpgi_fasta: bool) -> List[Dict[str, str]]:
    rows = list(prod.iter_internal_catalog(catalog))
    genome_by_id = {row["genome_id"]: row for row in genomes}
    if len(genome_by_id) != len(genomes):
        raise ValueError("Validated genome mapping contains duplicate genome_id values")
    if len({row["label"] for row in rows}) != len(rows):
        raise ValueError("Internal catalog contains duplicate labels")
    if {row["label"] for row in rows} != set(genome_by_id):
        raise ValueError("Internal catalog and validated genome mapping contain different genome_id values")
    expected_pathbed = os.path.realpath(pathbed_dir)
    for row in rows:
        genome = genome_by_id[row["label"]]
        if row["graph_sample"] != row["label"] or row["graph_hap"] != "0":
            raise ValueError(f"Internal catalog key is not normalized for {row['label']}")
        if row["hal_genome"] != genome["hal_genome"]:
            raise ValueError(f"Internal catalog HAL genome mismatch for {row['label']}")
        if os.path.realpath(row["path_bed_dir"]) != expected_pathbed:
            raise ValueError(f"Internal catalog pathBED directory mismatch for {row['label']}")
        if not Path(row["bed"]).is_file() or not Path(row["sv_tsv"]).is_file():
            raise FileNotFoundError(f"Prepared BED or SV input is absent for {row['label']}")
        if require_cpgi_fasta and not Path(row["cpgi_fa"]).is_file():
            raise FileNotFoundError(row["cpgi_fa"])
    kinds = [row["kind"] for row in rows]
    if kinds.count("reference_primary") != 1 or kinds.count("reference_comparison") > 1:
        raise ValueError("Internal catalog must contain one primary reference and at most one comparison reference")
    return rows


def validate_mapping_rows(genomes: Sequence[Dict[str, str]], contigs: Sequence[Dict[str, str]]) -> Dict[Tuple[str, str], Dict[str, str]]:
    genome_ids = {row["genome_id"] for row in genomes}
    result: Dict[Tuple[str, str], Dict[str, str]] = {}
    seen_gfa = set()
    for row in contigs:
        if row["mapping_status"] != "confirmed":
            raise ValueError(f"Unconfirmed contig mapping at input line {row['_line_number']}")
        if row["genome_id"] not in genome_ids:
            raise ValueError(f"Unknown genome_id in contig mapping: {row['genome_id']}")
        key = (row["genome_id"], row["contig_id"])
        if key in result:
            raise ValueError(f"Duplicate mapped contig: {key}")
        if row["gfa_path_id"] in seen_gfa:
            raise ValueError(f"GFA path is assigned more than once: {row['gfa_path_id']}")
        seen_gfa.add(row["gfa_path_id"])
        result[key] = row
    return result


def validate_pathbed(pathbed_dir: str, mapped: Dict[Tuple[str, str], Dict[str, str]], rows) -> Dict[str, int]:
    root = Path(pathbed_dir)
    observed = {}
    total_rows = 0
    total_span = 0
    for row in rows:
        key = (row["genome_id"], row["contig_id"])
        if key in observed:
            raise ValueError(f"Duplicate pathBED inventory key: {key}")
        if key not in mapped or row["gfa_path_id"] != mapped[key]["gfa_path_id"]:
            raise ValueError(f"pathBED inventory does not match confirmed mapping: {key}")
        output = root / row["output_file"]
        if not output.is_file() or output.stat().st_size == 0:
            raise FileNotFoundError(output)
        if len(row["sha256_uncompressed"]) != 64:
            raise ValueError(f"Invalid pathBED SHA-256 for {key}")
        node_rows = int(row["node_rows"])
        span_bp = int(row["span_bp"])
        if node_rows <= 0 or span_bp <= 0:
            raise ValueError(f"Invalid pathBED counts for {key}")
        observed[key] = row
        total_rows += node_rows
        total_span += span_bp
    if set(observed) != set(mapped):
        missing = sorted(set(mapped) - set(observed))
        extra = sorted(set(observed) - set(mapped))
        raise ValueError(f"pathBED/confirmed-mapping mismatch; missing={missing[:10]}, extra={extra[:10]}")
    summary = json.loads((root / "unfold_summary.json").read_text(encoding="utf-8"))
    if summary.get("implementation") != "PanCGI generic confirmed-mapping graph unfolding":
        raise ValueError("Unrecognized graph-unfold implementation")
    if int(summary.get("n_confirmed_paths", -1)) != len(observed):
        raise ValueError("Graph-unfold summary path count does not match inventory")
    if int(summary.get("n_pathbed_rows", -1)) != total_rows or int(summary.get("pathbed_span_bp", -1)) != total_span:
        raise ValueError("Graph-unfold summary totals do not match inventory")
    return {"n_pathbed_files": len(rows), "n_pathbed_rows": total_rows, "pathbed_span_bp": total_span}


def validate(args: argparse.Namespace) -> Dict[str, object]:
    genomes = mapping.read_tsv(args.genomes, mapping.GENOME_COLUMNS)
    contigs = mapping.read_tsv(args.contigs, mapping.CONTIG_COLUMNS)
    mapped = validate_mapping_rows(genomes, contigs)
    catalog = validate_catalog(args.catalog, genomes, args.pathbed_dir, args.require_cpgi_fasta)
    report: Dict[str, object] = {
        "status": "pass",
        "n_genomes": len(genomes),
        "n_contigs": len(contigs),
        "n_catalog_rows": len(catalog),
        "catalog_sha256": sha256_file(args.catalog),
        "mapping_status_required": "confirmed",
    }
    observed = read_exact_tsv(Path(args.pathbed_dir) / 'pathbed_inventory.tsv', PATHBED_COLUMNS)
    report.update(validate_pathbed(args.pathbed_dir, mapped, observed))
    expected_file = Path(args.contigs).parent / 'path_lengths.json'
    expected = json.loads(expected_file.read_text())
    if set(expected) != {r['gfa_path_id'] for r in observed}:
        raise ValueError('Unfolded paths differ from expected HAL sequences')
    if any(int(r['span_bp']) != expected[r['gfa_path_id']] for r in observed):
        raise ValueError('Unfolded path length differs from HAL sequence length')
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate prepared PanCGI inputs")
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--genomes", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--pathbed-dir", required=True)
    parser.add_argument("--require-cpgi-fasta", action="store_true")
    args = parser.parse_args()
    print(json.dumps(validate(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
