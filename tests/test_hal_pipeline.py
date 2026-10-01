from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pancgi_hal as hal


def test_prepare_bed4() -> None:
    with tempfile.TemporaryDirectory(prefix="pancgi_hal_prepare_") as directory:
        root = Path(directory)
        bed = root / "input.bed"
        out = root / "output.bed4"
        bed.write_text("ASM#1#ctgA\t10\t30\tCpG: 2\t20\t2\t12\t20\t60\t1\n", encoding="ascii")
        row = {
            "label": "sample_genome",
            "kind": "assembly",
            "graph_sample": "sample_genome",
            "graph_hap": "0",
            "bed": str(bed),
        }
        bed.write_text("sample_genome#0#contig_1\t10\t30\tCpG: 2\t20\t2\t12\t20\t60\t1\n", encoding="ascii")
        identifiers, digest = hal.prepare_bed4(row, out, {("sample_genome", "contig_1"): "hal_sequence_A"})
        assert identifiers == {"sample_genome#0#contig_1:10-30"}
        assert len(digest) == 64
        assert out.read_text(encoding="ascii") == "hal_sequence_A\t10\t30\tsample_genome#0#contig_1:10-30\n"


def test_parse_psl() -> None:
    with tempfile.TemporaryDirectory(prefix="pancgi_hal_psl_") as directory:
        path = Path(directory) / "output.psl"
        rows = [
            ["ASM#1#ctgA:10-30", "20", "0", "0", "0", "0", "0", "0", "0", "++", "ctgA", "100", "10", "30", "chr1", "1000", "100", "120", "1", "20,", "10,", "100,"],
            ["ASM#1#ctgB:20-40", "20", "0", "0", "0", "0", "0", "0", "0", "++", "ctgB", "100", "20", "40", "chr2", "1000", "200", "220", "1", "20,", "20,", "200,"],
            ["ASM#1#ctgB:20-40", "20", "0", "0", "0", "0", "0", "0", "0", "++", "ctgB", "100", "20", "40", "chr3", "1000", "300", "320", "1", "20,", "20,", "300,"],
        ]
        with path.open("wt", encoding="ascii", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(rows)
        records, counts, unknown = hal.parse_psl(path, {"ASM#1#ctgA:10-30", "ASM#1#ctgB:20-40"})
        assert records == 3
        assert counts["ASM#1#ctgA:10-30"] == 1
        assert counts["ASM#1#ctgB:20-40"] == 2
        assert unknown == []


if __name__ == "__main__":
    test_prepare_bed4()
    test_parse_psl()
    print("OK: PanCGI HAL pipeline tests passed")
