from __future__ import annotations

import gzip
import hashlib
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str((Path(__file__).resolve().parents[1] / "src" / "pancgi_app")))

import pancgi_graph_unfold as unfold
import pancgi_mapping as mapping


def read_text(path: Path) -> str:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="ascii") as handle:
            return handle.read()
    return path.read_text(encoding="ascii")


def write_contigs(path: Path, inventory: Path, statuses=None) -> None:
    rows = mapping.read_tsv(str(inventory), mapping.GFA_INVENTORY_COLUMNS)
    output = []
    for index, row in enumerate(rows):
        output.append({
            "genome_id": f"genome_{index + 1}",
            "contig_id": f"contig_{index + 1}",
            "gfa_path_id": row["gfa_path_id"],
            "hal_sequence": f"hal_{index + 1}",
            "cpgi_contig": f"bed_{index + 1}",
            "sv_contig": "",
            "mapping_status": statuses[index] if statuses else "confirmed",
            "mapping_method": "test",
            "mapping_note": "",
        })
    mapping.write_tsv(path, output, mapping.CONTIG_COLUMNS)


def expect_failure(function, text: str) -> None:
    try:
        function()
    except (ValueError, FileExistsError) as exc:
        assert text in str(exc), (text, str(exc))
    else:
        raise AssertionError(f"Expected failure containing {text!r}")


def test_w_paths_and_opaque_segments(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "w.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    out = root / "out"
    gfa.write_text(
        "H\tVN:Z:1.1\n"
        "S\tunitig-A\tAAA\n"
        "S\tsegment.2\tCC\n"
        "L\tunitig-A\t+\tsegment.2\t-\t0M\n"
        "W\tindividual name\tmaternal\tscaffold-A\t0\t5\t>unitig-A<segment.2\n",
        encoding="utf-8",
    )
    mapping.inventory_gfa(str(gfa), str(inventory), str(root / "inventory_work"))
    write_contigs(contigs, inventory)
    summary = unfold.unfold_gfa(str(gfa), str(out), contigs=str(contigs), gfa_inventory=str(inventory), compression="gzip")
    assert summary["segment_index_mode"] == "sqlite_opaque"
    assert read_text(out / "genome_1__contig_1.bed.gz") == "contig_1\t0\t3\t>0\ncontig_1\t3\t5\t<1\n"
    assert (out / "genome_1__contig_1.bed.gz").read_bytes() == bytes.fromhex(
        "1f8b08080000000000ff67656e6f6d655f315f5f636f6e7469675f312e6265642e677a2e7061727469616c00"
        "4acecf2bc94c8f37e434e034e6b433e04a86f18d394d396d0cb900000000ffff0300847437bd20000000"
    )
    assert (out / "index" / "segment_id_map.tsv.gz").is_file()


def test_p_path_name_is_not_parsed(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "p.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    out = root / "out"
    gfa.write_text(
        "H\tVN:Z:1.0\n"
        "S\t1\tAAA\n"
        "S\t2\tCC\n"
        "L\t1\t+\t2\t-\t0M\n"
        "P\topaque-path-without-delimiters\t1+,2-\t0M\n",
        encoding="utf-8",
    )
    mapping.inventory_gfa(str(gfa), str(inventory), str(root / "inventory_work"))
    write_contigs(contigs, inventory)
    summary = unfold.unfold_gfa(str(gfa), str(out), contigs=str(contigs), gfa_inventory=str(inventory), compression="none")
    assert summary["segment_index_mode"] == "dense_numeric"
    assert read_text(out / "genome_1__contig_1.bed") == "contig_1\t0\t3\t>1\ncontig_1\t3\t5\t<2\n"


def test_gfa2_o_path(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "o.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    out = root / "out"
    gfa.write_text(
        "H\tVN:Z:2.0\n"
        "S\ta\t3\tAAA\n"
        "S\tb\t2\tCC\n"
        "O\tarbitrary_name\ta+ b-\n",
        encoding="utf-8",
    )
    mapping.inventory_gfa(str(gfa), str(inventory), str(root / "inventory_work"))
    write_contigs(contigs, inventory)
    expect_failure(
        lambda: unfold.unfold_gfa(str(gfa), str(out), contigs=str(contigs), gfa_inventory=str(inventory), compression="none"),
        "unsupported coordinate assumptions",
    )


def test_confirmed_gate(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "p.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    gfa.write_text("S\t1\tAAA\nP\tpath\t1+\t*\n", encoding="utf-8")
    mapping.inventory_gfa(str(gfa), str(inventory), str(root / "inventory_work"))
    write_contigs(contigs, inventory, statuses=["candidate"])
    expect_failure(
        lambda: unfold.unfold_gfa(str(gfa), str(root / "out"), contigs=str(contigs), gfa_inventory=str(inventory)),
        "mapping_status=confirmed",
    )


def test_nonzero_overlap_fails(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "p.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    gfa.write_text("S\t1\tAAA\nS\t2\tCC\nP\tpath\t1+,2+\t1M\n", encoding="utf-8")
    mapping.inventory_gfa(str(gfa), str(inventory), str(root / "inventory_work"))
    write_contigs(contigs, inventory)
    expect_failure(
        lambda: unfold.unfold_gfa(str(gfa), str(root / "out"), contigs=str(contigs), gfa_inventory=str(inventory)),
        "unsupported coordinate assumptions",
    )


def test_unspecified_overlap_without_link_fails(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "p.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    gfa.write_text("S\t1\tAAA\nS\t2\tCC\nP\tpath\t1+,2+\t*\n", encoding="utf-8")
    mapping.inventory_gfa(str(gfa), str(inventory), str(root / "inventory_work"))
    write_contigs(contigs, inventory)
    expect_failure(
        lambda: unfold.unfold_gfa(str(gfa), str(root / "out"), contigs=str(contigs), gfa_inventory=str(inventory)),
        "missing_link",
    )


def test_deterministic_gzip(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "p.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    gfa.write_text("S\t1\tAAA\nP\tpath\t1+\t*\n", encoding="utf-8")
    mapping.inventory_gfa(str(gfa), str(inventory), str(root / "inventory_work"))
    write_contigs(contigs, inventory)
    for name in ["out1", "out2"]:
        unfold.unfold_gfa(str(gfa), str(root / name), contigs=str(contigs), gfa_inventory=str(inventory), compression="gzip")
    first = (root / "out1" / "genome_1__contig_1.bed.gz").read_bytes()
    second = (root / "out2" / "genome_1__contig_1.bed.gz").read_bytes()
    assert first == second
    assert hashlib.sha256(first).hexdigest() == hashlib.sha256(second).hexdigest()
    assert first == bytes.fromhex(
        "1f8b08080000000000ff67656e6f6d655f315f5f636f6e7469675f312e6265642e677a2e7061727469616c00"
        "4acecf2bc94c8f37e434e034e6b433e402000000ffff0300c8a0196910000000"
    )


def test_repeated_node_occurrences(root: Path) -> None:
    root.mkdir(parents=True)
    gfa = root / "w.gfa"
    inventory = root / "gfa.tsv"
    contigs = root / "contigs.tsv"
    gfa.write_text("S\t1\tAAA\nS\t2\tCC\nW\ts\t1\tc\t0\t8\t>1>2>1\n", encoding="ascii")
    mapping.inventory_gfa(str(gfa), str(inventory))
    write_contigs(contigs, inventory)
    path_id = mapping.stable_path_id("W", ["s", "1", "c"])
    out = root / "out"
    summary = unfold.unfold_gfa(str(gfa), str(out), contigs=str(contigs), gfa_inventory=str(inventory),
                                compression="none", expected_lengths={path_id: 8})
    assert summary["n_pathbed_rows"] == 3
    assert (out / "genome_1__contig_1.bed").read_bytes() == b"contig_1\t0\t3\t>1\ncontig_1\t3\t5\t>2\ncontig_1\t5\t8\t>1\n"


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="pancgi_graph_unfold_") as directory:
        root = Path(directory)
        tests = [
            test_w_paths_and_opaque_segments,
            test_p_path_name_is_not_parsed,
            test_gfa2_o_path,
            test_confirmed_gate,
            test_nonzero_overlap_fails,
            test_unspecified_overlap_without_link_fails,
            test_deterministic_gzip,
            test_repeated_node_occurrences,
        ]
        for index, test in enumerate(tests):
            test(root / f"case_{index}")
    print("OK: PanCGI graph unfolding tests passed")
