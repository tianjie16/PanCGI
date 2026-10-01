from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pancgi_graph_unfold as unfold
import pancgi_mapping as mapping


class StreamingGraphTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.gfa = self.root / "graph.gfa"
        self.inventory = self.root / "gfa.tsv"
        self.contigs = self.root / "contigs.tsv"
        self.out = self.root / "out"

    def tearDown(self):
        self.temporary.cleanup()

    def prepare(self, text, selected=None):
        self.gfa.write_text(text, encoding="utf-8")
        result = mapping.inventory_gfa(str(self.gfa), str(self.inventory), str(self.root / "inventory"))
        self.rows = mapping.read_tsv(str(self.inventory), mapping.GFA_INVENTORY_COLUMNS)
        rows = [row for row in self.rows if selected is None or row["sequence_id"] in selected]
        self.selected = rows
        mapping.write_tsv(self.contigs, [dict(genome_id=f"g{i}", contig_id=f"c{i}", gfa_path_id=row["gfa_path_id"],
                                            mapping_status="confirmed") for i, row in enumerate(rows)], mapping.CONTIG_COLUMNS)
        return result

    def run_unfold(self, **kwargs):
        return unfold.unfold_gfa(str(self.gfa), str(self.out), contigs=str(self.contigs),
                                 gfa_inventory=str(self.inventory), compression="none", **kwargs)

    def assert_no_outputs(self):
        self.assertFalse(list(self.out.glob("*.bed")))
        self.assertFalse(list(self.out.glob("*.bed.gz")))
        self.assertFalse((self.out / "pathbed_inventory.tsv").exists())
        self.assertFalse((self.out / "unfold_summary.json").exists())

    def test_inventory_is_one_metadata_pass_without_sql_or_lengths(self):
        text = "H\tVN:Z:1.1\nS\t1\tAAA\nS\t2\tCC\nW\ts\t1\tc\t0\t5\t>1>2\nP\tp\t1+,2+\t*\n"
        original = mapping.open_text
        opens = []

        def tracked(path, mode="rt"):
            if Path(path) == self.gfa:
                opens.append(path)
            return original(path, mode)

        with patch.object(mapping, "open_text", tracked), patch.object(mapping.sqlite3, "connect", side_effect=AssertionError("Inventory must not use SQL")), patch.object(mapping, "segment_length", side_effect=AssertionError("Inventory must not resolve lengths")):
            result = self.prepare(text)
        self.assertEqual(len(opens), 1)
        self.assertEqual(result["source_read_passes"], 1)
        self.assertEqual(result["segment_index_builds"], 0)
        self.assertEqual(result["link_index_builds"], 0)
        self.assertEqual(self.inventory.read_text().splitlines()[0].split("\t"), mapping.GFA_INVENTORY_COLUMNS)
        metadata = json.loads(mapping.inventory_metadata_path(str(self.inventory)).read_text())
        self.assertEqual(metadata["source_fingerprint"], mapping.source_fingerprint(str(self.gfa)))
        self.assertEqual(metadata["node_stats"], dict(n_segments=2, max_numeric=2, all_numeric=True, dense_numeric=True))
        w = next(row for row in self.rows if row["record_type"] == "W")
        p = next(row for row in self.rows if row["record_type"] == "P")
        self.assertEqual((w["path_length_bp"], w["source_lines"]), ("5", "4"))
        self.assertEqual((p["path_length_bp"], p["end0"], p["overlap_status"]), ("", "", "requires_link_verification"))

    def test_unfold_builds_one_index_and_reads_source_once(self):
        self.prepare("S\t1\tAAA\nS\t2\tCC\nW\ts\t1\tc\t0\t5\t>1<2\nP\tnot-selected\t1+,2+\t*\n", selected={"c"})
        with patch.object(unfold, "open_text", wraps=unfold.open_text) as opened, patch.object(unfold, "SegmentIndex", wraps=unfold.SegmentIndex) as indexed, patch.object(mapping, "LinkStore", side_effect=AssertionError("Unselected P must not build links")):
            result = self.run_unfold()
        self.assertEqual(opened.call_count, 1)
        self.assertEqual(indexed.call_count, 1)
        self.assertEqual(result["source_read_passes"], 1)
        self.assertEqual(result["segment_index_builds"], 1)
        self.assertEqual(result["link_index_builds"], 0)
        self.assertEqual(result["selected_records_spooled"], 1)
        spools = list((self.out / "index" / "selected_paths").iterdir())
        self.assertEqual(len(spools), 1)
        self.assertEqual(spools[0].read_text(), "3\tW\ts\t1\tc\t0\t5\t>1<2\n")
        self.assertEqual((self.out / "g0__c0.bed").read_bytes(), b"c0\t0\t3\t>1\nc0\t3\t5\t<2\n")

    def test_explicit_zero_and_w_do_not_query_topology(self):
        self.prepare("S\t1\tAAA\nS\t2\tCC\nL\t1\t+\t2\t+\t2M\nP\tp\t1+,2+\t0M\nW\ts\t1\tc\t0\t5\t>1>2\n")
        with patch.object(mapping, "LinkStore", side_effect=AssertionError("Explicit zero and W must not query links")):
            result = self.run_unfold()
        self.assertEqual(result["n_pathbed_rows"], 4)
        self.assertEqual(next(row for row in self.rows if row["record_type"] == "P")["overlap_status"], "zero")

    def test_w_errors_are_deferred_to_unfold(self):
        cases = [
            ("S\t1\tAAA\n", ">1>2", 5, "undefined segment"),
            ("S\t1\tAAA\n", ">1", 4, "span mismatch"),
            ("S\t1\tAAA\n", "invalid>1", 3, "Malformed GFA W walk"),
            ("S\t1\t*\n", ">1", 3, "no sequence or LN:i tag"),
            ("S\t1\t*\tLN:i:0\n", ">1", 3, "Sequence-bearing GFA segments are required"),
            ("S\t1\t*\tLN:i:-1\n", ">1", 3, "Sequence-bearing GFA segments are required"),
            ("S\t1\t\n", ">1", 3, "Invalid GFA segment"),
        ]
        for nodes, walk, end, message in cases:
            with self.subTest(message=message):
                self.out = self.root / f"out{end}-{len(walk)}-{len(nodes)}"
                self.prepare(nodes + f"W\ts\t1\tc\t0\t{end}\t{walk}\n")
                with self.assertRaisesRegex(ValueError, message):
                    self.run_unfold()
                self.assert_no_outputs()

    def test_duplicate_s_identifiers_still_reject(self):
        for i, name in enumerate(("1", "opaque-node")):
            with self.subTest(name=name):
                self.out = self.root / f"duplicate{i}"
                self.prepare(f"S\t{name}\tAAA\nS\t{name}\tCC\nW\ts\t1\tc\t0\t3\t>{name}\n")
                with self.assertRaisesRegex(ValueError, "Duplicate GFA segment ID"):
                    self.run_unfold()
                self.assert_no_outputs()

    def test_missing_sequence_with_declared_length(self):
        self.prepare("S\t1\t*\tLN:i:3\nW\ts\t1\tc\t0\t3\t>1\n")
        with self.assertRaisesRegex(ValueError, 'Sequence-bearing GFA segments'):
            self.run_unfold()
        self.assert_no_outputs()

    def test_repeated_occurrences_in_p_and_w_fragments(self):
        self.prepare("S\t1\tAAA\nS\t2\tCC\nP\tp\t1+,2+,1-\t0M,0M\nW\ts\t1\tc\t0\t5\t>1>2\nW\ts\t1\tc\t5\t8\t<1\n")
        result = self.run_unfold(expected_lengths={row["gfa_path_id"]: 8 for row in self.rows})
        self.assertEqual(result["n_pathbed_rows"], 6)
        for i in range(2):
            self.assertEqual((self.out / f"g{i}__c{i}.bed").read_text(), f"c{i}\t0\t3\t>1\nc{i}\t3\t5\t>2\nc{i}\t5\t8\t<1\n")

    def test_w_fragment_gaps_and_overlaps_fail_inventory(self):
        for start in (2, 4):
            with self.subTest(start=start), self.assertRaisesRegex(ValueError, "Non-contiguous W fragments"):
                self.prepare(f"S\t1\tAAA\nS\t2\tCC\nW\ts\t1\tc\t0\t3\t>1\nW\ts\t1\tc\t{start}\t{start+2}\t>2\n")

    def test_fragment_continuity_is_checked_during_unfold(self):
        self.prepare("S\t1\tAAA\nS\t2\tCC\nW\ts\t1\tc\t0\t3\t>1\nW\ts\t1\tc\t3\t5\t>2\n")
        original = unfold.build_index_and_spool

        def alter_spool(*args):
            result = original(*args)
            path = self.out / "index" / "selected_paths" / self.rows[0]["gfa_path_id"]
            path.write_text(path.read_text().replace("\t3\t5\t>2", "\t4\t6\t>2"))
            return result

        with patch.object(unfold, "build_index_and_spool", alter_spool), self.assertRaisesRegex(ValueError, "Non-contiguous W fragments"):
            self.run_unfold()
        self.assert_no_outputs()

    def test_p_unspecified_overlaps_reject_missing_nonzero_ambiguous_links(self):
        link = "L\t1\t+\t2\t+\t0M\n"
        cases = [("", "missing_link"), (link.replace("0M", "1M"), "link_overlap_unsupported"),
                 (link.replace("0M", "*"), "link_overlap_unsupported"),
                 (link + link, "multiple_matching_links"),
                 (link + "L\t2\t-\t1\t-\t0M\n", "multiple_matching_links"),
                 ("L\t1\t+\t2\t-\t0M\n", "missing_link")]
        for i, (records, message) in enumerate(cases):
            with self.subTest(links=records):
                self.out = self.root / f"links{i}"
                self.prepare("S\t1\tAAA\nS\t2\tCC\n" + records + "P\tp\t1+,2+\t*\n")
                self.assertEqual(self.rows[0]["overlap_status"], "requires_link_verification")
                with self.assertRaisesRegex(ValueError, message):
                    self.run_unfold()
                self.assert_no_outputs()

    def test_p_unspecified_reverse_link_is_resolved_in_same_source_pass(self):
        self.prepare("P\tp\t1+,2+\t*\nS\t1\tAAA\nS\t2\tCC\nL\t2\t-\t1\t-\t0M\n")
        with patch.object(unfold, "open_text", wraps=unfold.open_text) as opened, patch.object(mapping, "LinkStore", wraps=mapping.LinkStore) as linked:
            result = self.run_unfold()
        self.assertEqual(opened.call_count, 1)
        self.assertEqual(linked.call_count, 1)
        self.assertEqual(result["link_index_builds"], 1)
        self.assertEqual((self.out / "g0__c0.bed").read_text(), "c0\t0\t3\t>1\nc0\t3\t5\t>2\n")

    def test_singleton_p_unspecified_needs_no_links(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        with patch.object(mapping, "LinkStore", side_effect=AssertionError("Singleton needs no links")):
            self.run_unfold()

    def test_nonzero_and_jump_paths_reject(self):
        self.prepare("S\t1\tAAA\nS\t2\tCC\nP\tp\t1+,2+\t1M\n")
        with self.assertRaisesRegex(ValueError, "unsupported coordinate assumptions"):
            self.run_unfold()
        with self.assertRaisesRegex(ValueError, "jump separators"):
            self.prepare("S\t1\tAAA\nS\t2\tCC\nP\tp\t1+;2+\t*\n")
        with self.assertRaisesRegex(ValueError, "jump separators"):
            list(unfold.iter_p_tokens("1+;2+"))

    def test_numeric_sparse_and_opaque_ids_preserve_original_rules(self):
        for i, (names, mode, tokens) in enumerate([
            (("0", "3"), "dense_numeric", (">0", "<3")),
            (("100000", "200000"), "sqlite_opaque", (">0", "<1")),
            (("01", "unitig-A"), "sqlite_opaque", (">0", "<1")),
        ]):
            with self.subTest(names=names):
                self.out = self.root / f"ids{i}"
                self.prepare(f"S\t{names[0]}\tAAA\nS\t{names[1]}\tCC\nW\ts\t1\tc\t0\t5\t>{names[0]}<{names[1]}\n")
                self.assertEqual(self.run_unfold()["segment_index_mode"], mode)
                self.assertEqual((self.out / "g0__c0.bed").read_text(), f"c0\t0\t3\t{tokens[0]}\nc0\t3\t5\t{tokens[1]}\n")
                if mode == "sqlite_opaque":
                    with gzip.open(self.out / "index" / "segment_id_map.tsv.gz", "rt") as handle:
                        self.assertEqual(handle.read(), f"internal_node_id\tgfa_segment_id\tlength_bp\n0\t{names[0]}\t3\n1\t{names[1]}\t2\n")

    def test_source_change_before_unfold_fails_without_read(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        self.gfa.write_text(self.gfa.read_text() + "#changed\n")
        with patch.object(unfold, "open_text", side_effect=AssertionError("Must reject before opening source")), self.assertRaisesRegex(ValueError, "source changed"):
            self.run_unfold()

    def test_fingerprint_default_stays_strict_and_portable_checks_other_fields(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        actual = mapping.source_fingerprint(str(self.gfa))
        old = dict(actual, device=actual['device'] + 1000)
        with self.assertRaisesRegex(ValueError, 'source changed'):
            mapping.check_source_fingerprint(str(self.gfa), old)
        self.assertEqual(mapping.check_source_fingerprint(str(self.gfa), old, cross_node=True), actual)
        for key in ('path', 'inode', 'size', 'mtime_ns', 'ctime_ns'):
            altered = dict(old)
            altered[key] = altered[key] + ('x' if key == 'path' else 1)
            with self.subTest(field=key), self.assertRaisesRegex(ValueError, 'source changed'):
                mapping.check_source_fingerprint(str(self.gfa), altered, cross_node=True)
        for altered in (None, {}, dict(old, extra=1), dict(old, device='69'), {k: v for k, v in old.items() if k != 'device'}):
            with self.subTest(fingerprint=altered), self.assertRaisesRegex(ValueError, 'source changed'):
                mapping.check_source_fingerprint(str(self.gfa), altered, cross_node=True)

    def test_cross_node_unfold_preserves_inventory_and_pathbed(self):
        self.prepare("S\t1\tAAA\nS\t2\tCC\nW\ts\t1\tc\t0\t5\t>1<2\n")
        sidecar = mapping.inventory_metadata_path(str(self.inventory))
        metadata = json.loads(sidecar.read_text())
        metadata['source_fingerprint']['device'] += 1000
        sidecar.write_text(json.dumps(metadata))
        before = sidecar.read_bytes()
        result = self.run_unfold()
        self.assertEqual(result['source_read_passes'], 1)
        self.assertEqual(result['segment_index_builds'], 1)
        self.assertEqual(result['inventory_source_fingerprint'], metadata['source_fingerprint'])
        self.assertEqual(result['run_source_fingerprint'], mapping.source_fingerprint(str(self.gfa)))
        self.assertEqual(sidecar.read_bytes(), before)
        self.assertEqual((self.out / 'g0__c0.bed').read_bytes(), b'c0\t0\t3\t>1\nc0\t3\t5\t<2\n')

    def test_device_change_within_unfold_is_rejected_at_all_checkpoints(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        original = mapping.source_fingerprint
        for checkpoint in (2, 3, 4):
            with self.subTest(checkpoint=checkpoint):
                self.out = self.root / f'device_change_{checkpoint}'
                count = 0
                def changing(path):
                    nonlocal count
                    count += 1
                    state = original(path)
                    if count >= checkpoint:
                        state['device'] += 1
                    return state
                with patch.object(mapping, 'source_fingerprint', changing), \
                     self.assertRaisesRegex(ValueError, 'source changed'):
                    self.run_unfold()
                self.assertEqual(count, checkpoint)
                self.assert_no_outputs()

    def test_source_change_during_scan_rejects_before_publication(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        original = unfold.SegmentIndex.add

        def change_source(index, *args):
            original(index, *args)
            with self.gfa.open("at") as handle:
                handle.write("#changed\n")

        with patch.object(unfold.SegmentIndex, "add", change_source), self.assertRaisesRegex(ValueError, "source changed"):
            self.run_unfold()
        self.assert_no_outputs()

    def test_source_change_after_index_rejects_before_publication(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        original = unfold.write_rows

        def change_source(*args):
            original(*args)
            with self.gfa.open("at") as handle:
                handle.write("#changed\n")

        with patch.object(unfold, "write_rows", change_source), self.assertRaisesRegex(ValueError, "source changed"):
            self.run_unfold()
        self.assert_no_outputs()

    def test_missing_sidecar_does_not_trigger_fallback_scan(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        mapping.inventory_metadata_path(str(self.inventory)).unlink()
        with patch.object(unfold, "open_text", side_effect=AssertionError("No fallback scan")), self.assertRaisesRegex(ValueError, "Missing or invalid"):
            self.run_unfold()

    def test_mismatched_cached_statistics_fail(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        sidecar = mapping.inventory_metadata_path(str(self.inventory))
        metadata = json.loads(sidecar.read_text())
        metadata["node_stats"]["n_segments"] = 2
        sidecar.write_text(json.dumps(metadata))
        with self.assertRaisesRegex(ValueError, "statistics mismatched"):
            self.run_unfold()
        self.assert_no_outputs()

    def test_all_hal_lengths_checked_before_any_output_is_promoted(self):
        self.prepare("S\t1\tAAA\nS\t2\tCC\nP\tone\t1+\t*\nP\ttwo\t2+\t*\n")
        lengths = {row["gfa_path_id"]: (3 if row["raw_path_name"] == "one" else 2) for row in self.rows}
        lengths[self.rows[-1]["gfa_path_id"]] += 1
        with self.assertRaisesRegex(ValueError, "HAL span mismatch"):
            self.run_unfold(expected_lengths=lengths)
        self.assert_no_outputs()
        self.assertTrue((self.out / "g0__c0.bed.partial").exists())

    def test_hal_length_mapping_must_be_complete_and_integer(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        key = self.rows[0]["gfa_path_id"]
        for lengths in ({}, {key: 3, "extra": 3}, {key: True}, {key: 3.0}, {key: -1}):
            with self.subTest(lengths=lengths), self.assertRaisesRegex(ValueError, "Expected lengths"):
                self.run_unfold(expected_lengths=lengths)

    def test_progress_counters_during_indexing(self):
        self.prepare("S\t1\tAAA\n" + "#\n" * 100000 + "P\tp\t1+\t*\n")
        self.run_unfold()
        events = [json.loads(line) for line in (self.out / "progress.jsonl").read_text().splitlines()]
        self.assertTrue(any(event["phase"] == "index" and event["event"] == "progress" and event["source_lines"] == 100000 for event in events))
        self.assertTrue(any(event["event"] == "path_complete" for event in events))
        self.assertEqual(events[-1]["event"], "complete")
        self.assertEqual(events[-1]["source_read_passes"], 1)

    def test_gzip_source_is_read_once_and_spool_has_line_numbers(self):
        self.gfa = self.root / "graph.gfa.gz"
        with gzip.open(self.gfa, "wt") as handle:
            handle.write("S\t1\tAAA\nP\tp\t1+\t*\n")
        mapping.inventory_gfa(str(self.gfa), str(self.inventory))
        key = mapping.stable_path_id("P", ["p"])
        mapping.write_tsv(self.contigs, [dict(genome_id="g", contig_id="c", gfa_path_id=key, mapping_status="confirmed")], mapping.CONTIG_COLUMNS)
        with patch.object(unfold, "open_text", wraps=unfold.open_text) as opened:
            self.run_unfold()
        self.assertEqual(opened.call_count, 1)
        self.assertEqual((self.out / "index" / "selected_paths" / key).read_text(), "2\tP\tp\t1+\t*\n")

    def test_expected_lengths_cli(self):
        self.prepare("S\t1\tAAA\nP\tp\t1+\t*\n")
        expected = self.root / "lengths.json"
        expected.write_text(json.dumps({self.rows[0]["gfa_path_id"]: 3}))
        result = subprocess.run([sys.executable, "-B", str(Path(unfold.__file__)), "--gfa", str(self.gfa),
                                 "--contigs", str(self.contigs), "--gfa-inventory", str(self.inventory),
                                 "--out-dir", str(self.out), "--expected-lengths", str(expected), "--compression", "none"],
                                capture_output=True, text=True, check=True)
        self.assertTrue(json.loads(result.stdout)["expected_lengths_checked"])

    def test_mapping_review_only_suggests_unique_literal_matches(self):
        self.prepare("S\t1\tAAA\nW\tsample\t1\texact\t0\t3\t>1\n"
                     "W\tsample\t1\tshared\t0\t3\t>1\nW\tsample\t2\tshared\t0\t3\t>1\n"
                     "W\tsample\t1\tseq\t0\t3\t>1\nP\tunique\t1+\t*\nP\tcommon\t1+\t*\n"
                     "P\tsample#hap#near\t1+\t*\n")
        hal = self.root / "hal.tsv"
        pairs = [("sample", "exact"), ("sample", "shared"), ("sample", "unique"),
                 ("a", "common"), ("b", "common"), ("sample", "near"), ("sample", "Exact"), ("other", "seq")]
        mapping.write_tsv(hal, [dict(hal_genome=g, hal_sequence=s, length_bp=99) for g, s in pairs], mapping.HAL_INVENTORY_COLUMNS)
        output = self.root / "review"
        mapping.make_public_templates(argparse.Namespace(hal_inventory=str(hal), gfa_inventory=str(self.inventory), out_dir=str(output)))
        with (output / "paths.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(row["gfa_record_type"] == "" for row in rows))
        from pancgi_contract import PATHS
        self.assertEqual(list(rows[0]), PATHS)
        self.assertEqual(mapping.read_tsv(str(output / 'hal_sequences.tsv'), mapping.HAL_INVENTORY_COLUMNS),
                         mapping.read_tsv(str(hal), mapping.HAL_INVENTORY_COLUMNS))
        with (output / "mapping_review.tsv").open() as handle:
            review = {(row["hal_genome"], row["hal_sequence"]): row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(review[("sample", "exact")]["candidate_status"], "exact_literal_candidate")
        self.assertEqual(review[("sample", "exact")]["candidate_path_length_bp"], "3")
        self.assertEqual(review[("sample", "unique")]["candidate_status"], "exact_literal_candidate")
        self.assertEqual(review[("sample", "unique")]["candidate_path_length_bp"], "")
        self.assertEqual(review[("sample", "shared")]["candidate_status"], "multiple_exact_literal_matches")
        for key in pairs:
            if key not in {("sample", "exact"), ("sample", "unique")}:
                self.assertEqual(review[key]["candidate_gfa_path_id"], "")
            self.assertEqual(review[key]["sequence_identity_checked"], "false")


if __name__ == "__main__":
    unittest.main()
