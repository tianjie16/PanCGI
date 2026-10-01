from __future__ import annotations

import gzip
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pancgi_contract as contract
import pancgi_hal as hal
import pancgi_hal_runtime as runtime


class FastaHeaderTests(unittest.TestCase):
    def test_exact_names_and_overlapping_wrapped_intervals(self):
        records = [("chr 1#part:2", 1, 7, "first"), ("chr 1#part:2", 3, 5, "overlap"),
                   ("chr", 0, 2, "short_name")]
        sequences, observed = runtime.extract_intervals_from_fasta_stream(
            io.StringIO(">chr 1#part:2\r\nacgt\r\nacgt\r\n>chr\nTT\n"), records)
        self.assertEqual(sequences, {"first": "CGTACG", "overlap": "TA", "short_name": "TT"})
        self.assertEqual(observed, {"first": 6, "overlap": 2, "short_name": 2})

    def test_does_not_guess_name_from_first_token(self):
        sequences, observed = runtime.extract_intervals_from_fasta_stream(
            io.StringIO(">chr1 description\nACGT\n"), [("chr1", 0, 4, "id")])
        self.assertEqual(sequences, {"id": ""})
        self.assertEqual(observed, {"id": 0})

    def test_malformed_headers_rejected(self):
        for header in (">", "> ", "> chr1", ">chr1 ", ">chr\t1", ">chr\x001",
                       ">chr\x7f1", " >chr1", "\t>chr1"):
            with self.subTest(header=header), self.assertRaisesRegex(ValueError, "malformed header"):
                runtime.extract_intervals_from_fasta_stream(
                    io.StringIO(header + "\nACGT\n"), [("chr1", 0, 4, "id")])

    def test_duplicate_headers_rejected_even_for_unrequested_contigs(self):
        for name in ("chr1", "other", "chr1 exact name"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "duplicate header"):
                runtime.extract_intervals_from_fasta_stream(
                    io.StringIO(f">{name}\nACGT\n>{name}\nACGT\n"), [("chr1", 0, 4, "id")])

    def test_sequence_before_header_rejected(self):
        with self.assertRaisesRegex(ValueError, "before a header"):
            runtime.extract_intervals_from_fasta_stream(io.StringIO("ACGT\n"), [])


class HalOverheadTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="pancgi_hal_overhead_")
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        environment = patch.dict(os.environ, {"PANCGI_HAL_RUNTIME": "docker"})
        environment.start()
        self.addCleanup(environment.stop)
        self.bed = self.root / "input.bed"
        self.bed.write_text("G#0#C\t1\t5\tCpG: 1\t4\t1\t2\t50\t50\t1\n", encoding="ascii")
        self.hal_path = self.root / "input.hal"
        self.hal_path.write_bytes(b"HAL fixture")
        self.fasta = self.root / "out.fa.gz"
        self.row = dict(label="G", kind="assembly", graph_sample="G", graph_hap="0",
                        hal_genome="source", bed=str(self.bed), cpgi_fa=str(self.fasta))
        self.contigs = {("G", "C"): "hal_sequence_A"}
        self.fid = "G#0#C:1-5"
        self.psl_text = "\t".join([self.fid, "4", "0", "0", "0", "0", "0", "0", "0", "++",
                                   "hal_sequence_A", "6", "1", "5", "chr1", "100", "10", "14",
                                   "1", "4,", "1,", "10,"]) + "\n"
        self.psl = self.root / "psl" / "G.to_primary.psl"
        self.metadata = self.root / "metadata" / "G.json"
        self.log = self.root / "logs" / "G.hal2fasta.log"

    def extract(self):
        return runtime.extract_cpgi_fasta_from_hal(
            self.row, str(self.hal_path), "docker", "image", "hal2fasta", str(self.log), self.contigs)

    def project(self, overwrite=False):
        return hal.run_one(self.row, self.hal_path, "target", self.root, "docker", "image",
                           "image-id", "halLiftover", overwrite, self.contigs)

    def signature(self, image="image", image_id="image-id"):
        return hal.signature(self.row, "bed-digest", self.hal_path, "target", image, image_id, "halLiftover")

    def generated_psl(self, command, **kwargs):
        Path(command[-1]).write_text(self.psl_text * 2, encoding="ascii")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    def test_empty_plain_and_gzip_inputs_skip_all_hal_execution(self):
        for suffix in (".bed", ".bed.gz"):
            for output_suffix in (".fa", ".fa.gz"):
                with self.subTest(bed=suffix, fasta=output_suffix):
                    bed = self.root / ("empty" + suffix)
                    payload = b"" if suffix == ".bed" else b"# no intervals\n\n"
                    bed.write_bytes(gzip.compress(payload) if suffix.endswith(".gz") else payload)
                    output = self.root / ("empty" + output_suffix)
                    output.write_bytes(b"stale output")
                    self.row.update(bed=str(bed), cpgi_fa=str(output))
                    with patch.object(runtime, "container_command") as command, \
                            patch.object(runtime.subprocess, "Popen") as process:
                        result = self.extract()
                    command.assert_not_called()
                    process.assert_not_called()
                    self.assertEqual(result, dict(label="G", records=0, missing=[], output=str(output)))
                    data = output.read_bytes()
                    self.assertEqual(gzip.decompress(data) if output_suffix.endswith(".gz") else data, b"")
                    self.assertIn("no HAL FASTA extraction required", self.log.read_text())

    def test_nonempty_extraction_preserves_command_and_does_not_reread_output(self):
        process = Mock(stdout=io.StringIO(">hal_sequence_A\nacgtac\n"))
        process.wait.return_value = 0
        original_open = Path.open

        def checked_open(path, mode="r", *args, **kwargs):
            if path == self.fasta and "r" in mode:
                raise AssertionError("Generated FASTA must not be reread")
            return original_open(path, mode, *args, **kwargs)

        with patch.object(runtime.subprocess, "Popen", return_value=process) as popen, \
                patch.object(Path, "open", new=checked_open):
            result = self.extract()
        self.assertEqual(popen.call_args.args[0][-3:], ["--upper", str(self.hal_path), "source"])
        self.assertEqual(result, dict(label="G", records=1, missing=[], output=str(self.fasta)))
        self.assertEqual(gzip.decompress(self.fasta.read_bytes()), b">G#0#C:1-5\nCGTA\n")
        process.kill.assert_not_called()
        process.wait.assert_called_once_with()

    def test_incomplete_intervals_do_not_publish_fasta(self):
        for stream, length in ((">hal_sequence_A\nAC\n", 1), (">other\nACGTAC\n", 0),
                               (">hal_sequence_A description\nACGTAC\n", 0)):
            with self.subTest(stream=stream):
                self.fasta.write_bytes(b"existing output")
                process = Mock(stdout=io.StringIO(stream))
                process.wait.return_value = 0
                with patch.object(runtime.subprocess, "Popen", return_value=process):
                    result = self.extract()
                self.assertEqual(result["missing"], [(self.fid, 4, length)])
                self.assertEqual(self.fasta.read_bytes(), b"existing output")

    def test_nonzero_hal_exit_does_not_publish_complete_intervals(self):
        process = Mock(stdout=io.StringIO(">hal_sequence_A\nACGTAC\n"))
        process.wait.return_value = 7
        with patch.object(runtime.subprocess, "Popen", return_value=process), \
                self.assertRaisesRegex(RuntimeError, "exit code 7"):
            self.extract()
        self.assertFalse(self.fasta.exists())

    def test_invalid_fasta_kills_and_reaps_hal_without_publishing(self):
        for text in (">\nACGTAC\n", ">hal_sequence_A\nACGTAC\n>hal_sequence_A\nACGTAC\n"):
            with self.subTest(text=text):
                process = Mock(stdout=io.StringIO(text))
                with patch.object(runtime.subprocess, "Popen", return_value=process), \
                        self.assertRaises(ValueError):
                    self.extract()
                process.kill.assert_called_once_with()
                process.wait.assert_called_once_with()
                self.assertFalse(self.fasta.exists())

    def test_invalid_bed_is_not_treated_as_empty(self):
        self.bed.write_text("invalid BED\n", encoding="ascii")
        with patch.object(runtime.subprocess, "Popen") as process, self.assertRaises(ValueError):
            self.extract()
        process.assert_not_called()
        self.assertFalse(self.fasta.exists())

    def test_generated_and_reused_psl_are_parsed_once_each(self):
        with patch.object(hal.subprocess, "run", side_effect=self.generated_psl) as run, \
                patch.object(hal, "parse_psl", wraps=hal.parse_psl) as parse:
            generated = self.project()
        parse.assert_called_once()
        self.assertIn(".partial.", parse.call_args.args[0].name)
        self.assertEqual(run.call_args.args[0][-8:-1],
                         ["--bedType", "4", "--outPSLWithName", str(self.hal_path), "source",
                          str(self.root / "bed4" / "G.bed4"), "target"])
        with patch.object(hal.subprocess, "run") as run, \
                patch.object(hal, "parse_psl", wraps=hal.parse_psl) as parse:
            reused = self.project()
        run.assert_not_called()
        parse.assert_called_once_with(self.psl, {self.fid})
        self.assertEqual(generated["status"], "generated")
        self.assertEqual(reused["status"], "validated_existing")
        self.assertEqual({k: v for k, v in generated.items() if k != "status"},
                         {k: v for k, v in reused.items() if k != "status"})
        self.assertEqual((reused["psl_records"], reused["mapped_ids"], reused["multimap_ids"]), (2, 1, 1))
        self.assertEqual(reused["multimapped"], [(self.fid, 2)])

    def test_empty_projection_skips_hal_on_generation_reuse_and_overwrite(self):
        self.bed.write_bytes(b"")
        with patch.object(runtime, "container_command") as command, patch.object(hal.subprocess, "run") as run:
            generated = self.project()
            reused = self.project()
            self.metadata.write_text("{}", encoding="ascii")
            with self.assertRaisesRegex(RuntimeError, "does not match current inputs"):
                self.project()
            self.psl.write_text("stale PSL", encoding="ascii")
            overwritten = self.project(overwrite=True)
        command.assert_not_called()
        run.assert_not_called()
        self.assertEqual(self.psl.read_bytes(), b"")
        self.assertEqual((generated["status"], reused["status"], overwritten["status"]),
                         ("generated", "validated_existing", "generated"))
        for result in (generated, reused, overwritten):
            for key in ("input_intervals", "psl_records", "mapped_ids", "unmapped_ids", "multimap_ids"):
                self.assertEqual(result[key], 0)

    def test_psl_validation_rejects_malformed_and_unknown_rows_before_publication(self):
        for text in ("bad\tPSL\n", self.psl_text.replace(self.fid, "unknown")):
            with self.subTest(text=text):
                self.psl_text = text
                with patch.object(hal.subprocess, "run", side_effect=self.generated_psl), \
                        self.assertRaises(ValueError):
                    self.project()
                self.assertFalse(self.psl.exists())
                self.assertFalse(self.metadata.exists())

    def test_reused_psl_still_requires_structural_validation(self):
        with patch.object(hal.subprocess, "run", side_effect=self.generated_psl):
            self.project()
        for text in ("bad\tPSL\n", self.psl_text.replace(self.fid, "unknown")):
            with self.subTest(text=text):
                self.psl.write_text(text, encoding="ascii")
                with patch.object(hal.subprocess, "run") as run, self.assertRaises(ValueError):
                    self.project()
                run.assert_not_called()

    def test_liftover_failure_or_missing_output_does_not_publish(self):
        for code in (1, 0):
            with self.subTest(code=code):
                result = SimpleNamespace(returncode=code, stdout="", stderr="failed")
                with patch.object(hal.subprocess, "run", return_value=result), self.assertRaises(RuntimeError):
                    self.project()
                self.assertFalse(self.psl.exists())
                self.assertFalse(self.metadata.exists())

    def test_signature_uses_complete_contract_file_state_without_hal_hash(self):
        with patch.object(hal, "sha256_file", side_effect=AssertionError("No HAL hash")), \
                patch.object(contract, "digest", side_effect=AssertionError("No HAL hash")):
            signature = self.signature()
        self.assertEqual(signature["hal_file_state"], contract.file_state(self.hal_path))
        self.assertEqual(len(signature["hal_file_state"]), 5)
        self.assertIsNone(signature["hal_liftover_identity"])

    def test_device_inode_and_ctime_changes_invalidate_reuse(self):
        with patch.object(hal.subprocess, "run", side_effect=self.generated_psl):
            self.project()
        state = contract.file_state(self.hal_path)
        for index in (0, 1, 4):
            changed = state.copy()
            changed[index] += 1
            with self.subTest(field=index), patch.object(contract, "file_state", return_value=changed), \
                    patch.object(hal.subprocess, "run") as run, \
                    self.assertRaisesRegex(RuntimeError, "does not match current inputs"):
                self.project()
            run.assert_not_called()

    def test_replacement_hal_with_same_size_and_mtime_changes_signature(self):
        original = self.signature()
        stat = self.hal_path.stat()
        replacement = self.root / "replacement.hal"
        replacement.write_bytes(self.hal_path.read_bytes())
        os.utime(replacement, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        os.replace(replacement, self.hal_path)
        updated = self.signature()
        self.assertEqual(original["hal_file_state"][2:4], updated["hal_file_state"][2:4])
        self.assertNotEqual(original, updated)

    def test_incomplete_metadata_requires_explicit_overwrite(self):
        with patch.object(hal.subprocess, "run", side_effect=self.generated_psl):
            self.project()
        metadata = json.loads(self.metadata.read_text())
        metadata.pop("hal_file_state")
        metadata.pop("hal_liftover_identity")
        metadata.update(hal_size=self.hal_path.stat().st_size, hal_mtime_ns=self.hal_path.stat().st_mtime_ns)
        self.metadata.write_text(json.dumps(metadata), encoding="ascii")
        with patch.object(hal.subprocess, "run") as run, \
                self.assertRaisesRegex(RuntimeError, "does not match current inputs"):
            self.project()
        run.assert_not_called()
        with patch.object(hal.subprocess, "run", side_effect=self.generated_psl) as run:
            self.assertEqual(self.project(overwrite=True)["status"], "generated")
        run.assert_called_once()

    def test_native_identity_resolves_symlink_and_tracks_executable_state(self):
        executable = self.root / "halLiftover-real"
        executable.write_bytes(b"executable fixture")
        alias = self.root / "halLiftover"
        alias.symlink_to(executable)
        with patch.dict(os.environ, {"PANCGI_HAL_RUNTIME": "native"}), \
                patch.object(runtime.shutil, "which", return_value=str(alias)) as which:
            original = self.signature(image="", image_id="native")
            self.assertEqual(original["hal_liftover_identity"],
                             {"path": str(executable), "file_state": contract.file_state(executable)})
            executable.write_bytes(b"changed executable fixture")
            self.assertNotEqual(original, self.signature(image="", image_id="native"))
        self.assertEqual(which.call_count, 2)

    def test_native_identity_rejects_missing_executable_and_docker_image(self):
        with patch.dict(os.environ, {"PANCGI_HAL_RUNTIME": "native"}), \
                patch.object(runtime.shutil, "which", return_value=None):
            with self.assertRaises(FileNotFoundError):
                self.signature(image="", image_id="native")
            with self.assertRaises(ValueError):
                self.signature()


if __name__ == "__main__":
    unittest.main()
