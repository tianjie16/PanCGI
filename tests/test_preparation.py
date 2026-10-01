import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pancgi_mapping as mapping
import pancgi_preparation as preparation
from pancgi_contract import digest, file_state


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.gfa, self.hal = self.root / 'graph.gfa', self.root / 'alignment.hal'
        self.gfa.write_text('S\t1\tAAA\nP\tp\t1+\t*\n')
        self.hal.write_bytes(b'only file metadata is needed here')
        self.prepared = self.root / 'prepared'
        self.prepared.mkdir()
        mapping.inventory_gfa(str(self.gfa), str(self.prepared / 'gfa_paths.tsv'))
        for name in ('hal_sequences.tsv', 'hal_genomes.tsv'):
            (self.prepared / name).write_text('inventory\n')
        report = dict(schema_version=1, inventory_schema_version=2, status='complete',
            sources=dict(gfa=str(self.gfa), hal=str(self.hal)),
            source_states={str(p): file_state(p) for p in (self.gfa, self.hal)},
            files={name: digest(self.prepared / name) for name in preparation.FILES})
        (self.prepared / 'preparation.json').write_text(json.dumps(report))
        self.args = SimpleNamespace(prepared=self.prepared, gfa=self.gfa, hal=self.hal, out_dir=self.root / 'run')

    def tearDown(self):
        self.temp.cleanup()

    def test_reuses_inventories_without_source_reads(self):
        original = preparation.digest
        def guard(path):
            self.assertNotIn(Path(path), (self.gfa, self.hal))
            return original(path)
        with patch.object(preparation, 'digest', guard), patch.object(mapping, 'inventory_gfa', side_effect=AssertionError('no rescan')):
            result = preparation.validate_and_copy(self.args)
        self.assertEqual(result['gfa_source_scans'], 0)
        for name in preparation.FILES:
            self.assertEqual((self.args.out_dir / name).read_bytes(), (self.prepared / name).read_bytes())

    def test_changed_source_rejects(self):
        self.hal.write_bytes(b'different alignment')
        with self.assertRaisesRegex(ValueError, 'changed after preparation'):
            preparation.validate_and_copy(self.args)

    def test_changed_inventory_rejects(self):
        (self.prepared / 'gfa_paths.tsv').write_text('corrupt')
        with self.assertRaisesRegex(ValueError, 'inventory changed'):
            preparation.validate_and_copy(self.args)

    def test_wrong_input_rejects(self):
        self.args.gfa = self.hal
        with self.assertRaisesRegex(ValueError, 'does not belong'):
            preparation.validate_and_copy(self.args)

    def set_old_device(self):
        sidecar = self.prepared / 'gfa_paths.tsv.json'
        metadata = json.loads(sidecar.read_text())
        metadata['source_fingerprint']['device'] += 1000
        sidecar.write_text(json.dumps(metadata))
        receipt = self.prepared / 'preparation.json'
        report = json.loads(receipt.read_text())
        for state in report['source_states'].values():
            state[0] += 1000
        report['files']['gfa_paths.tsv.json'] = digest(sidecar)
        receipt.write_text(json.dumps(report))
        return report

    def test_cross_node_device_only_change_preserves_originals(self):
        report = self.set_old_device()
        before = {p.name: p.read_bytes() for p in self.prepared.iterdir() if p.is_file()}
        with patch.object(mapping, 'inventory_gfa', side_effect=AssertionError('no source rescan')), \
             patch.object(mapping, 'inventory_hal', side_effect=AssertionError('no HAL rescan')):
            result = preparation.validate_and_copy(self.args)
        self.assertEqual(result['original_source_states'], report['source_states'])
        self.assertEqual(result['run_source_states'], {str(p): file_state(p) for p in (self.gfa, self.hal)})
        self.assertTrue(result['within_run_full_state_checked'])
        self.assertEqual(result['preparation_sha256'], digest(self.prepared / 'preparation.json'))
        self.assertEqual(json.loads((self.args.out_dir / 'preparation_reuse.json').read_text()), result)
        for name, data in before.items():
            self.assertEqual((self.prepared / name).read_bytes(), data)
            if name in (*preparation.FILES, 'preparation.json'):
                self.assertEqual((self.args.out_dir / name).read_bytes(), data)

    def test_portable_reuse_rejects_each_other_state_field(self):
        report = self.set_old_device()
        receipt = self.prepared / 'preparation.json'
        for path in (self.gfa, self.hal):
            for index in range(1, 5):
                with self.subTest(path=path.name, index=index):
                    altered = json.loads(json.dumps(report))
                    altered['source_states'][str(path)][index] += 1
                    receipt.write_text(json.dumps(altered))
                    with self.assertRaisesRegex(ValueError, 'changed after preparation'):
                        preparation.validate_and_copy(self.args)

    def test_portable_reuse_rejects_malformed_state(self):
        report = self.set_old_device()
        for state in (None, [], [1, 2, 3, 4], [True, 2, 3, 4, 5], ['1', 2, 3, 4, 5]):
            with self.subTest(state=state):
                altered = json.loads(json.dumps(report))
                altered['source_states'][str(self.hal)] = state
                (self.prepared / 'preparation.json').write_text(json.dumps(altered))
                with self.assertRaisesRegex(ValueError, 'changed after preparation'):
                    preparation.validate_and_copy(self.args)

    def test_portable_reuse_rejects_corrupt_inventory(self):
        self.set_old_device()
        (self.prepared / 'hal_sequences.tsv').write_text('changed inventory')
        with self.assertRaisesRegex(ValueError, 'inventory changed'):
            preparation.validate_and_copy(self.args)

    def test_every_source_state_change_during_copy_is_rejected(self):
        self.set_old_device()
        for path in (self.gfa, self.hal):
            for index in range(5):
                with self.subTest(path=path.name, index=index):
                    self.args.out_dir = self.root / f'changed_{path.name}_{index}'
                    counts = {}
                    def changing(p):
                        p = Path(p)
                        counts[p] = counts.get(p, 0) + 1
                        state = file_state(p)
                        if p == path and counts[p] > 1:
                            state[index] += 1
                        return state
                    with patch.object(preparation, 'file_state', changing), \
                         self.assertRaisesRegex(RuntimeError, 'changed during preparation reuse'):
                        preparation.validate_and_copy(self.args)
                    self.assertFalse((self.args.out_dir / 'preparation_reuse.json').exists())

    def test_receipt_change_during_copy_is_rejected(self):
        original = preparation.shutil.copyfile
        def corrupt(source, target):
            result = original(source, target)
            if Path(source).name == 'preparation.json':
                Path(target).write_text('{}')
            return result
        with patch.object(preparation.shutil, 'copyfile', corrupt), \
             self.assertRaisesRegex(RuntimeError, 'changed while copying'):
            preparation.validate_and_copy(self.args)
        self.assertFalse((self.args.out_dir / 'preparation_reuse.json').exists())
