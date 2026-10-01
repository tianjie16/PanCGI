import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json

import pancgi_mapping as mapping
import pancgi_contract as contract
from pancgi_contract import GENOMES, PATHS, SV, validate_tables


class MappingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root/'empty.bed').write_text('')
        self.table('empty.sv', SV, [])
        self.genomes = [dict(hal_genome=g, role=r, cpgi_bed='empty.bed', sv_tsv='empty.sv')
                        for g, r in [('REF', 'primary_reference'), ('unrelated-species', 'sample')]]
        self.paths = [dict(hal_genome=g, hal_sequence=c, gfa_record_type='P', gfa_path_name=f'opaque/{i}')
                      for i, (g, c) in enumerate([('REF', 'alpha'), ('unrelated-species', 'scaf')])]
        self.hal = [dict(hal_genome=p['hal_genome'], hal_sequence=p['hal_sequence'], length_bp=100,
                        top_segments=1, bottom_segments=1) for p in self.paths]
        self.gfa = [dict(zip(mapping.GFA_INVENTORY_COLUMNS, [f'path{i}', 'P', p['gfa_path_name'], '', '', p['gfa_path_name'],
                    0, 100, 100, 1, 2, 'zero_single_segment'])) for i, p in enumerate(self.paths)]

    def tearDown(self):
        self.tmp.cleanup()

    def table(self, name, columns, rows):
        mapping.write_tsv(self.root/name, rows, columns)

    def validate(self):
        self.table('genomes', GENOMES, self.genomes)
        self.table('paths', PATHS, self.paths)
        self.table('hal', mapping.HAL_INVENTORY_COLUMNS, self.hal)
        self.table('gfa', mapping.GFA_INVENTORY_COLUMNS, self.gfa)
        return validate_tables(*(str(self.root/n) for n in ['genomes', 'paths', 'hal', 'gfa']),
                               hal_only_exclusions=getattr(self, 'exclusions', None))

    def add_hal_only(self, genome='unrelated-species'):
        self.hal.append(dict(hal_genome=genome, hal_sequence='no-graph', length_bp=100,
                            top_segments=1, bottom_segments=1))
        self.exclusions = self.root/'exclude.tsv'
        self.table('exclude.tsv', contract.HAL_ONLY, [dict(hal_genome=genome, hal_sequence='no-graph')])

    def test_declared_unused_hal_only_passes_without_reducing_paths(self):
        self.add_hal_only()
        genomes, paths, hal = self.validate()
        self.assertEqual(len(genomes), 2)
        self.assertEqual(len(paths), 2)
        self.assertEqual(len(hal), 3)

    def test_hal_only_not_automatically_excluded(self):
        self.add_hal_only()
        self.exclusions = None
        with self.assertRaisesRegex(ValueError, 'Missing mappings'):
            self.validate()

    def test_hal_only_cgi_rejected(self):
        self.add_hal_only()
        (self.root/'sample.bed').write_text('no-graph\t0\t100\tcgi\t100\t10\t50\t20\t50\t0.8\n')
        self.genomes[1]['cpgi_bed'] = 'sample.bed'
        with self.assertRaisesRegex(ValueError, 'HAL-only exclusion has CGI'):
            self.validate()

    def test_hal_only_sv_assembly_rejected(self):
        self.add_hal_only()
        self.table('sample.sv', SV, [dict(zip(SV, ['alpha', 10, 10, 'x', 'INS', 'no-graph', 30, 40, 10]))])
        self.genomes[1]['sv_tsv'] = 'sample.sv'
        with self.assertRaisesRegex(ValueError, 'HAL-only exclusion has SV'):
            self.validate()

    def test_hal_only_sv_reference_rejected(self):
        self.add_hal_only('REF')
        self.table('sample.sv', SV, [dict(zip(SV, ['no-graph', 10, 20, 'x', 'DEL', 'scaf', 30, 30, -10]))])
        self.genomes[1]['sv_tsv'] = 'sample.sv'
        with self.assertRaisesRegex(ValueError, 'HAL-only exclusion has SV'):
            self.validate()

    def test_exclusion_of_mapped_sequence_rejected(self):
        self.add_hal_only()
        self.table('exclude.tsv', contract.HAL_ONLY, [dict(hal_genome='REF', hal_sequence='alpha')])
        with self.assertRaisesRegex(ValueError, 'without a mapping'):
            self.validate()

    def test_duplicate_or_unknown_exclusion_rejected(self):
        self.add_hal_only()
        row = dict(hal_genome='unrelated-species', hal_sequence='no-graph')
        for records in ([row, row], [dict(row, hal_sequence='unknown')], [dict(row, hal_genome='unselected')]):
            with self.subTest(records=records):
                self.table('exclude.tsv', contract.HAL_ONLY, records)
                with self.assertRaisesRegex(ValueError, 'unique selected HAL'):
                    self.validate()

    def test_whole_genome_cannot_be_excluded(self):
        self.paths.pop()
        self.exclusions = self.root/'exclude.tsv'
        self.table('exclude.tsv', contract.HAL_ONLY, [dict(hal_genome='unrelated-species', hal_sequence='scaf')])
        with self.assertRaisesRegex(ValueError, 'Every selected genome'):
            self.validate()

    def test_hal_only_exclusions_locked_and_reported(self):
        self.add_hal_only()
        self.validate()
        (self.root/'source.gfa').write_text('source')
        (self.root/'source.hal').write_text('source')
        args = SimpleNamespace(gfa=self.root/'source.gfa', hal=self.root/'source.hal',
            genomes=self.root/'genomes', paths=self.root/'paths', gfa_inventory=self.root/'gfa',
            hal_inventory=self.root/'hal', out_dir=self.root/'validated', hal_only_exclusions=self.exclusions)
        self.assertEqual(contract.prepare(args)['hal_only_excluded_n'], 1)
        lock = contract.verify_lock(self.root/'validated/inputs.lock.json')
        row = lock['hal_only_exclusions'][0]
        self.assertEqual(row['length_bp'], 100)
        self.assertEqual(row['cgi_rows'] + row['sv_assembly_rows'] + row['sv_reference_rows'], 0)
        self.assertIn(str(self.exclusions), lock['files'])
        self.assertIn('hal_only_exclusions.tsv', lock['internal_files'])
        (self.root/'validated/hal_only_exclusions.tsv').write_text('changed')
        with self.assertRaisesRegex(RuntimeError, 'identity mapping changed'):
            contract.verify_lock(self.root/'validated/inputs.lock.json')

    def test_explicit_names_not_inferred(self):
        genomes, paths, hal = self.validate()
        self.assertEqual(genomes[1]['hal_genome'], 'unrelated-species')

    def test_no_missing_paths_even_without_cgi(self):
        self.paths.pop()
        with self.assertRaises(ValueError):
            self.validate()

    def test_reserved_genotype_headers(self):
        self.genomes[1]['hal_genome'] = 'locus_id'
        with self.assertRaisesRegex(ValueError, 'reserved output column'):
            self.validate()

    def test_no_shared_path_assignment(self):
        self.paths[1]['gfa_path_name'] = self.paths[0]['gfa_path_name']
        with self.assertRaises(ValueError):
            self.validate()

    def test_no_length_adaptation(self):
        self.gfa[0]['end0'] = 99
        with self.assertRaises(ValueError):
            self.validate()

    def test_p_length_deferred_to_actual_unfold(self):
        self.gfa[0]['end0'] = ''
        self.gfa[0]['path_length_bp'] = ''
        self.validate()

    def test_contract_no_whole_graph_or_hal_read(self):
        self.validate()
        graph = self.root/'source.gfa'
        hal = self.root/'source.hal'
        graph.write_bytes(b'input integrity supplied by user')
        hal.write_bytes(b'not opened by metadata contract')
        original = contract.digest

        def bounded_digest(path):
            if Path(path).resolve() in {graph.resolve(), hal.resolve()}:
                raise AssertionError('Whole graph or HAL hashing is not an input validation step')
            return original(path)

        args = SimpleNamespace(gfa=str(graph), hal=str(hal), genomes=str(self.root/'genomes'),
            paths=str(self.root/'paths'), gfa_inventory=str(self.root/'gfa'),
            hal_inventory=str(self.root/'hal'), out_dir=str(self.root/'validated'))
        with patch.object(contract, 'digest', bounded_digest), patch('subprocess.Popen', side_effect=AssertionError('No sequence export')):
            contract.prepare(args)
        lock = contract.verify_lock(self.root/'validated/inputs.lock.json')
        self.assertFalse(lock['sequence_identity_checked'])
        self.assertEqual(lock['selected_haplotype_n'], 1)
        lengths = json.loads((self.root/'validated/path_lengths.json').read_text())
        self.assertEqual(lengths, {'path0': 100, 'path1': 100})
        for item in lock['logical_inputs']:
            if item['role'] in ('gfa', 'hal'):
                self.assertIsNone(item['sha256'])
                self.assertEqual(item['integrity_method'], 'file_metadata')
        (self.root/'validated/path_lengths.json').write_text('{}\n')
        with self.assertRaisesRegex(RuntimeError, 'identity mapping changed'):
            contract.verify_lock(self.root/'validated/inputs.lock.json')

    def test_unknown_sv_does_not_disappear(self):
        self.table('empty.sv', SV, [dict(zip(SV, ['alpha', 10, 20, 'x', 'DUP', 'alpha', 30, 40, 10]))])
        with self.assertRaises(ValueError):
            self.validate()

    def test_reference_sv_cannot_be_silently_ignored(self):
        self.table('empty.sv', SV, [dict(zip(SV, ['alpha', 10, 10, 'x', 'INS', 'alpha', 30, 40, 10]))])
        with self.assertRaisesRegex(ValueError, 'Reference-role SV'):
            self.validate()

    def test_no_name_based_template_linking(self):
        from argparse import Namespace
        self.validate()
        mapping.make_public_templates(Namespace(hal_inventory=str(self.root/'hal'),
            gfa_inventory=str(self.root/'gfa'), out_dir=str(self.root/'templates')))
        with (self.root/'templates/paths.tsv').open() as handle:
            self.assertTrue(all(r['gfa_record_type'] == '' for r in csv.DictReader(handle, delimiter='\t')))

    def test_w_literals_and_mixed_record_types(self):
        self.gfa[1].update(record_type='W', raw_path_name='', w_sample_id='any#name',
                           w_haplotype_index='maternal-custom', sequence_id='arbitrary/scaffold')
        self.paths[1].update(gfa_record_type='W', gfa_path_name='', gfa_sample='any#name',
                             gfa_haplotype='maternal-custom', gfa_sequence='arbitrary/scaffold')
        self.assertEqual(self.validate()[1][1]['gfa_path_id'], 'path1')
        self.paths[1]['gfa_haplotype'] = 'different'
        with self.assertRaisesRegex(ValueError, 'Unknown literal'):
            self.validate()

    def test_p_does_not_accept_w_fields(self):
        self.paths[0]['gfa_sample'] = 'not-used'
        with self.assertRaisesRegex(ValueError, 'W fields must be empty'):
            self.validate()


if __name__ == '__main__':
    unittest.main()
