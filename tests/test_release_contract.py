import csv
import gzip
import json
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import contextlib
import io
from pathlib import Path

import cpgi_sv_annot as sv
import cpgi_nr_prod as prod
from pancgi_contract import GENOMES, PATHS, SV, read_exact, read_bed, storage_id, validate_sv, digest
from pancgi_inputs import prepare_inputs
from pancgi_genotyping.caller import evaluate, propagate
import pancgi_mapping as mapping


class ReleaseContractTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def test_exact_header_and_width(self):
        for text in ('x\ty\n1\t2\n', 'hal_genome\trole\tcpgi_bed\tsv_tsv\textra\n',
                     '\t'.join(GENOMES) + '\na\tsample\tb\ts\textra\n'):
            with self.assertRaises(ValueError):
                read_exact(self.write('bad.tsv', text), GENOMES)

    def test_empty_inputs_explicit(self):
        path = self.write('empty.tsv', '\t'.join(SV) + '\n')
        self.assertEqual(validate_sv(path, {'q': 100}, {'r': 100}), [])
        self.assertEqual(list(read_bed(self.write('empty.bed', ''), {'q': 100})), [])
        with self.assertRaises(FileNotFoundError):
            read_exact(self.root / 'absent', SV, True)

    def test_sv_coordinate_contract(self):
        header = '\t'.join(SV) + '\n'
        for text in ('r\t10\t10\ti\tINS\tq\t20\t30\t10\n', 'r\t10\t20\td\tDEL\tq\t30\t30\t-10\n'):
            self.assertEqual(len(validate_sv(self.write('sv', header + text), {'q': 100}, {'r': 100})), 1)
        for text in ('r\t0\t0\ti\tINS\tq\t20\t30\t10\n', 'r\t10\t20\ti\tINS\tq\t20\t30\t10\n',
                     'r\t10\t20\td\tDEL\tq\t30\t31\t-10\n', 'r\t10\t20\tx\tINV\tq\t30\t40\t10\n',
                     'r\t10\t20\tx\tDUP\tq\t30\t40\t10\n'):
            with self.assertRaises(ValueError):
                validate_sv(self.write('sv', header + text), {'q': 100}, {'r': 100})

    def index(self, records):
        fields = ['ID', 'contig', 'start', 'end', 'SVTYPE', '#CHROM', 'POS', 'VCF_END', 'SVLEN']
        path = self.root / 'internal.tsv'
        with path.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(fields)
            writer.writerows(records)
        return sv.build_sv_index(str(path))

    def test_del_junction_strict_and_length(self):
        index = self.index([['d', 'q', 500, 500, 'DEL', 'r', 100, 200, -100]])
        for start, end, expected in [(400, 600, 1), (500, 600, 0), (400, 500, 0)]:
            result = sv.annotate_feature_overlaps('q', start, end, index)
            self.assertEqual(result['sv_overlap_n'], expected)
            if expected:
                detail = json.loads(result['sv_overlap_detail_json'])[0]
                self.assertEqual((detail['span_bp'], detail['svlen_abs'], detail['overlap_bp']), (0, 100, 0))
                self.assertIsNone(detail['overlap_pct_sv'])
                self.assertEqual(result['sv_overlap_pct_sv'], 'NA')

    def test_arbitrary_interval_retrieval(self):
        index = self.index([['long', 'q', 0, 1000, 'INS', 'r', 10, 10, 1000], ['short', 'q', 100, 150, 'INS', 'r', 20, 20, 50]])
        result = sv.annotate_feature_overlaps('q', 500, 600, index)
        self.assertEqual(result['sv_ins_longest_id'], 'long')
        self.assertEqual(prod.query_ref_intervals_overlap({'r': [(1, 1000, 'long'), (101, 150, 'short')]}, 'r', 501, 600), [(1, 1000, 'long')])

    def test_all_ins_reference_longest_matching(self):
        index = self.index([['long', 'q', 100, 200, 'INS', 'r', 900, 900, 100],
                            ['mid', 'q', 300, 380, 'INS', 'r', 120, 120, 80],
                            ['small', 'q', 400, 450, 'INS', 'r', 220, 220, 50]])
        feat = sv.annotate_feature_overlaps('q', 50, 500, index)
        ref = {'r': [(100, 150, 'reference-A'), (200, 250, 'reference-B')]}
        result = prod.choose_sv_site_ref_candidate(feat, ref, {})
        self.assertEqual(feat['sv_ins_longest_id'], 'long')
        self.assertEqual(result['reference_ins_id'], 'mid')
        self.assertEqual(result['best_ref_fid'], 'reference-A')
        self.assertEqual(result['candidate_n'], 2)

    def test_anchor_base_inclusive(self):
        feat = sv.annotate_feature_overlaps('q', 0, 50, self.index([['i', 'q', 10, 40, 'INS', 'r', 100, 100, 30]]))
        self.assertIsNotNone(prod.choose_sv_site_ref_candidate(feat, {'r': [(50, 100, 'a')]}, {}))
        self.assertIsNone(prod.choose_sv_site_ref_candidate(feat, {'r': [(101, 200, 'b')]}, {}))

    def test_sv_length_required_at_both_readers(self):
        header = '\t'.join(SV) + '\n'
        for kind, pos, end, start, stop, invalid in [('INS', 10, 10, 20, 30, ['', '.', '0', '-10', '10.0', '10,20']),
                                                    ('DEL', 10, 20, 30, 30, ['', '.', '0', '10', '-10.0', '-10,-20'])]:
            for length in invalid:
                with self.subTest(kind=kind, length=length):
                    row = f'r\t{pos}\t{end}\ti\t{kind}\tq\t{start}\t{stop}\t{length}\n'
                    with self.assertRaises(ValueError):
                        validate_sv(self.write('sv', header + row), {'q':100}, {'r':100})
                    with self.assertRaises(ValueError):
                        self.index([['i', 'q', start, stop, kind, 'r', pos, end, length]])
        with self.assertRaises(ValueError):
            validate_sv(self.write('old.tsv', '\t'.join(SV[:-1])+'\nr\t10\t10\ti\tINS\tq\t20\t30\n'), {'q':100}, {'r':100})
        with self.assertRaises(ValueError):
            sv.build_sv_index(str(self.write('old_internal.tsv', 'ID\tcontig\tstart\tend\tSVTYPE\t#CHROM\tPOS\tVCF_END\ni\tq\t20\t30\tINS\tr\t10\t10\n')))

    def test_supplied_length_not_span_controls_ins_selection(self):
        index = self.index([['original_longer', 'q', 100, 605, 'INS', 'r', 120, 120, 543],
                            ['span_longer', 'q', 100, 610, 'INS', 'r', 220, 220, 510]])
        feat = sv.annotate_feature_overlaps('q', 200, 300, index)
        details = {r['id']:r for r in json.loads(feat['sv_overlap_detail_json'])}
        self.assertEqual((details['original_longer']['span_bp'], details['original_longer']['svlen']), (505, 543))
        self.assertAlmostEqual(details['original_longer']['overlap_pct_sv'], 100/505)
        self.assertEqual(feat['sv_ins_longest_id'], 'original_longer')
        result = prod.choose_sv_site_ref_candidate(feat, {'r': [(100,150,'ref-A'), (200,250,'ref-B')]}, {})
        self.assertEqual((result['reference_ins_id'], result['best_ref_fid']), ('original_longer', 'ref-A'))
        self.assertEqual(prod.sv_best_interval_for_type(feat, 'INS')['id'], 'original_longer')

    def test_supplied_del_length_not_reference_span(self):
        index = self.index([['d', 'q', 500, 500, 'DEL', 'r', 100, 200, -95]])
        detail = json.loads(sv.annotate_feature_overlaps('q', 400, 600, index)['sv_overlap_detail_json'])[0]
        self.assertEqual((detail['span_bp'], detail['svlen'], detail['svlen_abs']), (0, -95, 95))

    def test_ins_contact_does_not_create_overlap(self):
        index = self.index([['i', 'q', 100, 200, 'INS', 'r', 10, 10, 100]])
        for start, end, expected in [(0,100,0), (200,300,0), (99,101,1), (199,201,1)]:
            self.assertEqual(sv.annotate_feature_overlaps('q', start, end, index)['sv_overlap_n'], expected)

    def test_preparation_preserves_supplied_length_and_coordinates(self):
        bed = self.write('empty.bed', '')
        ref_sv = self.write('ref.sv', '\t'.join(SV)+'\n')
        sample_sv = self.write('sample.sv', '\t'.join(SV)+'\nr\t120\t120\ti\tINS\tq\t100\t605\t543\nr\t100\t200\td\tDEL\tq\t800\t800\t-95\n')
        genomes = [dict.fromkeys(mapping.GENOME_COLUMNS, '') for _ in range(2)]
        for row, gid, role, source in zip(genomes, ['R','S'], ['primary_reference','sample'], [ref_sv,sample_sv]):
            row.update(genome_id=gid, role=role, hal_genome=gid, cpgi_bed=str(bed), sv_file=str(source))
        contigs = [dict.fromkeys(mapping.CONTIG_COLUMNS, '') for _ in range(2)]
        for row, gid, sequence, cid in zip(contigs, ['R','S'], ['r','q'], ['Rseq','Sseq']):
            row.update(genome_id=gid, hal_sequence=sequence, contig_id=cid)
        gf, cf = self.root/'genomes.tsv', self.root/'contigs.tsv'
        mapping.write_tsv(gf, genomes, mapping.GENOME_COLUMNS)
        mapping.write_tsv(cf, contigs, mapping.CONTIG_COLUMNS)
        hashes = {str(p):digest(p) for p in [bed,ref_sv,sample_sv]}
        self.write('inputs.lock.json', json.dumps(dict(status='validated', files=hashes, internal_files={p.name:digest(p) for p in [gf,cf]})))
        self.write('pathbed_inventory.tsv', '')
        output = self.root/'prepared'
        prepare_inputs(SimpleNamespace(genomes=gf, contigs=cf, out_dir=output, pathbed_dir=self.root))
        with (output/'S.sv.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([(r['start'],r['end'],r['POS'],r['VCF_END'],r['SVLEN']) for r in rows],
                         [('100','605','120','120','543'), ('800','800','100','200','-95')])
        self.assertEqual({str(p):digest(p) for p in [bed,ref_sv,sample_sv]}, hashes)

    def test_name_collision_resistance(self):
        values = ['a__b', 'a/b', 'a#b', 'a%b', 'a:b', 'a_b']
        self.assertEqual(len({storage_id('G', value) for value in values}), len(values))

    def test_genotype_propagation(self):
        for call in ('0', 'NA'):
            locus, alleles = propagate(True, {'a'}, ['a', 'b'], call)
            self.assertEqual((locus, alleles), ('1', {'a': '1', 'b': '0'}))
        self.assertEqual(propagate(False, set(), ['a', 'b'], '0'), ('0', {'a': '0', 'b': '0'}))
        self.assertEqual(propagate(False, set(), ['a', 'b'], 'NA'), ('NA', {'a': 'NA', 'b': 'NA'}))
        with self.assertRaises(ValueError):
            propagate(False, {'a'}, ['a'], '0')

    def test_evidence_failure_not_na(self):
        with self.assertRaises(ValueError):
            evaluate(dict(candidate_evidence=[]))
        with self.assertRaises(RuntimeError):
            evaluate(dict(candidate_evidence=[dict(model_fid='arbitrary', exact_complete=False)]))

    def test_missing_hal_not_unmapped(self):
        with self.assertRaises(FileNotFoundError):
            prod.parse_hal_psl_file(str(self.root/'not_present'), label='a', kind='assembly',
                min_coverage=.5, min_identity=0, ambig_identity_delta=.001, ambig_coverage_delta=.01,
                ambig_aligned_bp_delta=10)

    def test_saturated_alignment_is_error(self):
        with patch.object(prod.parasail, 'sg_stats_scan_32', return_value=SimpleNamespace(saturated=True)):
            with self.assertRaises(ArithmeticError):
                prod.parasail_semiglobal_identity('ACGT', 'ACGT')

    def test_main_cli_is_complete_without_old_export(self):
        output = io.StringIO()
        with patch('sys.argv', ['pancgi', '--help']), contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as ex:
            prod.main()
        self.assertEqual(ex.exception.code, 0)
        self.assertNotIn('export-observed-prod', output.getvalue())

    def test_retired_recovery_functions_absent(self):
        for name in ['_tmp_glob_pattern_for_final', '_tmp_suffix_for_path', 'find_resume_pair',
                     'iter_complete_lines_lenient', 'read_tsv_rows_lenient', '_allele_catalog_member_ids',
                     'validate_complete_partial_locus', 'recover_completed_prefix_from_partial',
                     'write_recovered_prefix_rows', 'edlib_global_identity']:
            self.assertFalse(hasattr(prod, name), name)
        self.assertFalse(hasattr(prod.base, 'symmetric_edlib_similarity'))
        self.assertFalse(hasattr(prod.base, 'ensure_parent'))

    def test_retired_cli_options_rejected(self):
        argv = ['pancgi', 'cluster-alleles-prod', '--features', 'x', '--locus-catalog', 'x',
                '--locus-members', 'x', '--out-allele', 'x', '--out-allele-members', 'x']
        for options in [['--resume'], ['--resume-from-allele', 'x'],
                        ['--resume-from-allele-members', 'x'], ['--long-limit', '4000'], ['--kmer', '9']]:
            error = io.StringIO()
            with self.subTest(options=options), patch('sys.argv', argv+options), contextlib.redirect_stderr(error):
                with self.assertRaises(SystemExit) as ex:
                    prod.main()
                self.assertEqual(ex.exception.code, 2)
                self.assertIn('unrecognized arguments', error.getvalue())

    def test_unsupported_artifact_backends_rejected(self):
        argv = ['pancgi', 'export-locus-artifacts-prod', '--features', 'x', '--locus-catalog', 'x',
                '--locus-members', 'x', '--allele-members', 'x', '--out-dir', 'x']
        for options in [['--seq-backend', 'edlib'], ['--seq-backend', 'auto'], ['--very-long-backend', 'auto']]:
            error = io.StringIO()
            with self.subTest(options=options), patch('sys.argv', argv+options), contextlib.redirect_stderr(error):
                with self.assertRaises(SystemExit) as ex:
                    prod.main()
                self.assertEqual(ex.exception.code, 2)
                self.assertIn('invalid choice', error.getvalue())

    def test_fresh_empty_allele_run_does_not_scan_old_outputs(self):
        out_allele, out_members = self.root/'alleles.tsv.gz', self.root/'members.tsv.gz'
        self.write('alleles.tsv.tmp.previous.gz', 'unrelated incomplete data')
        args = SimpleNamespace(out_allele=str(out_allele), out_allele_members=str(out_members),
                               threads=1, parallel_backend='process')
        messages = io.StringIO()
        with patch.object(prod, 'build_allele_tasks_from_args', return_value=([], {'task_loci': 0})), \
             patch.object(prod.glob, 'glob', side_effect=AssertionError('Old output discovery is forbidden')), \
             contextlib.redirect_stderr(messages):
            prod.cmd_cluster_alleles_prod(args)
        for path, header in [(out_allele, prod.ALLELE_CATALOG_HEADER), (out_members, prod.ALLELE_MEMBER_HEADER)]:
            with gzip.open(path, 'rt') as handle:
                self.assertEqual(list(csv.reader(handle, delimiter='\t')), [header])
        self.assertEqual(json.loads(messages.getvalue())['n_alleles'], 0)

    def test_empty_gt_cannot_become_na(self):
        import pandas as pd
        for cell in ['', '.', 'nan', '2']:
            with self.assertRaises(ValueError):
                prod.validate_strict_matrix_compat(pd.DataFrame({'g': [cell]}, index=['a']),
                    label_order=['g'], expected_ids=['a'], name='test')

    def test_corrupt_evidence_is_not_absence(self):
        with self.assertRaises(FileNotFoundError):
            prod.parse_sv_ins_tsv(str(self.root/'missing'))
        for value in ['', '.', '{bad', '{}', '[1]']:
            with self.assertRaises(ValueError):
                prod.parse_sv_overlap_detail_json({'sv_overlap_detail_json': value})
        with self.assertRaises(ValueError):
            prod.hal_fields_for_feature('invalid', {}, kind='reference_primary')
        self.assertEqual(prod.parse_sv_overlap_detail_json({'sv_overlap_detail_json': '[]'}), [])


if __name__ == '__main__':
    unittest.main()
