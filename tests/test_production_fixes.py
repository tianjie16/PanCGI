import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import cpgi_nr_prod as prod
import cpgi_sv_annot as sv
import pancgi_mapping as mapping
from pancgi_contract import read_bed
import wfa_longalign_wrapper as wfa
from pancgi_annotations import mechanism
from pancgi_results import relative_strand
from pancgi_merged_schema import merged_dataframe
import pancgi_hal_runtime as hal_runtime


class ProductionFixTests(unittest.TestCase):
    def test_native_runtime_is_explicit(self):
        with patch.dict('os.environ', {'PANCGI_HAL_RUNTIME':'native'}), patch('shutil.which', return_value='/tools/halStats'):
            self.assertEqual(hal_runtime.container_command('docker','',[], 'halStats',['--version']), ['/tools/halStats','--version'])
            with self.assertRaises(ValueError):
                hal_runtime.container_command('docker','image',[], 'halStats',[])

    def test_merged_numeric_null_schema(self):
        columns = ['allele_id','rep_hal_primary_start0','locus_primary_start_median1','hap']
        df = merged_dataframe([dict(zip(columns, ['None',10,10.5,'NA'])), dict(zip(columns, ['null','','','0']))], columns, ['hap'])
        self.assertEqual(str(df.dtypes['rep_hal_primary_start0']), 'Int64')
        self.assertEqual(str(df.dtypes['locus_primary_start_median1']), 'Float64')
        self.assertEqual(df['allele_id'].tolist(), ['None','null'])
        self.assertEqual(df['hap'].tolist(), ['NA','0'])
        self.assertTrue(df['rep_hal_primary_start0'].isna().iloc[1])

    def test_psl_relative_strand(self):
        for raw, expected in [('', ''), ('+', '+'), ('-', '-'), ('++', '+'), ('--', '+'), ('+-', '-'), ('-+', '-')]:
            self.assertEqual(relative_strand(raw), expected)
        with self.assertRaises(ValueError):
            relative_strand('invalid')

    def test_complete_mechanism_vocabulary(self):
        ins = lambda c: dict(svtype='INS', **{'class':c})
        deletion = dict(svtype='DEL', **{'class':'contains_del'})
        for events, label in [([], 'non-SV'), ([ins('inside_ins')], 'INS-internal'),
                               ([ins('partial_ins')], 'INS-junction'), ([ins('contains_ins')], 'INS-spanning'),
                               ([deletion], 'DEL-junction'), ([ins('contains_ins'),deletion], 'INS+DEL'),
                               ([ins('contains_ins'),ins('partial_ins'),ins('inside_ins')], 'INS-internal')]:
            self.assertEqual(mechanism(events)['mechanism_class'], label)
            self.assertEqual(mechanism(list(reversed(events)))['mechanism_class'], label)
        self.assertEqual(mechanism([], 'primary_reference')['mechanism_class'], '')
        self.assertEqual(mechanism([], 'primary_reference')['mechanism_status'], 'not_applicable_reference_role')

    def test_32bit_long_identity(self):
        self.assertEqual(prod.parasail_semiglobal_identity('ACGT'*4500, 'ACGT'*4500), 1.0)

    def test_parasail_short_precision_equivalence(self):
        matrix = prod.parasail_matrix(2, -3)
        for a, b in [('ACGT'*30, 'ACGT'*30), ('ACGT'*30, 'ACGT'*10+'AAAA'+'ACGT'*20),
                     ('ACGT'*30, 'ACGT'*29), ('AACGTACGT', 'ACGTTCGT')]:
            for x, y in [(a, b), (b, a)]:
                old = prod.parasail.sg_stats_scan_16(x, y, 5, 2, matrix)
                self.assertFalse(old.saturated)
                self.assertEqual(prod.parasail_semiglobal_identity(x, y), old.matches/old.length)

    def test_wfa_declared_cigar_convention(self):
        for a, b in [('ACGT', 'ACAGT'), ('ACAGT', 'ACGT'), ('ACTT', 'ACGT')]:
            result = wfa.alignment_identity(a, b)
            self.assertEqual(result['status'], 0)
            self.assertEqual(result['matches'], 4 if len(a) != len(b) else 3)
            self.assertEqual(result['alignment_length'], max(len(a), len(b)))

    def test_bad_wfa_cigar_rejects(self):
        for a, b, cigar in [('AAAA', 'TTTT', '4='), ('AAAA', 'AAAA', '4X'),
                             ('A', 'AA', '1M1I'), ('A', 'A', '2M'), ('A', 'A', '1S'),
                             ('A', 'A', '0M'), ('AA', 'AA', '1M')]:
            with self.subTest(cigar=cigar), self.assertRaises((ValueError, RuntimeError)):
                wfa.replay_cigar_identity(a, b, cigar)

    def test_failed_wfa_does_not_emit_identity(self):
        with patch.object(wfa, 'WavefrontAligner') as factory:
            factory.return_value.status = -1
            with self.assertRaisesRegex(RuntimeError, 'status'):
                wfa.alignment_identity('A', 'A')

    def test_literal_sv_ids_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'sv.tsv'
            path.write_text('ID\tVARID\tcontig\tstart\tend\tSVTYPE\t#CHROM\tPOS\tVCF_END\tSVLEN\n'
                            'None\tNone\tc\t10\t30\tINS\tr\t20\t20\t20\n'
                            'null\tnull\tc\t10\t30\tINS\tr\t20\t20\t20\n')
            result = sv.annotate_feature_overlaps('c', 15, 25, sv.build_sv_index(str(path)))
            self.assertEqual(result['sv_overlap_n'], 2)
            self.assertEqual(set(result['sv_overlap_ids'].split(';')), {'None', 'null'})

    def test_bed_metric_consistency(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'cgi.bed'
            for values, valid in [('2\t12\t20\t60', True), ('0\t0\t99\t99', False),
                                  ('2\t12\t20.01\t60', False)]:
                path.write_text(f'c\t0\t20\tid\t20\t{values}\t1\n')
                if valid:
                    self.assertEqual(len(list(read_bed(path, {'c':20}))), 1)
                else:
                    with self.assertRaises(ValueError):
                        list(read_bed(path, {'c':20}))

    def test_inventory_does_not_certify_topology(self):
        for path_record in ['P\tx\t1+,2+\t0M\n', 'P\tx\t1+,2+\t*\n',
                            'W\ts\t1\tc\t0\t4\t>1>2\n',
                            'W\ts\t1\tc\t0\t2\t>1\nW\ts\t1\tc\t2\t4\t>2\n']:
            for edge in ['', 'L\t1\t+\t2\t+\t0M\n', 'L\t1\t+\t2\t-\t0M\n']:
                with tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    (root/'g.gfa').write_text('S\t1\tAC\nS\t2\tGT\n'+edge+path_record)
                    mapping.inventory_gfa(str(root/'g.gfa'), str(root/'paths.tsv'), str(root/'index'))
                    row = mapping.read_tsv(str(root/'paths.tsv'), mapping.GFA_INVENTORY_COLUMNS)[0]
                    self.assertTrue(row['gfa_path_id'].startswith('GP_'))
                    self.assertFalse(list(root.rglob('*.sqlite')))

    def test_parquet_no_type_retry(self):
        import pandas as pd
        with tempfile.TemporaryDirectory() as tmp, patch.object(pd.DataFrame, 'to_parquet', side_effect=ValueError('bad type')) as save:
            with self.assertRaises(ValueError):
                prod.write_optional_parquet(pd.DataFrame({'x': [1, 'bad']}), str(Path(tmp)/'out.parquet'))
            self.assertEqual(save.call_count, 1)

    def test_default_production_profile(self):
        argv = ['pancgi', 'cluster-alleles-prod', '--features', 'x', '--locus-catalog', 'x',
                '--locus-members', 'x', '--out-allele', 'x', '--out-allele-members', 'x']
        with patch('sys.argv', argv), patch.object(prod, 'cmd_cluster_alleles_prod') as command:
            prod.main()
        args = command.call_args.args[0]
        self.assertEqual((args.identity, args.min_len_ratio, args.allele_cluster_mode), (.8, .8, 'allpairs_clique'))


if __name__ == '__main__':
    unittest.main()
