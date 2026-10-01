import argparse
import ast
import copy
import io
import json
import marshal
import random
import subprocess
import sys
import tempfile
import unittest
import zlib
from pathlib import Path

import cpgi_nr_prod as prod
import pancgi_anchor_disk as disk
import test_anchor_bounded as fixture
from test_anchor_bounded import ROOT


class AnchorDiskTests(unittest.TestCase):
    def test_disk_lookup_all_global_weights_and_seed_values(self):
        case = fixture.AnchorBoundedTests()
        case.setUp()
        try:
            reader = case.reader()
            args = argparse.Namespace(in_members=str(case.members), anchor_max_record_bytes=1024*1024)
            path = case.root/'lookup.sqlite'
            metrics = {}
            disk.build_lookup(path, args, reader, metrics)
            lookup = disk.Lookup(path)
            n, frequencies = prod.collect_shingle_df(str(case.db), 3, 8)
            weights = disk.KeyMap(lookup, 'weights', 512)
            expected = prod.shingle_weights(n, frequencies)
            self.assertEqual({k:v.hex() for k,v in weights.items()}, {k:v.hex() for k,v in expected.items()})
            self.assertEqual(list(row[0] for row in lookup.rows('SELECT fid FROM tasks ORDER BY sk COLLATE PYKEY')),
                             sorted(reader, key=prod.anchor_task_sort_key))
            expected_refs = {fid:rec for fid,rec in reader.items() if rec['kind'] == 'reference_primary'}
            self.assertEqual(dict(disk.ReferenceMap(lookup, reader)), expected_refs)
            expected_intervals = prod.build_ref_interval_index(expected_refs)
            intervals = disk.IntervalMap(lookup, 1024*1024)
            for chrom in expected_intervals:
                self.assertEqual(intervals.get(chrom), expected_intervals[chrom])
            self.assertLessEqual(weights.cache.peak, 512)
            lookup.close()
            reader.close()
        finally:
            case.tearDown()

    def test_cli_group_limit_preserves_failure_and_no_outputs(self):
        case = fixture.AnchorBoundedTests()
        case.setUp()
        try:
            out = case.root/'failure'
            out.mkdir()
            result = subprocess.run([sys.executable, str(ROOT/'cpgi_nr_prod.py'), 'anchor-loci-primary',
                '--features', str(case.db), '--in-locus', str(case.loci), '--in-members', str(case.members),
                '--out-locus', str(out/'locus.tsv.gz'), '--out-members', str(out/'members.tsv.gz'),
                '--anchor-group-bytes', '512'], capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            report = json.loads((out/'locus.tsv.gz.anchor.json').read_text())
            self.assertEqual(report['status'], 'failed')
            self.assertIn('MemoryError', report['error'])
            self.assertFalse((out/'locus.tsv.gz').exists())
            self.assertFalse((out/'members.tsv.gz').exists())
            self.assertTrue((Path(report['work_dir'])/'disk_metrics.json').is_file())
        finally:
            case.tearDown()


    def test_frequency_parser_all_boundaries(self):
        pairs = {'a': 1234567890123456789, 'escape\\\"': 1, '\u03b1\U0001f600': 9,
                 'long' * 300: 37, 'empty': 0}
        raw = json.dumps(pairs, ensure_ascii=False)
        for width in (1, 2, 3, 7, 64, len(raw)):
            self.assertEqual(dict(disk.object_pairs((raw[i:i+width] for i in range(0, len(raw), width)), 16384)), pairs)
        self.assertEqual(list(disk.object_pairs(['{', ' ', '}'], 10)), [])
        compressed = io.BytesIO(zlib.compress(raw.encode()))
        self.assertEqual(dict(disk.object_pairs(disk.compressed_text(compressed), 16384)), pairs)

    def test_frequency_parser_rejects_corruption_and_limits(self):
        for raw in ('', '{', '{"a":', '{"a":1,}', '{"a":1}x', '{"a":true}',
                    '{"a":1.5}', '{"a":1e2}', '{1:2}', '[]', '{"a":1 "b":2}'):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                list(disk.object_pairs(iter(raw), 128))
        with self.assertRaises(MemoryError):
            list(disk.object_pairs(iter('{"' + 'x' * 1000 + '":1}'), 128))
        for raw in (zlib.compress(b'{}')[:-1], zlib.compress(b'{}') + b'junk'):
            with self.assertRaises(ValueError):
                list(disk.compressed_text(io.BytesIO(raw)))

    def test_python_sort_preserves_large_integers_and_ties(self):
        keys = [['chr2', -4, 10**80, 'a'], ['chr1', 2, 10**80 + 1, 'b'],
                ['chr1', 2, 10**80, 'b'], ['chr1', 2, 10**80, 'a']]
        with tempfile.TemporaryDirectory() as root:
            db = disk.connect(Path(root) / 'sort.sqlite')
            db.execute('CREATE TABLE t(sk TEXT)')
            db.executemany('INSERT INTO t VALUES (?)', [(disk.key_text(k),) for k in keys])
            self.assertEqual([json.loads(r[0]) for r in db.execute('SELECT sk FROM t ORDER BY sk COLLATE PYKEY')], sorted(keys))
            db.close()

    @staticmethod
    def rows(seed, count=90):
        rng = random.Random(seed)
        rows = []
        for i in range(count):
            start = rng.choice((100, 1000, 10000)) + rng.randrange(300)
            end = start + rng.randrange(-20, 180)
            chrom = rng.choice(('chr1', 'chr2', 'chr10'))
            row = dict(fid=str(i), seed_locus_id=rng.choice(('', 'x', '2', '0002', '99')),
                       anchor_type='nonref', anchor_key='NR::' + str(i), sv_overlap_detail_json='[]',
                       hal_multimap_label='unique', hal_query_cov='0.8',
                       hal_primary_chr=chrom, hal_primary_start1=start, hal_primary_end1=end)
            if i % 4 == 0:
                row['sv_overlap_detail_json'] = json.dumps([dict(id='s', chrom=chrom,
                    pos1=start, vcf_end1=end, svtype='INS', svlen=10)])
            if i % 7 == 0:
                row['hal_multimap_label'] = 'ambiguous_tie'
            if i % 11 == 0:
                row['anchor_type'] = 'primary_ref'
            rows.append(row)
        rng.shuffle(rows)
        return rows

    def test_nonref_clustering_exact_without_global_rows(self):
        for seed in range(20):
            for window in (-1, 0, 100):
                with self.subTest(seed=seed, window=window), tempfile.TemporaryDirectory() as root:
                    args = argparse.Namespace(nonref_site_window_bp=window, hal_nonref_min_coverage=.5,
                        anchor_max_record_bytes=1024*1024, anchor_group_bytes=1024*1024)
                    rows = self.rows(seed)
                    expected = copy.deepcopy(rows)
                    prod.cluster_nonref_site_rows(expected, args)
                    db = disk.create_spool(Path(root) / 'assign.sqlite')
                    for ordinal, row in enumerate(rows, 1):
                        disk.add_assignment(db, ordinal, row, args)
                    metrics = {}
                    disk.cluster_sites(db, args, metrics)
                    actual = [marshal.loads(r[0]) for r in db.execute('SELECT payload FROM assignments ORDER BY ordinal')]
                    self.assertEqual(actual, expected)
                    self.assertLessEqual(metrics['site_component_peak_accounted_bytes'], args.anchor_group_bytes)
                    db.close()

    def test_indivisible_group_fails_without_truncation(self):
        args = argparse.Namespace(nonref_site_window_bp=100, hal_nonref_min_coverage=.5,
            anchor_max_record_bytes=1024*1024, anchor_group_bytes=5000)
        rows = [dict(fid=str(i), anchor_key='NR::1', anchor_type='nonref', seed_locus_id='1',
            sv_overlap_detail_json='[]', hal_multimap_label='unique', hal_query_cov='1', hal_primary_chr='chr1',
            hal_primary_start1=100, hal_primary_end1=200) for i in range(50)]
        with tempfile.TemporaryDirectory() as root:
            db = disk.create_spool(Path(root) / 'assign.sqlite')
            for ordinal, row in enumerate(rows, 1):
                disk.add_assignment(db, ordinal, row, args)
            with self.assertRaisesRegex(MemoryError, 'no cluster truncation'):
                disk.cluster_sites(db, args, {})
            self.assertEqual(db.execute('SELECT count(*) FROM assignments').fetchone()[0], len(rows))
            db.close()

    def test_empty_assignments_keep_headers(self):
        with tempfile.TemporaryDirectory() as root:
            args = argparse.Namespace(out_members=str(Path(root)/'members.tsv'),
                out_locus=str(Path(root)/'locus.tsv'), anchor_group_bytes=1024)
            db = disk.create_spool(Path(root) / 'assign.sqlite')
            metrics = {}
            disk.aggregate_and_write(db, args, {}, {}, metrics, {})
            self.assertEqual(metrics['locus_count'], 0)
            self.assertEqual(len(Path(args.out_locus).read_text().splitlines()), 1)
            db.close()
