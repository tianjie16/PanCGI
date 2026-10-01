import argparse
import ast
import gzip
import importlib.util
import json
import marshal
import os
import sqlite3
import unittest
import zlib
from pathlib import Path
from unittest.mock import patch

import cpgi_nr_prod as prod
import pancgi_anchor_bounded as bounded
import pancgi_anchor_disk as disk
import pancgi_anchor_tables as tables
import test_anchor_bounded as fixture
from test_anchor_bounded import ROOT


class AnchorMetadataTests(unittest.TestCase):
    def setUp(self):
        self.case = fixture.AnchorBoundedTests()
        self.case.setUp()

    def tearDown(self):
        self.case.tearDown()

    def mutate_raw(self, transform):
        fid = self.case.records[0]['fid']
        with sqlite3.connect(self.case.db) as db:
            raw = db.execute('SELECT raw FROM polish_feature_store WHERE fid=?', (fid,)).fetchone()[0]
            db.execute('UPDATE polish_feature_store SET raw=? WHERE fid=?', (transform(raw), fid))
        return fid

    def test_stream_is_ordered_complete_and_never_prepares_or_fills_cache(self):
        reader = self.case.reader(cache=1)
        try:
            with patch.object(prod, 'prep_feature', side_effect=AssertionError('unexpected full preparation')):
                actual = list(reader.iter_metadata())
            expected = [(r['fid'], {k:r[k] for k in bounded.METADATA_FIELDS if k in r})
                        for r in self.case.records]
            self.assertEqual(actual, expected)
            self.assertEqual(reader.loads, 0)
            self.assertEqual(reader.raw_loads, len(expected))
            self.assertEqual(reader.metadata_reads, len(expected))
            self.assertEqual(reader.cache_bytes, 0)
            self.assertEqual(reader[self.case.records[0]['fid']], prod.prep_feature(self.case.records[0], 3, 8))
            self.assertEqual(reader.loads, 1)
            reader.check()
        finally:
            reader.close()

    def test_stream_is_one_ordered_query_without_per_fid_selects(self):
        reader = self.case.reader()
        statements = []
        try:
            reader.conn.set_trace_callback(statements.append)
            list(reader.iter_metadata())
            scans = [s for s in statements if 'FROM polish_feature_store' in s]
            self.assertEqual(len(scans), 1)
            self.assertIn('ORDER BY ordinal', scans[0])
            self.assertNotIn('WHERE fid=', scans[0])
        finally:
            reader.close()

    def test_corruption_rejected_by_both_access_paths(self):
        original = zlib.compress(json.dumps(self.case.records[0]).encode())
        wrong = dict(self.case.records[0], fid='wrong')
        examples = [(zlib.compress(json.dumps(wrong).encode()), ValueError),
                    (original + b'trailing', ValueError), (original[:-1], MemoryError),
                    (zlib.compress(b'{bad json}'), ValueError), (b'not zlib', zlib.error)]
        fid = self.case.records[0]['fid']
        for payload, error in examples:
            with self.subTest(error=error, payload=payload[:12]):
                self.mutate_raw(lambda _: payload)
                for metadata in (True, False):
                    reader = self.case.reader()
                    try:
                        with self.assertRaises(error):
                            if metadata:
                                list(reader.iter_metadata())
                            else:
                                reader[fid]
                    finally:
                        reader.close()

    def test_size_limits_apply_to_stream_and_point_access(self):
        fid = self.case.records[0]['fid']
        for payload, limit, message in [(b'x' * 4096, 40, 'compressed'),
                (zlib.compress(json.dumps(self.case.records[0]).encode()), 40, 'raw record')]:
            self.mutate_raw(lambda _: payload)
            for metadata in (True, False):
                reader = self.case.reader(max_record=limit)
                try:
                    with self.assertRaisesRegex(MemoryError, message):
                        list(reader.iter_metadata()) if metadata else reader[fid]
                finally:
                    reader.close()

    def test_missing_null_zero_and_coordinate_semantics(self):
        fid = self.case.records[0]['fid']
        variants = [dict(primary_chr='', primary_start1=None, primary_end1='.'),
                    dict(primary_chr='chr10', primary_start1='001', primary_end1=30),
                    dict(primary_chr=None, primary_start1=0, primary_end1=0, nonref_only=None)]
        for i, fields in enumerate(variants):
            record = dict(self.case.records[0], **fields)
            if i == 0:
                record.pop('nonref_only', None)
            self.mutate_raw(lambda _: zlib.compress(json.dumps(record).encode()))
            reader = self.case.reader()
            try:
                meta = next(reader.iter_metadata())[1]
                self.assertEqual(prod.ref_interval1_from_feat(meta), prod.ref_interval1_from_feat(reader[fid]))
                self.assertEqual('nonref_only' in meta, 'nonref_only' in record)
                self.assertEqual(meta, {k:record[k] for k in bounded.METADATA_FIELDS if k in record})
            finally:
                reader.close()

    def test_stream_end_rechecks_source_state(self):
        reader = self.case.reader()
        try:
            iterator = reader.iter_metadata()
            next(iterator)
            state = self.case.source.stat()
            os.utime(self.case.source, ns=(state.st_atime_ns, state.st_mtime_ns+1000000))
            with self.assertRaisesRegex(ValueError, 'Source features changed'):
                list(iterator)
        finally:
            reader.close()


    def test_member_writer_only_consumes_projected_fields(self):
        tree = ast.parse((ROOT/'pancgi_anchor_tables.py').read_text())
        func = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'member_values')
        used = set()
        for n in ast.walk(func):
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and isinstance(n.func.value, ast.Name) and n.func.value.id == 'feat':
                self.assertEqual(n.func.attr, 'get')
                self.assertIsInstance(n.args[0], ast.Constant)
                used.add(n.args[0].value)
            if isinstance(n, ast.Subscript) and isinstance(n.value, ast.Name):
                self.assertNotEqual(n.value.id, 'feat')
        self.assertEqual(used, {'primary_chr', 'primary_start1', 'primary_end1', 'nonref_only'})
        self.assertTrue(used <= set(bounded.METADATA_FIELDS))

    def test_serial_parallel_outputs_and_metadata_access(self):
        baseline = ROOT
        expected = self.case.cli(baseline/'cpgi_nr_prod.py', 'serial', workers=1)
        actual = self.case.cli(ROOT/'cpgi_nr_prod.py', 'parallel', workers=4,
                               extra=['--anchor-cache-bytes', '1', '--anchor-max-pending', '2'])
        self.assertEqual(expected, actual)
        out = self.case.root/'parallel'
        receipt = json.loads((out/'locus.tsv.gz.anchor.json').read_text())
        metrics = receipt['disk_pipeline']
        self.assertEqual(metrics['lookup_prepared_loads'], 0)
        self.assertEqual(metrics['member_output_prepared_loads'], 0)
        self.assertEqual(metrics['metadata_records'], len(self.case.records))
        self.assertEqual(metrics['member_output_records'], len(self.case.records))
        events = [json.loads(line) for line in Path(receipt['progress_path']).read_text().splitlines()]
        for phase in ('metadata', 'lookup_sort', 'assignment', 'site_clustering', 'locus_science', 'member_output'):
            self.assertTrue(any(e['phase'] == phase and e['event'] == 'complete' for e in events), phase)
        metadata = next(e for e in events if e['phase'] == 'metadata' and e['event'] == 'complete')
        self.assertEqual(metadata['parent_feature_access']['loads'], 0)
        begin = next(e for e in events if e['phase'] == 'member_output' and e['event'] == 'start')
        end = next(e for e in events if e['phase'] == 'member_output' and e['event'] == 'complete')
        self.assertEqual(begin['parent_feature_access']['raw_loads'], end['parent_feature_access']['raw_loads'])


if __name__ == '__main__':
    unittest.main()
