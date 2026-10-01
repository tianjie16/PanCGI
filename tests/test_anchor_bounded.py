import argparse
import ast
import csv
import gzip
import hashlib
import json
import multiprocessing
import os
import random
import sqlite3
import subprocess
import sys
import tempfile
import time
import unittest
import zlib
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import cpgi_nr_prod as prod
import pancgi_anchor_bounded as bounded
import pancgi_features as shared
from pancgi_contract import file_state
from test_locus_polish_exact import make_feat, write_jsonl_gz, write_tsv


ROOT = Path(__file__).resolve().parents[1]
SCALING_EVIDENCE = []


def delayed(value):
    time.sleep(.003 * (7 - value % 7))
    return value


def killed(value):
    if value == 2:
        os._exit(23)
    return value


def worker_read(store, fid, connection):
    try:
        record = store[fid]
        store.check()
        connection.send((record, store.pid, store.statistics()))
    finally:
        connection.close()


class AnchorBoundedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / 'features.jsonl.gz'
        self.db = self.root / 'features.sqlite'
        self.records = self.fixture()
        self.prepare(self.records)

    def tearDown(self):
        self.temp.cleanup()

    def fixture(self):
        rows = []
        def add(name, kind='assembly', start='', end='', chrom='chr1', nodes=(1, 2, 3), hal=True):
            fid = f'{name}#0#ctg:0-30'
            row = make_feat(fid, name, kind, nodes,
                chr1=chrom if start != '' else '', start1=start, end1=end)
            row.update(prod.hal_empty_fields('no_hal'))
            row.update(sv_ins_n=0, sv_ins_detail_json='[]', sv_overlap_n=0,
                       sv_overlap_detail_json='[]')
            if start != '' and hal:
                row.update(hal_multimap_label='unique', hal_primary_chr=chrom,
                    hal_primary_start0=start-1, hal_primary_end0=end,
                    hal_primary_start1=start, hal_primary_end1=end,
                    hal_query_cov='0.750000', hal_identity='1.000000',
                    hal_strand='-', hal_admit_rule='unique')
            rows.append(row)
            return row
        add('PRIMARY1', 'reference_primary', 101, 200)
        add('PRIMARY2', 'reference_primary', 151, 250, nodes=(1, 2, 4))
        add('COMPARISON', 'reference_comparison', 101, 130)
        add('HAL_SINGLE', start=101, end=120)
        add('HAL_MULTI', start=160, end=190)
        add('HAL_LOW_COVERAGE', start=1200, end=1230)['hal_query_cov'] = '0.499999'
        add('HAL_AT_COVERAGE', start=1500, end=1530)['hal_query_cov'] = '0.500000'
        add('HAL_AMBIGUOUS', start=101, end=130)['hal_multimap_label'] = 'ambiguous_tie'
        
        add('GRAPH_ONLY', start=101, end=130, hal=False)
        add('NO_SEED')
        for name, position, svtype, svlen, hal in [
            ('INS_REF', 180, 'INS', 60, False),
            ('INS_NR', 3010, 'INS', 80, False),
            ('INS_NR2', 3080, 'INS', 40, False),
            ('DEL_HAL', 4000, 'DEL', -90, True),
            ('DEL_NO_HAL', 4050, 'DEL', -90, False),
        ]:
            row = add(name, start=position, end=position+30, hal=hal)
            detail = dict(id=name, chrom='chr1', pos1=position, vcf_end1=position+10,
                svtype=svtype, svlen=svlen, svlen_abs=abs(svlen), overlap_bp=20, span_bp=30)
            row.update(sv_overlap_n=1, sv_overlap_types=svtype,
                sv_overlap_ids=name, sv_overlap_detail_json=json.dumps([detail]))
            if svtype == 'INS':
                row.update(sv_ins_n=1, sv_ins_detail_json=json.dumps([detail]),
                    sv_ins_longest_chrom='chr1', sv_ins_longest_site1=position,
                    sv_ins_longest_start1=position, sv_ins_longest_end1=position,
                    sv_ins_longest_len=abs(svlen), sv_ins_longest_id=name)
        
        for i, pos in enumerate((7000, 7130, 7260, 7391, 8000)):
            add(f'NR_{i}', start=pos, end=pos+30, nodes=(20, 21, 22))
        add('NR_CHR2', start=7000, end=7030, chrom='chr2')
        
        for i in range(12):
            add(f'SEED_{i}', nodes=(30, 31, 32 if i % 2 else 33))
        return rows

    def prepare(self, records):
        write_jsonl_gz(self.source, records)
        shared.build(self.source, self.db, 3, 8)
        self.loci = self.root / 'seed.tsv.gz'
        self.members = self.root / 'members.tsv.gz'
        write_tsv(self.loci, ['locus_id', 'lead_fid'], [['1', records[0]['fid']]])
        memberships = [[str(200 if r['label'].startswith('SEED_') else i+1), r['fid']]
            for i, r in enumerate(self.records) if r['label'] != 'NO_SEED']
        write_tsv(self.members, ['locus_id', 'fid'], memberships)

    def reader(self, cache=18000, max_record=16*1024*1024):
        return bounded.FeatureMap(self.db, 3, 8, cache, max_record)

    def cli(self, script, name, workers=1, backend='process', seed='17', extra=()):
        out = self.root / name
        out.mkdir()
        command = [sys.executable, str(script), 'anchor-loci-primary',
            '--features', str(self.db), '--in-locus', str(self.loci),
            '--in-members', str(self.members), '--out-locus', str(out / 'locus.tsv.gz'),
            '--out-members', str(out / 'members.tsv.gz'), '--threads', str(workers),
            '--parallel-backend', backend, '--max-medoid-members', '4'] + list(extra)
        env = dict(os.environ, PYTHONHASHSEED=seed)
        result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=120)
        (out / 'stderr.log').write_text(result.stderr)
        self.assertEqual(result.returncode, 0, result.stderr)
        return [gzip.decompress((out / p).read_bytes()) for p in ('locus.tsv.gz', 'members.tsv.gz')]

    def test_all_fields_and_full_global_weights_exact(self):
        store = self.reader()
        try:
            self.assertEqual(list(store), [r['fid'] for r in self.records])
            for r in self.records:
                self.assertEqual(store[r['fid']], prod.prep_feature(r, 3, 8))
            self.assertEqual(prod.collect_shingle_df(str(self.db), 3, 8),
                             prod.collect_shingle_df(str(self.source), 3, 8))
            store.check()
            self.assertLessEqual(store.peak_cache_bytes, store.cache_limit)
            self.assertGreater(store.evictions + store.oversize_uncached, 0)
        finally:
            store.close()

    def test_outputs_exact_serial_serial_thread_process_and_input_order(self):
        expected = self.cli(ROOT / 'cpgi_nr_prod.py', 'serial')
        for workers, backend in [(1, 'process'), (2, 'thread'), (2, 'process'), (8, 'process')]:
            with self.subTest(workers=workers, backend=backend):
                actual = self.cli(ROOT / 'cpgi_nr_prod.py', f'parallel{workers}{backend}',
                    workers, backend, extra=['--anchor-cache-bytes', '18000', '--anchor-max-pending', '3'])
                self.assertEqual(expected, actual)
        self.assertIn(b'hal_multi_ref_best', expected[1])
        self.assertIn(b'sv_insert_multi_ref_best', expected[1])
        self.assertIn(b'nonref_projected_site_cluster', expected[1])
        self.assertIn(b'no_primary_ref_overlap', expected[1])
        self.db.unlink()
        shuffled = list(self.records)
        random.Random(29).shuffle(shuffled)
        self.prepare(shuffled)
        actual = self.cli(ROOT / 'cpgi_nr_prod.py', 'shuffled', 2,
            extra=['--anchor-cache-bytes', '1', '--anchor-max-pending', '1'])
        self.assertEqual(expected, actual)


        for seed in ('0', '999'):
            serial = self.cli(ROOT / 'cpgi_nr_prod.py', 'serial_seed' + seed, seed=seed)
            actual = self.cli(ROOT / 'cpgi_nr_prod.py', 'parallel_seed' + seed, 2, seed=seed)
            self.assertEqual(serial, actual)

    def test_missing_identity_and_raw_limit_fail_closed(self):
        store = self.reader(max_record=40)
        try:
            with self.assertRaises(KeyError):
                store['missing']
            with self.assertRaises(MemoryError):
                store[self.records[0]['fid']]
            self.assertEqual(store.loads, 0)
        finally:
            store.close()

    def test_cache_residency_does_not_grow_with_feature_count(self):
        for count in (64, 8192):
            source, index = self.root / f'scale{count}.gz', self.root / f'scale{count}.sqlite'
            def rows():
                for i in range(count):
                    row = dict(self.records[-1], fid=f'SCALE{i}#0#ctg:0-16384',
                               source_seq='ACGT' * 4096, source_seq_len=16384)
                    yield row
            write_jsonl_gz(source, rows())
            shared.build(source, index, 3, 8)
            store = bounded.FeatureMap(index, 3, 8, 65536, 65536)
            try:
                seen = 0
                for fid in store:
                    self.assertEqual(store[fid]['source_seq'], 'ACGT' * 4096)
                    seen += 1
                    self.assertLessEqual(store.cache_bytes, 65536)
                self.assertEqual(seen, count)
                self.assertEqual(store.loads, count)
                self.assertGreater(store.evictions, 0)
                self.assertLessEqual(store.peak_cache_records, 3)
                store.check()
                SCALING_EVIDENCE.append(dict(records=count, **store.statistics()))
            finally:
                store.close()

    def test_state_parameter_and_device_checks_not_bypassed(self):
        with self.assertRaisesRegex(ValueError, 'parameter mismatch'):
            bounded.FeatureMap(self.db, 4, 8, 1000, 100000)
        store = self.reader()
        try:
            stat = self.source.stat()
            os.utime(self.source, ns=(stat.st_atime_ns, stat.st_mtime_ns+1000000000))
            with self.assertRaisesRegex(ValueError, 'Source features changed'):
                store.check()
        finally:
            store.close()
        with sqlite3.connect(self.db) as db:
            meta = json.loads(db.execute('SELECT payload FROM shared_metadata').fetchone()[0])
            meta['source_state'] = file_state(self.source)
            meta['source_state'][0] += 1
            db.execute('UPDATE shared_metadata SET payload=?', (json.dumps(meta),))
        with self.assertRaisesRegex(ValueError, 'Source features changed'):
            self.reader()

    def test_index_raw_corruption_and_count_mismatch_rejected(self):
        fid = self.records[0]['fid']
        wrong = dict(self.records[0], fid='wrong')
        with sqlite3.connect(self.db) as db:
            db.execute('UPDATE polish_feature_store SET raw=? WHERE fid=?',
                (zlib.compress(json.dumps(wrong).encode()), fid))
        store = self.reader()
        try:
            with self.assertRaisesRegex(ValueError, 'identity mismatch'):
                store[fid]
            with sqlite3.connect(self.db) as db:
                db.execute('DELETE FROM polish_feature_store WHERE fid=?', (fid,))
            with self.assertRaisesRegex(RuntimeError, 'index changed'):
                store.check()
        finally:
            store.close()
        with self.assertRaisesRegex(ValueError, 'count mismatch'):
            self.reader()

    def test_worker_local_connection_and_parent_cache_independent(self):
        store = self.reader()
        try:
            fid = self.records[0]['fid']
            store[fid]
            context = multiprocessing.get_context('fork')
            parent, child = context.Pipe()
            process = context.Process(target=worker_read, args=(store, fid, child))
            process.start()
            child.close()
            self.assertTrue(parent.poll(15))
            record, pid, stats = parent.recv()
            process.join(15)
            self.assertEqual(process.exitcode, 0)
            self.assertNotEqual(pid, os.getpid())
            self.assertEqual(stats['loads'], 1)
            self.assertEqual(record, store[fid])
            self.assertEqual(store.pid, os.getpid())
            store.check()
        finally:
            store.close()

    def test_ordered_submission_bound_and_worker_death(self):
        for backend in ('thread', 'process'):
            args = argparse.Namespace(threads=2, anchor_max_pending=3, parallel_backend=backend,
                                      mp_start_method='fork', maxtasksperchild=0)
            seen = []
            def tasks():
                for i in range(30):
                    seen.append(i)
                    yield i
            output = []
            for result in bounded.ordered_assignments(delayed, tasks(), args):
                output.append(result)
                self.assertLessEqual(len(seen) - len(output), 2)
            self.assertEqual(output, list(range(30)))
        with self.assertRaises(BrokenProcessPool):
            list(bounded.ordered_assignments(killed, range(10), args))
        self.assertEqual(multiprocessing.active_children(), [])

    def test_publication_waits_for_input_recheck_and_no_overwrite(self):
        args = argparse.Namespace(features=str(self.db), in_locus=str(self.loci), in_members=str(self.members),
            out_locus=str(self.root / 'out.tsv.gz'), out_members=str(self.root / 'out.members.tsv.gz'),
            anchor_resource_report='', anchor_k=3, max_mid_anchors=8, anchor_cache_bytes=18000,
            anchor_max_record_bytes=100000, anchor_max_pending=1, threads=1)
        def changed(staged, store):
            Path(staged.out_locus).write_bytes(b'locus')
            Path(staged.out_members).write_bytes(b'members')
            with self.members.open('ab') as handle:
                handle.write(b'changed')
        with self.assertRaisesRegex(RuntimeError, 'Anchor input changed'):
            bounded.run(args, changed)
        self.assertFalse(Path(args.out_locus).exists())
        self.assertFalse(Path(args.out_members).exists())
        receipt = json.loads(Path(args.out_locus + '.anchor.json').read_text())
        self.assertEqual(receipt['status'], 'failed')
        with self.assertRaises(FileExistsError):
            bounded.run(args, changed)


if __name__ == '__main__':
    unittest.main()
