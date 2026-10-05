import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import tempfile
import unittest
import zlib
from contextlib import closing
from pathlib import Path
from unittest.mock import patch

import pancgi_anchor_bounded as anchor
import pancgi_features as store
from pancgi_contract import file_state
from test_locus_polish_exact import make_feat


class FeatureIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / 'features.jsonl.gz'
        self.index = self.root / 'features.sqlite'
        self.records = [make_feat(f'm{i}', 'sample', 'assembly', nodes)
                        for i, nodes in enumerate(([1, 2, 3], [1, 2, 4], [7, 8, 9]), 1)]
        with gzip.open(self.source, 'wt', encoding='utf-8', newline='') as handle:
            for rec in self.records:
                handle.write(json.dumps(rec) + '\n')
        store.build(self.source, self.index, 3, 8)
        store._verify_source.cache_clear()

    def tearDown(self):
        store._verify_source.cache_clear()
        self.temp.cleanup()

    def metadata(self):
        with sqlite3.connect(self.index) as db:
            return json.loads(db.execute('SELECT payload FROM shared_metadata').fetchone()[0])

    def update(self, change):
        with sqlite3.connect(self.index) as db:
            info = json.loads(db.execute('SELECT payload FROM shared_metadata').fetchone()[0])
            change(info)
            db.execute('UPDATE shared_metadata SET payload=?', (json.dumps(info),))

    def different_filesystem(self, unhashed=False):
        def change(info):
            info['source_state'][0] += 100
            if unhashed:
                info['schema_version'] = 1
                info.pop('source_sha256')
        self.update(change)

    def rewrite_records(self, records):
        with gzip.open(self.source, 'wt', encoding='utf-8', newline='') as handle:
            for rec in records:
                handle.write(json.dumps(rec) + '\n')

    def touch(self, path):
        state = path.stat()
        os.utime(path, ns=(state.st_atime_ns, state.st_mtime_ns + 1000000000))

    def test_build_hashes_exact_compressed_bytes(self):
        info = self.metadata()
        self.assertEqual(info['schema_version'], 2)
        self.assertEqual(info['source_sha256'], hashlib.sha256(self.source.read_bytes()).hexdigest())

    def test_plain_text_hash_crlf_blank_lines_and_unicode(self):
        source, index = self.root / 'plain.jsonl', self.root / 'plain.sqlite'
        records = [dict(self.records[0], source_seq='ACGT', genome='sample_\u03b1')]
        raw = ('\r\n' + json.dumps(records[0], ensure_ascii=False) + '\r\n\r\n').encode()
        source.write_bytes(raw)
        info = store.build(source, index, 3, 8)
        self.assertEqual(info['source_sha256'], hashlib.sha256(raw).hexdigest())
        self.assertEqual(list(store.iter_raw(index)), records)

    def test_cross_device_hash_checked_once_and_source_and_index_unchanged(self):
        self.different_filesystem()
        before = {p: (file_state(p), p.read_bytes()) for p in (self.source, self.index)}
        with patch.object(store, 'digest', wraps=store.digest) as digest:
            self.assertEqual(list(store.iter_raw(self.index)), self.records)
            self.assertEqual(list(store.selected_raw(self.index, ['m3', 'm1'])), ['m1', 'm3'])
            self.assertEqual(store.frequencies(self.index, 3, 8)[0], 3)
            reader = anchor.FeatureMap(self.index, 3, 8, 65536, 65536)
            try:
                reader['m1']
                reader.check()
                reader.check()
            finally:
                reader.close()
            self.assertEqual(digest.call_count, 1)
        self.assertEqual(before, {p: (file_state(p), p.read_bytes()) for p in before})

    def test_matching_environment_does_not_rescan_source(self):
        with patch.object(store, 'digest', side_effect=AssertionError('unnecessary scan')):
            self.assertEqual(list(store.iter_raw(self.index)), self.records)

    def test_same_bytes_new_inode_and_timestamps_accepted(self):
        copy = self.root / 'replacement.gz'
        shutil.copyfile(self.source, copy)
        copy.replace(self.source)
        self.touch(self.source)
        self.assertEqual(list(store.iter_raw(self.index)), self.records)

    def test_wrong_content_rejected_without_rebuilding(self):
        original = self.index.read_bytes()
        self.rewrite_records([dict(self.records[0], kind='changed'), *self.records[1:]])
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            store.frequencies(self.index, 3, 8)
        self.assertEqual(original, self.index.read_bytes())

    def test_same_size_content_change_with_restored_mtime_rejected(self):
        original_state = self.source.stat()
        raw = bytearray(self.source.read_bytes())
        raw[-5] ^= 1
        self.source.write_bytes(raw)
        os.utime(self.source, ns=(original_state.st_atime_ns, original_state.st_mtime_ns))
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            list(store.iter_raw(self.index))

    def test_missing_hash_rejected_without_unhashed_downgrade(self):
        self.update(lambda info: info.pop('source_sha256'))
        with self.assertRaisesRegex(ValueError, 'SHA256'):
            store.frequencies(self.index, 3, 8)

    def test_unknown_schema_rejected(self):
        self.update(lambda info: info.update(schema_version=99))
        with self.assertRaisesRegex(ValueError, 'Unsupported'):
            store.frequencies(self.index, 3, 8)

    def test_unhashed_device_difference_compares_all_records_once(self):
        self.different_filesystem(unhashed=True)
        original = self.index.read_bytes()
        with patch.object(store, '_compare_source_records', wraps=store._compare_source_records) as compare:
            self.assertEqual(list(store.iter_raw(self.index)), self.records)
            self.assertEqual(store.frequencies(self.index, 3, 8)[0], 3)
            self.assertEqual(compare.call_count, 1)
        self.assertEqual(original, self.index.read_bytes())

    def test_unhashed_changed_source_order_rejected(self):
        self.different_filesystem(unhashed=True)
        self.rewrite_records(self.records[::-1])
        with self.assertRaisesRegex(ValueError, 'content mismatch'):
            list(store.iter_raw(self.index))

    def test_unhashed_missing_record_rejected(self):
        self.different_filesystem(unhashed=True)
        self.rewrite_records(self.records[:-1])
        with self.assertRaisesRegex(ValueError, 'count mismatch'):
            list(store.iter_raw(self.index))

    def test_unhashed_extra_record_rejected(self):
        self.different_filesystem(unhashed=True)
        self.rewrite_records(self.records + [dict(self.records[0], fid='extra')])
        with self.assertRaisesRegex(ValueError, 'count mismatch'):
            list(store.iter_raw(self.index))

    def test_unhashed_fid_corruption_rejected(self):
        self.different_filesystem(unhashed=True)
        with sqlite3.connect(self.index) as db:
            db.execute("UPDATE polish_feature_store SET fid='wrong' WHERE ordinal=1")
        with self.assertRaisesRegex(ValueError, 'content mismatch'):
            list(store.iter_raw(self.index))

    def test_unhashed_trailing_compressed_payload_rejected(self):
        self.different_filesystem(unhashed=True)
        with sqlite3.connect(self.index) as db:
            raw = db.execute('SELECT raw FROM polish_feature_store WHERE ordinal=1').fetchone()[0]
            db.execute('UPDATE polish_feature_store SET raw=? WHERE ordinal=1', (raw + zlib.compress(b'extra'),))
        with self.assertRaisesRegex(ValueError, 'content mismatch'):
            list(store.iter_raw(self.index))

    def test_source_mutation_during_hash_rejected(self):
        self.different_filesystem()
        digest = store.digest
        def mutate(path):
            value = digest(path)
            self.touch(self.source)
            return value
        with patch.object(store, 'digest', side_effect=mutate):
            with self.assertRaisesRegex(ValueError, 'changed during identity verification'):
                list(store.iter_raw(self.index))

    def test_source_mutation_during_unhashed_check_rejected(self):
        self.different_filesystem(unhashed=True)
        compare = store._compare_source_records
        def mutate(db, info):
            compare(db, info)
            self.touch(self.source)
        with patch.object(store, '_compare_source_records', side_effect=mutate):
            with self.assertRaisesRegex(ValueError, 'changed during identity verification'):
                list(store.iter_raw(self.index))

    def test_cached_verification_does_not_hide_later_change(self):
        self.different_filesystem()
        list(store.iter_raw(self.index))
        self.rewrite_records(self.records[::-1])
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            list(store.iter_raw(self.index))

    def test_device_change_within_reader_still_rejected(self):
        with closing(store.connect(self.index)) as db:
            store.metadata(db)
            state = file_state
            def changed(path):
                result = state(path)
                if str(path) == str(self.source):
                    result[0] += 1
                return result
            with patch.object(store, 'file_state', side_effect=changed):
                with self.assertRaisesRegex(ValueError, 'changed during reading'):
                    store.metadata(db)

    def test_source_touched_during_iteration_rejected(self):
        rows = store.iter_raw(self.index)
        next(rows)
        self.touch(self.source)
        with self.assertRaisesRegex(ValueError, 'changed during reading'):
            list(rows)

    def test_index_modified_during_iteration_rejected(self):
        rows = store.iter_raw(self.index)
        next(rows)
        self.touch(self.index)
        with self.assertRaisesRegex(RuntimeError, 'index changed'):
            list(rows)

    def test_anchor_check_rejects_source_change_after_cache_hit(self):
        self.different_filesystem()
        reader = anchor.FeatureMap(self.index, 3, 8, 65536, 65536)
        try:
            reader['m1']
            self.touch(self.source)
            with self.assertRaisesRegex(ValueError, 'changed during reading'):
                reader.check()
        finally:
            reader.close()

    def test_unhashed_wal_change_invalidates_verification_cache(self):
        self.different_filesystem(unhashed=True)
        with closing(sqlite3.connect(self.index)) as writer:
            writer.execute('PRAGMA journal_mode=WAL')
            writer.execute('PRAGMA wal_autocheckpoint=0')
            list(store.iter_raw(self.index))
            writer.execute("UPDATE polish_feature_store SET fid='wrong' WHERE ordinal=1")
            writer.commit()
            with self.assertRaisesRegex(ValueError, 'content mismatch'):
                list(store.iter_raw(self.index))


if __name__ == '__main__':
    unittest.main()
