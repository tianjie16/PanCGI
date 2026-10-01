import gzip
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import cpgi_nr_prod as prod
import pancgi_features as store
from test_locus_polish_exact import make_feat


class SharedFeatureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source, self.db = self.root / 'features.jsonl.gz', self.root / 'features.sqlite'
        self.records = [make_feat(f'm{i}', 'sample', 'assembly', nodes)
            for i, nodes in enumerate(([1,2,3], [1,2,4], [1,2,3], [7,8,9]), 1)]
        with gzip.open(self.source, 'wt') as f:
            for rec in self.records:
                f.write(json.dumps(rec) + '\n')
        store.build(self.source, self.db, 3, 8)

    def tearDown(self):
        self.temp.cleanup()

    def test_weights_and_original_order_exact(self):
        before = prod.collect_shingle_df(str(self.source), 3, 8)
        after = prod.collect_shingle_df(str(self.db), 3, 8)
        self.assertEqual(before, after)
        self.assertEqual(prod.shingle_weights(*before), prod.shingle_weights(*after))
        self.assertEqual(list(prod.iter_feature_jsonl(self.db)), self.records)
        self.assertEqual(prod.load_prepped_feature_index(str(self.source), 3, 8), prod.load_prepped_feature_index(str(self.db), 3, 8))

    def test_subset_keeps_full_universe(self):
        before = prod.load_prepped_feature_index_selected(str(self.source), {'m4', 'm1'}, 3, 8)
        after = prod.load_prepped_feature_index_selected(str(self.db), {'m4', 'm1'}, 3, 8)
        self.assertEqual(list(before), list(after))
        self.assertEqual(before, after)
        self.assertEqual(store.frequencies(self.db, 3, 8)[0], 4)

    def test_polish_payload_unchanged(self):
        args = SimpleNamespace(features=str(self.db), anchor_k=3, max_mid_anchors=8)
        self.assertEqual(prod.load_or_build_polish_feature_store(args)[0], str(self.db))
        with prod.PolishSQLiteFeatureStore(str(self.db)) as reader:
            for rec in self.records:
                p = prod.prep_feature(rec, 3, 8)
                expected = prod._polish_store_unpack_record(rec['fid'], prod._polish_store_pack_record(p))
                self.assertEqual(reader.get(rec['fid']), expected)

    def test_mismatched_parameters_fail(self):
        with self.assertRaisesRegex(ValueError, 'parameter mismatch'):
            store.frequencies(self.db, 5, 8)

    def test_source_change_fails_not_rebuilt(self):
        self.source.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'changed'):
            store.frequencies(self.db, 3, 8)
