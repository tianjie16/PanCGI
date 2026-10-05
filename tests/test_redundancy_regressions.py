import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pancgi_pathbed as pathbed
from pancgi_genotyping import pipeline
from pancgi_validation_store import ValidationStore


class RedundancyTests(unittest.TestCase):
    def test_sequential_features_do_not_load_contig_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / 'path.bed'
            bed.write_text('ctg\t0\t100\t>1\nctg\t100\t180\t<2\nctg\t180\t280\t>1\n')
            features = [(30, 140, 'first'), (200, 250, 'second')]
            expected = pathbed._map_features_on_rows('ctg', features, list(pathbed.iter_pathbed(str(bed))), flank_bp=50, flank_max_steps=3)
            with patch.object(pathbed, 'load_pathbed_rows_cached', side_effect=AssertionError('No contig cache')):
                actual = pathbed.map_features_on_contig('ctg', features, str(bed), flank_bp=50, flank_max_steps=3)
            self.assertEqual(actual, expected)

    def test_observed_target_does_not_build_occurrence_index(self):
        with tempfile.TemporaryDirectory() as directory:
            db = Path(directory) / 'models.sqlite'
            sqlite3.connect(db).close()
            task = (dict(label='g', hal_genome='sample'), str(db), ['L1'], {'L1': ['A1', 'A2']}, {'L1'}, {'A1'}, directory)
            with patch.object(pipeline, 'build_target', side_effect=AssertionError('No unobserved locus')):
                label, loci, alleles = pipeline.run_target(task)
            self.assertEqual(loci, {'L1': '1'})
            self.assertEqual(alleles, {'A1': '1', 'A2': '0'})

    def test_disk_index_preserves_order_and_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as directory:
            store = ValidationStore(Path(directory) / 'data.sqlite')
            try:
                values = [{'id': 'b', 'value': '1'}, {'id': 'a', 'value': '2'}]
                rows = store.add(iter(values), 'id', lambda row: None)
                self.assertEqual(list(rows), values)
                self.assertEqual(list(rows.index), ['b', 'a'])
                self.assertEqual(rows.index['b'], values[0])
                with self.assertRaisesRegex(ValueError, 'Duplicate id'):
                    store.add(iter(values * 2), 'id', lambda row: None)
            finally:
                store.close()
