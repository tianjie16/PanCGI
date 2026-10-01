import unittest
from urllib.parse import unquote

from pancgi_identifiers import allele_id_map, representative_allele_id


class PublicIdentifierTests(unittest.TestCase):
    def member(self, genome='HG002_hap1', contig='chr1', start=950, end=1350):
        return dict(hal_genome=genome, hal_sequence=contig, start0=start, end0=end)

    def test_literal_haplotype_and_interval(self):
        self.assertEqual(representative_allele_id(self.member()), 'HG002_hap1_chr1:950:1350')

    def test_same_contig_in_different_genomes(self):
        one = representative_allele_id(self.member())
        two = representative_allele_id(self.member(genome='HG002_hap2'))
        self.assertNotEqual(one, two)

    def test_names_are_not_split_or_shortened(self):
        row = self.member(genome='sample_1_hap_2', contig='sample_1_hap_2_contig_3')
        self.assertEqual(representative_allele_id(row), 'sample_1_hap_2_sample_1_hap_2_contig_3:950:1350')

    def test_reserved_characters(self):
        row = self.member(genome='g:1% x', contig='c/#>2')
        self.assertEqual(representative_allele_id(row), 'g%3A1%25%20x_c%2F%23%3E2:950:1350')

    def test_literal_escape_text_is_distinct(self):
        one = representative_allele_id(self.member(contig='c:1'))
        two = representative_allele_id(self.member(contig='c%3A1'))
        self.assertNotEqual(one, two)

    def test_utf8_names_are_reversible(self):
        row = self.member(genome='g\u00e9', contig='c\u4e00')
        prefix = representative_allele_id(row).rsplit(':', 2)[0]
        self.assertEqual(prefix, 'g%C3%A9_c%E4%B8%80')
        self.assertEqual([unquote(part) for part in prefix.split('_')], ['g\u00e9', 'c\u4e00'])

    def test_zero_based_half_open_interval(self):
        self.assertEqual(representative_allele_id(self.member(start='0', end='1')), 'HG002_hap1_chr1:0:1')

    def test_invalid_coordinates_fail(self):
        for start, end in [(-1, 1), (2, 2), (3, 2), ('1.0', 2), (True, 2), ('\u0661', 2)]:
            with self.subTest(start=start, end=end), self.assertRaises(ValueError):
                representative_allele_id(self.member(start=start, end=end))

    def test_invalid_names_fail(self):
        for field in ('hal_genome', 'hal_sequence'):
            for value in ('', None, ' x', 'x ', 'x\x00y'):
                row = self.member()
                row[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    representative_allele_id(row)

    def test_mapping_preserves_order_and_input(self):
        members = {'f1': self.member(), 'f2': self.member(genome='A')}
        catalogue = [dict(allele_id='z', allele_rep_fid='f1'), dict(allele_id='a', allele_rep_fid='f2')]
        mapping = allele_id_map(catalogue, members)
        self.assertEqual(list(mapping), ['z', 'a'])
        self.assertEqual(list(mapping.values()), ['HG002_hap1_chr1:950:1350', 'A_chr1:950:1350'])
        self.assertEqual(catalogue[0]['allele_id'], 'z')
        self.assertEqual(members['f1']['hal_genome'], 'HG002_hap1')

    def test_underscore_collision_fails(self):
        members = {'f1': self.member('A_B', 'C'), 'f2': self.member('A', 'B_C')}
        catalogue = [dict(allele_id='a', allele_rep_fid='f1'), dict(allele_id='b', allele_rep_fid='f2')]
        with self.assertRaisesRegex(ValueError, 'collision'):
            allele_id_map(catalogue, members)

    def test_duplicate_internal_id_fails(self):
        row = dict(allele_id='a', allele_rep_fid='f1')
        with self.assertRaisesRegex(ValueError, 'duplicate internal'):
            allele_id_map([row, row], {'f1': self.member()})

    def test_duplicate_representative_fails(self):
        catalogue = [dict(allele_id='a', allele_rep_fid='f1'), dict(allele_id='b', allele_rep_fid='f1')]
        with self.assertRaisesRegex(ValueError, 'collision'):
            allele_id_map(catalogue, {'f1': self.member()})

    def test_unknown_representative_fails(self):
        with self.assertRaisesRegex(ValueError, 'Unknown allele representative'):
            allele_id_map([dict(allele_id='a', allele_rep_fid='absent')], {})


if __name__ == '__main__':
    unittest.main()
