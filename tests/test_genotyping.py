import copy
import unittest
from pancgi_genotyping.caller import evaluate


def source(name, offsets):
    return dict(model_fid=name, exact_complete=True, envelope_len=200,
        left_flank_len=1000000, right_flank_len=1000000,
        left_runs=[dict(source_start_bp=-1000000, source_end_bp=0, target_start1=o-1000000+1,
            target_end1=o, direction=1, chrom='arbitrary', placement_n=1) for o in offsets],
        right_runs=[dict(source_start_bp=0, source_end_bp=1000000, target_start1=o+201,
            target_end1=o+1000200, direction=1, chrom='arbitrary', placement_n=1) for o in offsets])


class ORTests(unittest.TestCase):
    def test_reliable_source_not_vetoed(self):
        for other in ([2000000], [2000000, 4000000]):
            result = evaluate(dict(candidate_evidence=[source('a', [1000000]), source('b', other)]))
            self.assertEqual(result['call'], '0')
            self.assertIn('a', result['certified_witnesses'])

    def test_tied_source_is_na(self):
        result = evaluate(dict(candidate_evidence=[source('a', [1000000, 4000000])]))
        self.assertEqual(result['call'], 'NA')

    def test_fixed_evidence_independent_of_labels(self):
        record = dict(candidate_evidence=[source('a', [1000000])])
        self.assertEqual(evaluate(record), evaluate(dict(record, truth='NA', sd_overlap=True)))

    def test_no_one_sided_rescue(self):
        model = source('a', [1000000])
        model['left_runs'] = []
        self.assertEqual(evaluate(dict(candidate_evidence=[model]))['call'], 'NA')


if __name__ == '__main__':
    unittest.main()
