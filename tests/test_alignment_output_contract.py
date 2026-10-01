import argparse
import json
import subprocess
import unittest
from unittest.mock import patch

import cpgi_nr_prod as core
import wfa_longalign_wrapper as wfa


class AlignmentOutputContractTests(unittest.TestCase):
    def parse(self, output):
        args = argparse.Namespace(very_long_external_template='align {seq1} {seq2}')
        result = subprocess.CompletedProcess('align', 0, output, '')
        with patch.object(core.subprocess, 'run', return_value=result) as call:
            identity = core.external_long_identity('ACGTACGT', 'ACGTTCGT', args)
        self.assertEqual(call.call_count, 1)
        return identity

    def test_numeric_identity(self):
        for value in (0, 0.125, 1):
            with self.subTest(value=value):
                self.assertEqual(self.parse(json.dumps({'identity': value})), value)

    def test_actual_wrapper_output(self):
        output = wfa.alignment_identity('ACGTACGT', 'ACGTTCGT')
        self.assertEqual(self.parse(json.dumps(output)), output['identity'])

    def test_malformed_output_rejected(self):
        invalid = [
            '', '0.5\t1\t2', '{"matches":1,"alignment_length":2}',
            '[]', 'null', '{"identity":"0.5"}', '{"identity":true}',
            '{"identity":NaN}', '{"identity":Infinity}',
            '{"identity":-0.1}', '{"identity":1.1}',
            '{"identity":0.5}\n{"identity":0.6}',
        ]
        for output in invalid:
            with self.subTest(output=output), self.assertRaises(RuntimeError):
                self.parse(output)

    def test_alignment_failure_is_not_retried(self):
        args = argparse.Namespace(very_long_external_template='align {seq1} {seq2}')
        with patch.object(core.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'align')) as call:
            with self.assertRaises(RuntimeError):
                core.external_long_identity('ACGT', 'ACGT', args)
        self.assertEqual(call.call_count, 1)


if __name__ == '__main__':
    unittest.main()
