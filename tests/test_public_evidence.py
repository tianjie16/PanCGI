import copy
import json
import tempfile
import unittest
from pathlib import Path

from pancgi_evidence import export_record, validate_record
from pancgi_public_parameters import STAGES, public_parameters


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.loci = {'1': dict(locus_type='Novel', primary_sequence='ref:chr', primary_start0='10', primary_end0='20', anchor_member_id='')}
        self.genomes = {'ref': dict(role='primary_reference'), 'sample': dict(role='sample')}
        self.members = {'M1': dict(member_id='M1', locus_id='1', allele_id='A1', hal_genome='sample', hal_sequence='sample:chr', start0='30', end0='40')}
        self.alleles = {'A1': dict(locus_id='1', representative_member_id='M1')}
        self.sequences = {('ref', 'ref:chr'), ('sample', 'sample:chr')}
        self.identities = {'target': dict(hal_genome='ref', hal_sequence='ref:chr')}
        decision = dict(call='0', reason='dominant_ordered_placement', best_score=1.0, runner_score=0.0,
                        placements=1, contig='target', direction=-1, lo=10.0, hi=10.0, model='fid')
        no_support = dict(call='NA', reason='no_supported_placement', best_score=0, runner_score=0, placements=0)
        self.internal = dict(locus_id='1', hal_genome='ref', genotype='0', evidence=dict(
            call='0', reason='at_least_one_certified_source', source_model_n=2, certified_source_n=1,
            supporting_sources=['fid'], certified_witnesses={'fid': decision},
            source_decisions={'fid': decision, 'main_reference_context:1': no_support}))

    def export(self, record=None):
        return export_record(record or self.internal, {'fid': self.members['M1']}, {'fid': ('1', 'A1')},
                             self.loci, self.identities, 'ref')

    def check(self, record):
        validate_record(record, self.loci, self.alleles, self.members, self.genomes, self.sequences, 0.02, {'1': ['A1']})

    def test_typed_sources_and_targets(self):
        record = self.export()
        self.check(record)
        self.assertEqual(record['evidence']['source_decisions'][0]['source']['member_id'], 'M1')
        self.assertEqual(record['evidence']['source_decisions'][1]['source']['type'], 'primary_reference_context')
        self.assertEqual(record['evidence']['source_decisions'][0]['placement']['hal_sequence'], 'ref:chr')
        self.assertNotIn('fid', json.dumps(record))
        self.assertNotIn('main_reference_context:', json.dumps(record))

    def test_observed_has_no_placeholder_call(self):
        record = self.export(dict(locus_id='1', hal_genome='sample', genotype='1', evidence=dict(call='0', reason='observed_member')))
        self.check(record)
        self.assertEqual(record['evidence'], dict(callability='not_evaluated_observed', reason='observed_member'))

    def test_unresolved_sources_preserved(self):
        record = copy.deepcopy(self.internal)
        ev = record['evidence']
        decision = ev['source_decisions']['fid']
        decision.update(call='NA', reason='competing_placements', placements=2, runner_score=1.0)
        ev.update(call='NA', reason='no_certified_source', certified_source_n=0, supporting_sources=[], certified_witnesses={})
        record['genotype'] = 'NA'
        self.check(self.export(record))

    def test_unknown_internal_source_rejected(self):
        record = copy.deepcopy(self.internal)
        record['evidence']['source_decisions']['unknown'] = record['evidence']['source_decisions'].pop('main_reference_context:1')
        with self.assertRaises(ValueError):
            self.export(record)

    def test_wrong_target_genome_rejected(self):
        self.identities['target']['hal_genome'] = 'sample'
        with self.assertRaises(ValueError):
            self.export()

    def test_internal_witness_mismatch_rejected(self):
        self.internal['evidence']['supporting_sources'] = []
        with self.assertRaises(ValueError):
            self.export()

    def test_public_evidence_rejects_mutations(self):
        mutations = [
            lambda r: r['evidence'].update(call='0'),
            lambda r: r['evidence'].update(certified_source_n=0),
            lambda r: r['evidence'].update(supporting_sources=[]),
            lambda r: r['evidence'].update(callability='unresolved'),
            lambda r: r.update(genotype='NA'),
            lambda r: r['evidence']['source_decisions'][0]['source'].update(member_id='unknown'),
            lambda r: r['evidence']['source_decisions'][0]['source'].update(start0=31),
            lambda r: r['evidence']['source_decisions'][0]['source'].update(start0=True),
            lambda r: r['evidence']['source_decisions'][0]['source'].update(hal_genome='ref'),
            lambda r: r['evidence']['source_decisions'][1]['source'].update(locus_id='other'),
            lambda r: r['evidence']['source_decisions'][1]['source'].update(hal_genome='sample'),
            lambda r: r['evidence']['source_decisions'][0]['placement'].update(hal_sequence='unknown'),
            lambda r: r['evidence']['source_decisions'][0]['placement'].update(direction=0),
            lambda r: r['evidence']['source_decisions'][0]['placement'].update(estimated_end0=9),
            lambda r: r['evidence']['source_decisions'][0].update(best_score=float('nan')),
            lambda r: r['evidence']['source_decisions'][0].update(runner_score=2),
            lambda r: r['evidence']['source_decisions'][0].update(placements=True),
            lambda r: r['evidence']['source_decisions'][1].update(best_score=1),
        ]
        for i, mutation in enumerate(mutations):
            with self.subTest(case=i):
                record = self.export()
                mutation(record)
                with self.assertRaises(ValueError):
                    self.check(record)

    def test_missing_source_even_with_consistent_counts_rejected(self):
        record = self.export()
        record['evidence']['source_decisions'].pop()
        record['evidence']['source_model_n'] = 1
        with self.assertRaises(ValueError):
            self.check(record)

    def test_duplicate_source_rejected(self):
        record = self.export()
        record['evidence']['source_decisions'].append(record['evidence']['source_decisions'][0])
        with self.assertRaises(ValueError):
            self.check(record)

    def test_nonrepresentative_source_rejected(self):
        self.alleles['A1']['representative_member_id'] = 'another'
        with self.assertRaises(ValueError):
            self.check(self.export())


class ParameterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.values = {cmd: dict(cmd=cmd, **{field: 1 for field in fields.split()}) for cmd, (_, fields) in STAGES.items()}
        self.values['cluster-alleles-prod'].update(
            locus_ids_file='', locus_ids='', locus_start_index=None, locus_end_index=None,
            locus_shard_count=0, locus_shard_index=0, locus_list_only=False,
            very_long_backend='external', very_long_external_template='"/private folder/venv/python" "/application/wfa_longalign_wrapper.py" --seq1 "{seq1}" --seq2 "{seq2}" --log "/work/log"')

    def export(self):
        for cmd, value in self.values.items():
            (self.root/(cmd + '.json')).write_text(json.dumps(value))
        return public_parameters(self.root)

    def test_only_typed_public_values_exported(self):
        self.values['make-cpgi-fasta'].update(hal='/private/source.hal', hal2fasta='C:\\tools\\hal2fasta.exe')
        output = self.export()
        self.assertEqual(output['allele_clustering']['very_long_backend'], 'pywfa')
        self.assertNotIn('/private', json.dumps(output))
        self.assertNotIn('C:', json.dumps(output))
        self.assertEqual(set(output), {stage for stage, _ in STAGES.values()})

    def test_nonnumeric_scientific_value_rejected(self):
        self.values['cluster-alleles-prod']['identity'] = '/private/path'
        with self.assertRaises(ValueError):
            self.export()

    def test_custom_template_rejected(self):
        self.values['cluster-alleles-prod']['very_long_external_template'] = 'other-program --input {seq1}'
        with self.assertRaises(ValueError):
            self.export()

    def test_unknown_parameter_rejected(self):
        self.values['cluster-alleles-prod']['new_threshold'] = 0.9
        with self.assertRaises(ValueError):
            self.export()

    def test_subset_not_described_as_full_run(self):
        self.values['cluster-alleles-prod']['locus_ids_file'] = '/private/subset'
        with self.assertRaises(ValueError):
            self.export()

    def test_missing_stage_rejected(self):
        self.values.pop('unfold-graph')
        with self.assertRaises(ValueError):
            self.export()


if __name__ == '__main__':
    unittest.main()
