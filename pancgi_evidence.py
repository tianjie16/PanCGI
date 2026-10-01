import math


def object_schema(properties):
    return dict(type='object', properties=properties, required=list(properties), additionalProperties=False)


TEXT = dict(type='string', minLength=1)
NUMBER = dict(type='number')
COUNT = dict(type='integer', minimum=0)
SOURCE_MEMBER = object_schema(dict(type=dict(const='member'), member_id=TEXT, hal_genome=TEXT,
                                   hal_sequence=TEXT, start0=COUNT, end0=COUNT))
SOURCE_CONTEXT = object_schema(dict(type=dict(const='primary_reference_context'), locus_id=TEXT,
                                    hal_genome=TEXT, hal_sequence=TEXT, start0=COUNT, end0=COUNT))
SOURCE = dict(oneOf=[SOURCE_MEMBER, SOURCE_CONTEXT])
PLACEMENT = object_schema(dict(hal_genome=TEXT, hal_sequence=TEXT, direction=dict(enum=[-1, 1]),
    estimated_start0=NUMBER, estimated_end0=NUMBER))
DECISION = object_schema(dict(source=SOURCE, callability=dict(enum=['callable', 'unresolved']),
    reason=dict(enum=['dominant_ordered_placement', 'competing_placements', 'no_supported_placement']),
    best_score=dict(type='number', minimum=0), runner_score=dict(type='number', minimum=0),
    placements=COUNT, placement=dict(oneOf=[dict(type='null'), PLACEMENT])))
OBSERVED = object_schema(dict(callability=dict(const='not_evaluated_observed'), reason=dict(const='observed_member')))
EVALUATED = object_schema(dict(callability=dict(enum=['callable', 'unresolved']),
    reason=dict(enum=['at_least_one_certified_source', 'no_certified_source']), source_model_n=COUNT,
    certified_source_n=COUNT, supporting_sources=dict(type='array', items=SOURCE, uniqueItems=True),
    source_decisions=dict(type='array', items=DECISION, minItems=1)))
EVIDENCE_SCHEMA = dict(**object_schema(dict(locus_id=TEXT, hal_genome=TEXT,
    genotype=dict(enum=['0', '1', 'NA']), evidence=dict(oneOf=[OBSERVED, EVALUATED]))),
    **{'$schema': 'https://json-schema.org/draft/2020-12/schema'})


def require(condition, message):
    if not condition:
        raise ValueError('Genotype evidence: ' + message)


def keys(value, expected):
    require(isinstance(value, dict) and set(value) == set(expected), 'unexpected fields')


def source_key(source):
    return (source['type'], source.get('member_id') if source['type'] == 'member' else source.get('locus_id'))


def export_record(record, members, member_assignments, loci, identities, primary):
    lid, genome, genotype = record['locus_id'], record['hal_genome'], record['genotype']
    evidence = record['evidence']
    if evidence['reason'] == 'observed_member':
        keys(evidence, ('call', 'reason'))
        require(genotype == '1' and evidence['call'] == '0', 'invalid observed record')
        return dict(locus_id=lid, hal_genome=genome, genotype=genotype,
                    evidence=dict(callability='not_evaluated_observed', reason='observed_member'))

    def source(identifier):
        if identifier in members:
            member = members[identifier]
            require(identifier in member_assignments and member_assignments[identifier][0] == lid,
                    'source member belongs to a different locus')
            return dict(type='member', member_id=member['member_id'], hal_genome=member['hal_genome'],
                        hal_sequence=member['hal_sequence'], start0=int(member['start0']), end0=int(member['end0']))
        require(identifier == 'main_reference_context:' + lid, 'unknown source identity')
        locus = loci[lid]
        require(locus['locus_type'] == 'Novel' and locus['primary_sequence'] and
                locus['primary_start0'] != '' and locus['primary_end0'] != '', 'missing reference context')
        return dict(type='primary_reference_context', locus_id=lid, hal_genome=primary,
                    hal_sequence=locus['primary_sequence'], start0=int(locus['primary_start0']), end0=int(locus['primary_end0']))

    keys(evidence, ('call', 'reason', 'source_model_n', 'certified_source_n', 'supporting_sources',
                    'certified_witnesses', 'source_decisions'))
    decisions = []
    certified = {name: decision for name, decision in evidence['source_decisions'].items() if decision['call'] == '0'}
    require(evidence['certified_witnesses'] == certified and evidence['supporting_sources'] == list(certified),
            'inconsistent internal certified sources')
    for identifier, decision in evidence['source_decisions'].items():
        placement = None
        expected = ['call', 'reason', 'best_score', 'runner_score', 'placements']
        if 'contig' in decision:
            expected += ['contig', 'direction', 'lo', 'hi', 'model']
            require(decision['model'] == identifier, 'source/placement identity mismatch')
            identity = identities[decision['contig']]
            require(identity['hal_genome'] == genome, 'placement belongs to a different target genome')
            placement = dict(hal_genome=genome, hal_sequence=identity['hal_sequence'], direction=decision['direction'],
                             estimated_start0=decision['lo'], estimated_end0=decision['hi'])
        keys(decision, expected)
        require(decision['call'] in ('0', 'NA'), 'invalid source call')
        decisions.append(dict(source=source(identifier), callability='callable' if decision['call'] == '0' else 'unresolved',
                              reason=decision['reason'], best_score=decision['best_score'], runner_score=decision['runner_score'],
                              placements=decision['placements'], placement=placement))
    require(evidence['call'] in ('0', 'NA') and genotype == evidence['call'], 'invalid unobserved call')
    return dict(locus_id=lid, hal_genome=genome, genotype=genotype, evidence=dict(
        callability='callable' if evidence['call'] == '0' else 'unresolved', reason=evidence['reason'],
        source_model_n=evidence['source_model_n'], certified_source_n=evidence['certified_source_n'],
        supporting_sources=[source(identifier) for identifier in evidence['supporting_sources']], source_decisions=decisions))


def validate_record(record, loci, alleles, members, genomes, sequences, margin, locus_alleles):
    keys(record, ('locus_id', 'hal_genome', 'genotype', 'evidence'))
    lid, genome, genotype = record['locus_id'], record['hal_genome'], record['genotype']
    require(lid in loci and genome in genomes and genotype in ('0', '1', 'NA'), 'unknown cell')
    evidence = record['evidence']
    require(isinstance(evidence, dict), 'invalid evidence object')
    if genotype == '1':
        require(evidence == dict(callability='not_evaluated_observed', reason='observed_member'), 'observed callability mismatch')
        return
    keys(evidence, ('callability', 'reason', 'source_model_n', 'certified_source_n', 'supporting_sources', 'source_decisions'))
    require(type(evidence['source_model_n']) is int and evidence['source_model_n'] > 0, 'invalid source count')
    require(type(evidence['certified_source_n']) is int and evidence['certified_source_n'] >= 0, 'invalid certified count')
    require(isinstance(evidence['source_decisions'], list) and isinstance(evidence['supporting_sources'], list), 'invalid source lists')

    def check_source(source):
        require(isinstance(source, dict), 'invalid source')
        kind = source.get('type')
        if kind == 'member':
            keys(source, SOURCE_MEMBER['required'])
            require(source['member_id'] in members, 'unknown source member')
            member = members[source['member_id']]
            require(member['locus_id'] == lid and member['allele_id'] in alleles, 'source member locus mismatch')
            require(source['member_id'] in (loci[lid]['anchor_member_id'], alleles[member['allele_id']]['representative_member_id']),
                    'source member is neither an anchor nor a representative')
            require(all(source[key] == member[key] for key in ('hal_genome', 'hal_sequence')), 'source member identity mismatch')
            expected = int(member['start0']), int(member['end0'])
        else:
            require(kind == 'primary_reference_context', 'unknown source type')
            keys(source, SOURCE_CONTEXT['required'])
            locus = loci[lid]
            require(source['locus_id'] == lid and locus['locus_type'] == 'Novel', 'reference context locus mismatch')
            require(source['hal_genome'] in genomes and genomes[source['hal_genome']]['role'] == 'primary_reference',
                    'reference context genome mismatch')
            require(source['hal_sequence'] == locus['primary_sequence'], 'reference context sequence mismatch')
            expected = int(locus['primary_start0']), int(locus['primary_end0'])
        require(all(type(source[key]) is int for key in ('start0', 'end0')) and
                (source['start0'], source['end0']) == expected and 0 <= expected[0] < expected[1], 'source coordinates mismatch')
        require((source['hal_genome'], source['hal_sequence']) in sequences, 'undeclared source sequence')

    seen, certified = set(), []
    for decision in evidence['source_decisions']:
        keys(decision, DECISION['required'])
        source = decision['source']
        check_source(source)
        key = source_key(source)
        require(key not in seen, 'duplicate source')
        seen.add(key)
        n, best, runner = decision['placements'], decision['best_score'], decision['runner_score']
        require(type(n) is int and n >= 0, 'invalid placement count')
        require(all(type(value) in (int, float) and math.isfinite(value) and value >= 0 for value in (best, runner)) and runner <= best,
                'invalid placement scores')
        accepted = n > 0 and best > 0 and best > runner and (best - runner) / best >= margin
        expected_reason = 'dominant_ordered_placement' if accepted else ('competing_placements' if n else 'no_supported_placement')
        require(decision['reason'] == expected_reason and decision['callability'] == ('callable' if accepted else 'unresolved'),
                'placement decision contradicts scores')
        if not n:
            require(best == runner == 0 and decision['placement'] is None, 'unsupported placement has coordinates or score')
        else:
            placement = decision['placement']
            keys(placement, PLACEMENT['required'])
            require(placement['hal_genome'] == genome and (genome, placement['hal_sequence']) in sequences,
                    'placement target identity mismatch')
            require(type(placement['direction']) is int and placement['direction'] in (-1, 1), 'invalid direction')
            coordinates = [placement['estimated_start0'], placement['estimated_end0']]
            require(all(type(value) in (int, float) and math.isfinite(value) for value in coordinates) and coordinates[0] <= coordinates[1],
                    'invalid estimated placement coordinates')
        if accepted:
            certified.append(source)
    expected_sources = {('member', alleles[aid]['representative_member_id']) for aid in locus_alleles[lid]}
    if loci[lid]['anchor_member_id']:
        expected_sources.add(('member', loci[lid]['anchor_member_id']))
    else:
        expected_sources.add(('primary_reference_context', lid))
    require(seen == expected_sources, 'missing or extra source models')
    require(evidence['source_model_n'] == len(seen) and evidence['certified_source_n'] == len(certified), 'source counts do not close')
    require(evidence['supporting_sources'] == certified, 'supporting sources disagree with source decisions')
    require(genotype == ('0' if certified else 'NA'), 'certified sources contradict genotype')
    require(evidence['callability'] == ('callable' if certified else 'unresolved') and
            evidence['reason'] == ('at_least_one_certified_source' if certified else 'no_certified_source'), 'cell callability mismatch')
