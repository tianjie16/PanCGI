import copy
import json
import math

from .directions import choose, placements, placements_sorted
from .disk_candidates import CandidateCache

CONFIG = dict(windows=[100, 1000, 10000], min_bp=50, min_fraction=0.05, margin=0.02,
              gap_penalty=0, repeat_power=0, span_tolerance=10000, placement_tolerance=1000,
              source_context_bp=1000000)


def combine_hard(by_scale, config):
    windows, tolerance = config['windows'], config['placement_tolerance']
    combined = []
    for large in by_scale[windows[-1]]:
        seed = large['best']
        scores = []
        for window in windows:
            matching = [g['best']['score'] for g in by_scale[window]
                        if g['best']['contig'] == seed['contig'] and g['best']['direction'] == seed['direction']
                        and abs(g['best']['lo'] - seed['lo']) <= tolerance
                        and abs(g['best']['hi'] - seed['hi']) <= tolerance]
            scores.append(max(matching, default=0))
        score = math.prod(scores) ** (1 / len(scores)) if all(s > 0 for s in scores) else 0
        if score > 0:
            combined.append(dict(seed, score=score))
    return placements(combined, tolerance)


def evaluate(record, scratch=None):
    config = CONFIG
    ids = [m['model_fid'] for m in record['candidate_evidence']]
    if not ids or not all(isinstance(n, str) and n for n in ids) or len(ids) != len(set(ids)):
        raise ValueError('Source model IDs must be nonempty and unique')
    if not all(m.get('exact_complete', False) for m in record['candidate_evidence']):
        raise RuntimeError('Incomplete evidence is a technical failure, not NA')
    by_source = {name: {} for name in sorted(ids)}
    params = (config['gap_penalty'], config['repeat_power'], config['min_bp'], config['min_fraction'], config['span_tolerance'])
    with CandidateCache(scratch) as cache:
        for window in config['windows']:
            store = cache.get(record, window, params)
            store.conn.execute('CREATE INDEX IF NOT EXISTS by_source ON candidates(model,score DESC,contig,direction,lo,hi,seq)')
            names = {r[0] for r in store.conn.execute('SELECT DISTINCT model FROM candidates')}
            if not names <= set(ids):
                raise RuntimeError('Unknown candidate source')
            seen = 0
            for name in sorted(ids):
                def records():
                    nonlocal seen
                    query = 'SELECT payload FROM candidates WHERE model=? ORDER BY score DESC,contig,direction,lo,hi,seq'
                    for (payload,) in store.conn.execute(query, (name,)):
                        seen += 1
                        yield json.loads(payload)
                by_source[name][window] = placements_sorted(records(), config['placement_tolerance'])
            if seen != store.unique_n:
                raise RuntimeError('Source grouping omitted or duplicated candidates')
    decisions = {name: choose(combine_hard(scales, config), config['margin']) for name, scales in by_source.items()}
    certified = {name: copy.deepcopy(d) for name, d in decisions.items() if d['call'] == '0'}
    return dict(call='0' if certified else 'NA', reason='at_least_one_certified_source' if certified else 'no_certified_source',
                source_model_n=len(decisions), certified_source_n=len(certified), supporting_sources=list(certified),
                certified_witnesses=certified, source_decisions=decisions)


def propagate(observed_locus, observed_alleles, allele_ids, evidence_call):
    if evidence_call not in ('0', 'NA') or not observed_alleles <= set(allele_ids):
        raise ValueError('Invalid genotype state or allele membership')
    if bool(observed_alleles) != bool(observed_locus):
        raise ValueError('Locus and allele observed membership disagree')
    locus = '1' if observed_locus else evidence_call
    alleles = {a: '1' if a in observed_alleles else ('NA' if locus == 'NA' else '0') for a in allele_ids}
    return locus, alleles
