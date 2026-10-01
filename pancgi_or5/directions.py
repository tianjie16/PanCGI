from dataclasses import dataclass
from collections import defaultdict
import math

import numpy as np


@dataclass(frozen=True)
class Anchor:
    contig: str
    direction: int
    s0: float
    s1: float
    t0: float
    t1: float
    multiplicity: int


def anchors(runs, side, window):
    result = []
    seen = set()
    low, high = (-window, 0) if side == 'left' else (0, window)
    for r in runs:
        s0, s1 = max(low, r['source_start_bp']), min(high, r['source_end_bp'])
        if s1 <= s0:
            continue
        direction = int(r['direction'])
        start, end = float(r['target_start1'] - 1), float(r['target_end1'])
        t0, t1 = (start, end) if direction == 1 else (-end, -start)
        if not math.isclose(t1 - t0, r['source_end_bp'] - r['source_start_bp'], abs_tol=1e-6):
            raise RuntimeError('an exact run has inconsistent source and target bp spans')
        t0 += s0 - r['source_start_bp']
        t1 -= r['source_end_bp'] - s1
        a = Anchor(str(r['chrom']), direction, s0, s1, t0, t1, int(r['placement_n']))
        if a.multiplicity < 1:
            raise RuntimeError('invalid exact-run multiplicity')
        if a not in seen:
            result.append(a)
            seen.add(a)
    return sorted(result, key=lambda a: (a.contig, a.direction, a.s0, a.s1, a.t0, a.t1))


def chain_scores(items, max_gap, gap_penalty, repeat_power=0, reverse=False, vectorized=True):
    if not items:
        return [], []
    transformed = [(i, -a.s1, -a.s0, -a.t1, -a.t0) if reverse else (i, a.s0, a.s1, a.t0, a.t1)
                   for i, a in enumerate(items)]
    transformed.sort(key=lambda x: (x[1], x[2], x[3], x[4], x[0]))
    arr = np.asarray([x[1:] for x in transformed], dtype=float)
    scores = np.zeros(len(items))
    support = np.zeros(len(items))
    for j, (original, s0, s1, t0, t1) in enumerate(transformed):
        weight = (s1 - s0) / items[original].multiplicity ** repeat_power
        scores[j] = weight
        support[j] = s1 - s0
        if not j:
            continue
        if vectorized:
            sg, tg = s0 - arr[:j, 1], t0 - arr[:j, 3]
            valid = (sg >= 0) & (tg >= 0) & (sg <= max_gap) & (tg <= max_gap)
            candidates = np.where(valid, scores[:j] - gap_penalty * np.abs(tg - sg), -np.inf)
            k = int(np.argmax(candidates))
            best = candidates[k]
        else:
            best, k = -np.inf, -1
            for i in range(j):
                sg, tg = s0 - arr[i, 1], t0 - arr[i, 3]
                if 0 <= sg <= max_gap and 0 <= tg <= max_gap:
                    value = scores[i] - gap_penalty * abs(tg - sg)
                    if value > best:
                        best, k = value, i
        if best > 0:
            scores[j] += best
            support[j] += support[k]
    restored_scores, restored_support = [0.0] * len(items), [0.0] * len(items)
    for j, (original, *_rest) in enumerate(transformed):
        restored_scores[original] = float(scores[j])
        restored_support[original] = float(support[j])
    return restored_scores, restored_support


def iter_pair_candidates(model, window, gap_penalty, repeat_power, min_bp, min_fraction, span_tolerance):
    if not model.get('exact_complete', False):
        raise RuntimeError('incomplete evidence cannot produce a genotype')
    left, right = anchors(model['left_runs'], 'left', window), anchors(model['right_runs'], 'right', window)
    groups = defaultdict(lambda: [[], []])
    for side, records in enumerate((left, right)):
        for a in records:
            groups[(a.contig, a.direction)][side].append(a)
    left_den = min(window, float(model['left_flank_len']))
    right_den = min(window, float(model['right_flank_len']))
    if left_den <= 0 or right_den <= 0:
        return
    env = float(model['envelope_len'])
    for (contig, direction), (ls, rs) in groups.items():
        if not ls or not rs:
            continue
        lscore, lbp = chain_scores(ls, window, gap_penalty, repeat_power)
        rscore, rbp = chain_scores(rs, window, gap_penalty, repeat_power, reverse=True)
        for i, l in enumerate(ls):
            if lbp[i] < min_bp or lscore[i] / left_den < min_fraction:
                continue
            for j, r in enumerate(rs):
                if rbp[j] < min_bp or rscore[j] / right_den < min_fraction:
                    continue
                target_gap = r.t0 - l.t1
                source_gap = env + r.s0 - l.s1
                if target_gap < 0 or abs(target_gap - source_gap) > span_tolerance:
                    continue
                left_est = l.t1 - l.s1
                right_est = r.t0 - r.s0
                if right_est < left_est:
                    continue
                a, b = lscore[i] / left_den, rscore[j] / right_den
                score = 2 * a * b / (a + b)
                lo, hi = (left_est, right_est) if direction == 1 else (-right_est, -left_est)
                yield {'contig': contig, 'direction': direction, 'lo': lo, 'hi': hi,
                                   'score': score, 'left_fraction': a, 'right_fraction': b,
                                   'left_support_bp': lbp[i], 'right_support_bp': rbp[j],
                                   'model': model['model_fid'], 'window': window}


def placements(candidates, tolerance):
    return placements_sorted(sorted(candidates, key=lambda x: (-x['score'], x['contig'], x['direction'], x['lo'], x['hi'], x['model'])), tolerance)


def placements_sorted(candidates, tolerance):
    if tolerance <= 0:
        raise ValueError('positive placement tolerance required')
    groups = []
    buckets = defaultdict(list)
    for c in candidates:
        bin0, bin1 = math.floor(c['lo'] / tolerance), math.floor(c['hi'] / tolerance)
        eligible = []
        for i in range(bin0 - 1, bin0 + 2):
            for j in range(bin1 - 1, bin1 + 2):
                eligible.extend(buckets.get((c['contig'], c['direction'], i, j), []))
        found = None
        for index in sorted(eligible):
            g = groups[index]
            if (c['contig'], c['direction']) != (g['best']['contig'], g['best']['direction']):
                continue
            if max(g['lo_max'], c['lo']) - min(g['lo_min'], c['lo']) <= tolerance and max(g['hi_max'], c['hi']) - min(g['hi_min'], c['hi']) <= tolerance:
                found = g
                break
        if found is None:
            buckets[(c['contig'], c['direction'], bin0, bin1)].append(len(groups))
            groups.append({'best': dict(c), 'lo_min': c['lo'], 'lo_max': c['lo'],
                           'hi_min': c['hi'], 'hi_max': c['hi'], 'models': {c['model']}, 'windows': {c['window']}})
        else:
            found['lo_min'], found['lo_max'] = min(found['lo_min'], c['lo']), max(found['lo_max'], c['lo'])
            found['hi_min'], found['hi_max'] = min(found['hi_min'], c['hi']), max(found['hi_max'], c['hi'])
            found['models'].add(c['model'])
            found['windows'].add(c['window'])
    return groups


def choose(groups, margin):
    if not groups:
        return {'call': 'NA', 'reason': 'no_supported_placement', 'best_score': 0, 'runner_score': 0, 'placements': 0}
    groups = sorted(groups, key=lambda g: -g['best']['score'])
    best = groups[0]['best']
    runner = groups[1]['best']['score'] if len(groups) > 1 else 0
    accepted = best['score'] > 0 and (best['score'] - runner) / best['score'] >= margin and best['score'] > runner
    return {'call': '0' if accepted else 'NA',
            'reason': 'dominant_ordered_placement' if accepted else 'competing_placements',
            'best_score': best['score'], 'runner_score': runner, 'placements': len(groups),
            'contig': best['contig'], 'direction': best['direction'], 'lo': best['lo'], 'hi': best['hi'],
            'model': best['model']}


