from __future__ import annotations

import bisect
import csv
import gzip
import json
import os
from itertools import accumulate
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from pancgi_contract import parse_sv_length


def open_text(path: str, mode: str = 'rt'):
    if str(path).endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def normalize_sv_contig_name(contig: object) -> str:
    return str(contig or '').strip()


def _clean_text(val: object) -> str:
    s = str(val or '').strip()
    if s in ('', '.', 'nan', 'NaN', 'None', 'null'):
        return ''
    return s


def _parse_int(val: object, default: object = '') -> object:
    try:
        s = str(val).strip()
        if s in ('', '.', 'nan', 'NaN', 'None', 'null'):
            return default
        return int(float(s))
    except Exception:
        return default


def _effective_interval0(start0: int, end0: int) -> Tuple[int, int]:
    s = int(start0)
    e = int(end0)
    if s < 0 or e < s:
        raise ValueError('Invalid assembly SV interval')
    return s, e


def _overlap_len0(a_start0: int, a_end0: int, b_start0: int, b_end0: int) -> int:
    return max(0, min(int(a_end0), int(b_end0)) - max(int(a_start0), int(b_start0)))


def build_sv_index(path: Optional[str], *, contig_coordinate_base: int = 0) -> Dict[str, Dict[str, object]]:
    path = str(path or '').strip()
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(path)
    if contig_coordinate_base != 0:
        raise ValueError('Internal SV intervals must be 0-based half-open')
    by_contig: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        required = {'ID', 'contig', 'start', 'end', 'SVTYPE', '#CHROM', 'POS', 'VCF_END', 'SVLEN'}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f'SV table is missing required columns: {missing}')
        for line_number, row in enumerate(reader, start=2):
            contig = normalize_sv_contig_name(row.get('contig'))
            if not contig:
                raise ValueError(f'SV table line {line_number} has an empty contig')
            start0 = _parse_int(row.get('start'), None)
            end0 = _parse_int(row.get('end'), None)
            if start0 is None:
                raise ValueError(f'SV table line {line_number} has an invalid start')
            if end0 is None:
                raise ValueError(f'SV table line {line_number} has an invalid end')


            svtype = _clean_text(row.get('SVTYPE')).upper()
            if svtype not in {'INS', 'DEL'}:
                raise ValueError(f'SV table line {line_number} requires INS or DEL')
            eff_start0, eff_end0 = _effective_interval0(int(start0), int(end0))
            if svtype == 'DEL':
                if start0 != end0:
                    raise ValueError('DEL requires a zero-length assembly junction')
                if int(row['POS']) >= int(row['VCF_END']):
                    raise ValueError('DEL requires VCF_END > POS')
            else:
                if start0 >= end0 or int(row['POS']) != int(row['VCF_END']):
                    raise ValueError('INS requires a positive assembly interval and one reference anchor')
            svlen = parse_sv_length(row['SVLEN'], svtype, f'{path}:{line_number}')
            svlen_abs = abs(svlen)
            rec = {
                'id': row['ID'],
                'varid': row.get('VARID', ''),
                'contig': contig,
                'asm_start0': int(start0),
                'asm_end0': int(end0),
                'eff_start0': eff_start0,
                'eff_end0': eff_end0,
                'span_bp': eff_end0 - eff_start0,
                'chrom': _clean_text(row['#CHROM']),
                'pos1': _parse_int(row.get('POS'), ''),
                'vcf_end1': _parse_int(row.get('VCF_END'), ''),
                'svtype': svtype,
                'svlen': svlen,
                'svlen_abs': svlen_abs,
                'gt': _clean_text(row.get('GT')),
                'assembly': _clean_text(row.get('assembly')),
                'TR': _clean_text(row.get('TR')),
                'CONFORMATION': _clean_text(row.get('CONFORMATION')),
                'SD': _clean_text(row.get('SD')),
                'ITYPE_N': _clean_text(row.get('ITYPE_N')),
                'DTYPE_N': _clean_text(row.get('DTYPE_N')),
                'FAM_N': _clean_text(row.get('FAM_N')),
            }
            by_contig[contig].append(rec)
    out: Dict[str, Dict[str, object]] = {}
    for contig, rows in by_contig.items():
        rows.sort(key=lambda r: (int(r['eff_start0']), int(r['eff_end0']), str(r.get('id') or ''), str(r.get('svtype') or '')))
        out[contig] = {
            'rows': rows,
            'starts': [int(r['eff_start0']) for r in rows],
            'prefix_max_end': list(accumulate((int(r['eff_end0']) for r in rows), max)),
        }
    return out


def _classify_overlap(feat_start0: int, feat_end0: int, sv: Dict[str, object]) -> str:
    svtype = str(sv.get('svtype') or '').lower()
    fs = int(feat_start0)
    fe = int(feat_end0)
    ss = int(sv['eff_start0'])
    se = int(sv['eff_end0'])
    if svtype == 'del':
        if not fs < ss < fe:
            raise ValueError('DEL junction must be strictly inside CGI')
        return 'contains_del'
    if fs >= ss and fe <= se:
        return f'inside_{svtype}'
    if ss >= fs and se <= fe:
        return f'contains_{svtype}'
    return f'partial_{svtype}'


def _fold_list(values: List[object], *, fmt: Optional[str] = None) -> str:
    out: List[str] = []
    for v in values:
        if v in ('', None):
            out.append('')
        elif fmt is not None:
            out.append(format(v, fmt))
        else:
            out.append(str(v))
    return ';'.join(out)


def _render_site(chrom: str, pos1: object) -> str:
    if not chrom or pos1 in ('', None):
        return ''
    return f'{chrom}:{pos1}'


def _render_interval(chrom: str, start1: object, end1: object) -> str:
    if not chrom or start1 in ('', None) or end1 in ('', None):
        return ''
    return f'{chrom}:{start1}-{end1}'


def _render_asm_interval(contig: str, start0: int, end0: int) -> str:
    return f'{contig}:{start0}-{end0}'


def fold_overlap_details(details: List[Dict[str, object]]) -> Dict[str, object]:
    details = list(details)
    if not details:
        return {
            'sv_overlap_n': 0,
            'sv_overlap_ids': '',
            'sv_overlap_types': '',
            'sv_overlap_classes': '',
            'sv_overlap_primary_sites': '',
            'sv_overlap_primary_intervals': '',
            'sv_overlap_asm_intervals': '',
            'sv_overlap_bp': '',
            'sv_overlap_pct_cpgi': '',
            'sv_overlap_pct_sv': '',
            'sv_overlap_TR': '',
            'sv_overlap_CONFORMATION': '',
            'sv_overlap_SD': '',
            'sv_overlap_ITYPE_N': '',
            'sv_overlap_DTYPE_N': '',
            'sv_overlap_FAM_N': '',
            'sv_overlap_detail_json': '[]',
        }
    return {
        'sv_overlap_n': len(details),
        'sv_overlap_ids': _fold_list([d.get('id', '') for d in details]),
        'sv_overlap_types': _fold_list([d.get('svtype', '') for d in details]),
        'sv_overlap_classes': _fold_list([d.get('class', '') for d in details]),
        'sv_overlap_primary_sites': _fold_list([_render_site(str(d.get('chrom') or ''), d.get('pos1')) for d in details]),
        'sv_overlap_primary_intervals': _fold_list([_render_interval(str(d.get('chrom') or ''), d.get('pos1'), d.get('vcf_end1')) for d in details]),
        'sv_overlap_asm_intervals': _fold_list([_render_asm_interval(str(d.get('contig') or ''), int(d.get('asm_start0') or 0), int(d.get('asm_end0') or 0)) for d in details]),
        'sv_overlap_bp': _fold_list([int(d.get('overlap_bp') or 0) for d in details]),
        'sv_overlap_pct_cpgi': _fold_list([float(d.get('overlap_pct_cpgi') or 0.0) for d in details], fmt='.6f'),
        'sv_overlap_pct_sv': ';'.join('NA' if d.get('overlap_pct_sv') is None else format(float(d['overlap_pct_sv']), '.6f') for d in details),
        'sv_overlap_TR': _fold_list([d.get('TR', '') for d in details]),
        'sv_overlap_CONFORMATION': _fold_list([d.get('CONFORMATION', '') for d in details]),
        'sv_overlap_SD': _fold_list([d.get('SD', '') for d in details]),
        'sv_overlap_ITYPE_N': _fold_list([d.get('ITYPE_N', '') for d in details]),
        'sv_overlap_DTYPE_N': _fold_list([d.get('DTYPE_N', '') for d in details]),
        'sv_overlap_FAM_N': _fold_list([d.get('FAM_N', '') for d in details]),
        'sv_overlap_detail_json': json.dumps(details, separators=(',', ':'), ensure_ascii=False),
    }


def annotate_feature_overlaps(contig: str, start0: int, end0: int, sv_index: Dict[str, Dict[str, object]]) -> Dict[str, object]:
    base_rec = fold_overlap_details([])
    base_rec.update({
        'sv_ins_overlap_class': 'outside_ins',
        'sv_ins_n': 0,
        'sv_ins_ids': '',
        'sv_ins_primary_sites': '',
        'sv_ins_primary_intervals': '',
        'sv_ins_asm_intervals': '',
        'sv_ins_source': '',
        'sv_ins_longest_id': '',
        'sv_ins_longest_chrom': '',
        'sv_ins_longest_site1': '',
        'sv_ins_longest_start1': '',
        'sv_ins_longest_end1': '',
        'sv_ins_longest_len': '',
        'sv_ins_overlap_bp': '',
        'sv_ins_overlap_pct_cpgi': '',
        'sv_ins_overlap_pct_sv': '',
        'sv_ins_detail_json': '[]',
    })
    key = normalize_sv_contig_name(contig)
    bundle = sv_index.get(key)
    if not bundle:
        return base_rec
    rows = bundle['rows']
    starts = bundle['starts']
    q_start0, q_end0 = _effective_interval0(int(start0), int(end0))
    i = bisect.bisect_left(starts, int(q_end0))
    overlaps: List[Dict[str, object]] = []
    j = i - 1
    while j >= 0:
        r = rows[j]
        if bundle['prefix_max_end'][j] <= int(q_start0):
            break
        if r['svtype'] == 'DEL':
            hit = q_start0 < r['asm_start0'] < q_end0
        else:
            hit = int(r['eff_start0']) < int(q_end0) and int(r['eff_end0']) > int(q_start0)
        if hit:
            overlaps.append(r)
        j -= 1
    if not overlaps:
        return base_rec
    uniq: Dict[Tuple[str, str, int, int], Dict[str, object]] = {}
    for r in overlaps:
        key2 = (str(r.get('id') or ''), str(r.get('svtype') or ''), int(r['eff_start0']), int(r['eff_end0']))
        uniq[key2] = r
    details: List[Dict[str, object]] = []
    cpgi_len = max(1, int(end0) - int(start0))
    for r in uniq.values():
        ov_bp = _overlap_len0(q_start0, q_end0, int(r['eff_start0']), int(r['eff_end0']))
        span_bp = int(r['span_bp'])
        detail = {
            'id': str(r.get('id') or ''),
            'varid': str(r.get('varid') or ''),
            'svtype': str(r.get('svtype') or ''),
            'class': _classify_overlap(int(start0), int(end0), r),
            'contig': key,
            'asm_start0': int(r['asm_start0']),
            'asm_end0': int(r['asm_end0']),
            'eff_start0': int(r['eff_start0']),
            'eff_end0': int(r['eff_end0']),
            'span_bp': span_bp,
            'chrom': str(r.get('chrom') or ''),
            'pos1': r.get('pos1', ''),
            'vcf_end1': r.get('vcf_end1', ''),
            'svlen': r['svlen'],
            'svlen_abs': r['svlen_abs'],
            'gt': str(r.get('gt') or ''),
            'assembly': str(r.get('assembly') or ''),
            'TR': str(r.get('TR') or ''),
            'CONFORMATION': str(r.get('CONFORMATION') or ''),
            'SD': str(r.get('SD') or ''),
            'ITYPE_N': str(r.get('ITYPE_N') or ''),
            'DTYPE_N': str(r.get('DTYPE_N') or ''),
            'FAM_N': str(r.get('FAM_N') or ''),
            'overlap_bp': ov_bp,
            'overlap_pct_cpgi': float(ov_bp) / float(cpgi_len),
            'overlap_pct_sv': float(ov_bp) / float(span_bp) if span_bp else None,
        }
        details.append(detail)
    details.sort(key=lambda d: (int(d.get('overlap_bp') or 0), int(d.get('svlen_abs') or 0), int(d.get('span_bp') or 0), str(d.get('id') or '')), reverse=True)
    out = fold_overlap_details(details)
    ins_details = [d for d in details if str(d.get('svtype') or '').upper() == 'INS']
    if not ins_details:
        out.update({
            'sv_ins_overlap_class': 'outside_ins',
            'sv_ins_n': 0,
            'sv_ins_ids': '',
            'sv_ins_primary_sites': '',
            'sv_ins_primary_intervals': '',
            'sv_ins_asm_intervals': '',
            'sv_ins_source': '',
            'sv_ins_longest_id': '',
            'sv_ins_longest_chrom': '',
            'sv_ins_longest_site1': '',
            'sv_ins_longest_start1': '',
            'sv_ins_longest_end1': '',
            'sv_ins_longest_len': '',
            'sv_ins_overlap_bp': '',
            'sv_ins_overlap_pct_cpgi': '',
            'sv_ins_overlap_pct_sv': '',
            'sv_ins_detail_json': '[]',
        })
        return out
    inside_any = any(str(d.get('class')) == 'inside_ins' for d in ins_details)
    partial_any = any(str(d.get('class')) == 'partial_ins' for d in ins_details)
    contains_any = any(str(d.get('class')) == 'contains_ins' for d in ins_details)
    if inside_any:
        ins_class = 'inside_ins'
    elif partial_any:
        ins_class = 'partial_ins'
    elif contains_any:
        ins_class = 'contains_ins'
    else:
        ins_class = 'partial_ins'
    longest = max(ins_details, key=lambda d: (int(d.get('svlen_abs') or 0), int(d.get('overlap_bp') or 0), int(d.get('span_bp') or 0), str(d.get('id') or '')))
    out.update({
        'sv_ins_overlap_class': ins_class,
        'sv_ins_n': len(ins_details),
        'sv_ins_ids': _fold_list([d.get('id', '') for d in ins_details]),
        'sv_ins_primary_sites': _fold_list([_render_site(str(d.get('chrom') or ''), d.get('pos1')) for d in ins_details]),
        'sv_ins_primary_intervals': _fold_list([_render_interval(str(d.get('chrom') or ''), d.get('pos1'), d.get('vcf_end1')) for d in ins_details]),
        'sv_ins_asm_intervals': _fold_list([_render_asm_interval(str(d.get('contig') or ''), int(d.get('asm_start0') or 0), int(d.get('asm_end0') or 0)) for d in ins_details]),
        'sv_ins_source': 'provided_sv_callset',
        'sv_ins_longest_id': str(longest.get('id') or ''),
        'sv_ins_longest_chrom': str(longest.get('chrom') or ''),
        'sv_ins_longest_site1': longest.get('pos1', ''),
        'sv_ins_longest_start1': longest.get('pos1', ''),
        'sv_ins_longest_end1': longest.get('vcf_end1', ''),
        'sv_ins_longest_len': longest.get('svlen_abs', ''),
        'sv_ins_overlap_bp': _fold_list([int(d.get('overlap_bp') or 0) for d in ins_details]),
        'sv_ins_overlap_pct_cpgi': _fold_list([float(d.get('overlap_pct_cpgi') or 0.0) for d in ins_details], fmt='.6f'),
        'sv_ins_overlap_pct_sv': _fold_list([float(d.get('overlap_pct_sv') or 0.0) for d in ins_details], fmt='.6f'),
        'sv_ins_detail_json': json.dumps(ins_details, separators=(',', ':'), ensure_ascii=False),
    })
    return out
