from __future__ import annotations

import argparse
import gc
import csv
import gzip
import glob
import json
import marshal
import math
import bisect
import os
import random
import multiprocessing as mp
import sqlite3
import statistics
import subprocess
import sys
import tempfile
import warnings
from collections import Counter, defaultdict, OrderedDict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import parasail
except Exception:
    parasail = None

import pancgi_core as base
import pancgi_pathbed as pbh
import cpgi_sv_annot as svh
import pancgi_graph_unfold as graph_unfold
import pancgi_hal_runtime as halrt
import pancgi_mapping as mapping


try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


_SV_INS_CACHE: Dict[str, Dict[str, Dict[str, object]]] = {}
_PAIRWISE2 = None


def load_pairwise2():
    global _PAIRWISE2
    if _PAIRWISE2 is not None:
        return _PAIRWISE2
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            from Bio import pairwise2 as module
    except Exception as exc:
        raise RuntimeError(
            'Biopython pairwise2 is required only for center-star MSA output; '
            'install biopython or use --dump-msa-backend none/mafft.'
        ) from exc
    _PAIRWISE2 = module
    return module

def open_text(path: str, mode: str = 'rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def ensure_parent(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def temp_output_path(path: str) -> str:
    if path.endswith('.gz'):
        return f"{path[:-3]}.tmp.{os.getpid()}.gz"
    return f"{path}.tmp.{os.getpid()}"


def finalize_output_path(tmp_path: str, final_path: str) -> None:
    ensure_parent(final_path)
    os.replace(tmp_path, final_path)


def write_optional_parquet(df: pd.DataFrame, path: str, *, index: bool = False) -> None:
    ensure_parent(path)
    df.to_parquet(path, index=index)


def json_compact(x: object) -> str:
    return json.dumps(x, separators=(',', ':'), ensure_ascii=False)


def json_pretty(x: object) -> str:
    return json.dumps(x, indent=2, ensure_ascii=False)


def dominant_value(values: Sequence[str]) -> Optional[str]:
    vals = [v for v in values if v not in (None, '', '.')]
    if not vals:
        return None
    return Counter(vals).most_common(1)[0][0]


def median_int(values: Sequence[int]) -> Optional[int]:
    vals = [int(v) for v in values if v not in (None, '', '.')]
    if not vals:
        return None
    return int(round(statistics.median(vals)))


def reverse_complement(seq: Optional[str]) -> Optional[str]:
    if not seq:
        return seq
    table = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(table)[::-1]


def flip_token(tok: str) -> str:
    if tok.startswith('>'):
        return '<' + tok[1:]
    if tok.startswith('<'):
        return '>' + tok[1:]
    raise ValueError(f'Bad token: {tok}')


def reverse_tokens(tokens: Sequence[str]) -> List[str]:
    return [flip_token(t) for t in reversed(tokens)]


def write_tsv_header(path: str, header: List[str]):
    ensure_parent(path)
    fh = open_text(path, 'wt')
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(header)
    return fh, w


def read_tsv(path: str) -> pd.DataFrame:


    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = list(reader)
        cols = list(reader.fieldnames or [])
    if not cols:
        return pd.DataFrame()
    return pd.DataFrame(rows, columns=cols)


def iter_catalog_labels(path: Optional[str]) -> List[str]:
    if not path or not os.path.exists(path):
        return []
    labels: List[str] = []
    seen = set()
    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            lab = str(row.get('label') or '').strip()
            if lab and lab not in seen:
                labels.append(lab)
                seen.add(lab)
    return labels


INTERNAL_CATALOG_COLUMNS = [
    'label',
    'kind',
    'graph_sample',
    'graph_hap',
    'hal_genome',
    'bed',
    'cpgi_fa',
    'path_bed_dir',
    'sv_tsv',
]


def validate_internal_catalog_header(fieldnames: Optional[Sequence[str]]) -> None:
    observed = list(fieldnames or [])
    if observed != INTERNAL_CATALOG_COLUMNS:
        raise ValueError(f'Invalid PanCGI publication catalog columns: {observed!r}')


def validate_internal_catalog_row(row: Dict[str, str], line_number: int) -> None:
    label = str(row.get('label') or '').strip()
    kind = str(row.get('kind') or '').strip()
    required_nonempty = ['label', 'kind', 'graph_sample', 'graph_hap', 'hal_genome', 'bed', 'cpgi_fa', 'path_bed_dir']
    missing = [column for column in required_nonempty if not str(row.get(column) or '').strip()]
    if kind == 'assembly' and not str(row.get('sv_tsv') or '').strip():
        missing.append('sv_tsv')
    if missing:
        raise ValueError(f'Catalog line {line_number} ({label or "unlabelled"}) has empty required fields: {sorted(set(missing))}')
    if kind not in {'reference_primary', 'reference_comparison', 'assembly'}:
        raise ValueError(f'Catalog line {line_number} ({label}) has unsupported kind={kind!r}')


def iter_internal_catalog(path: str) -> Iterator[Dict[str, str]]:
    with open_text(path, 'rt') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        validate_internal_catalog_header(reader.fieldnames)
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: str(value or '').strip() for key, value in row.items()}
            validate_internal_catalog_row(normalized, line_number)
            yield normalized


def resolve_parallel_backend(args: argparse.Namespace, default: str = 'process') -> str:
    backend = str(getattr(args, 'parallel_backend', default) or default).lower()
    if backend not in ('process', 'thread'):
        raise ValueError(f'Unknown parallel backend: {backend}')
    return backend


def resolve_mp_start_method(args: argparse.Namespace, default: str = 'fork') -> str:
    method = str(getattr(args, 'mp_start_method', default) or default).lower()
    if method not in ('fork', 'spawn', 'forkserver'):
        raise ValueError(f'Unknown multiprocessing start method: {method}')
    return method


def anchor_task_sort_key(fid: str) -> Tuple[str, str, str, int, int, str]:
    core = base.parse_fid(str(fid))
    return (str(core.get('sample') or ''), str(core.get('hap') or ''), str(core.get('contig') or ''), int(core.get('start0') or 0), int(core.get('end0') or 0), str(fid))


def fid_for_bed_row(chrom_field: str, start0: int, end0: int, kind: str) -> str:
    parts = str(chrom_field).split('#', 2)
    if len(parts) != 3 or any(not value for value in parts):
        raise ValueError(f'Internal CpGI BED key is invalid: {chrom_field!r}')
    return f'{chrom_field}:{start0}-{end0}'


def iter_cpgi_bed(path: str, kind: str) -> Iterator[Dict[str, object]]:
    with open_text(path, 'rt') as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith('#'):
                continue
            f = raw.rstrip('\n').split('\t')
            if len(f) != 10:
                raise ValueError(f'CpGI BED line requires exactly 10 tab-separated columns: {raw[:200]}')
            chrom = f[0]
            start0 = int(f[1])
            end0 = int(f[2])
            fid = fid_for_bed_row(chrom, start0, end0, kind)
            core = base.parse_fid(fid)
            yield {
                'fid': fid,
                'contig': str(core['contig']),
                'start0': start0,
                'end0': end0,
                'bed_name': f[3],
                'bed_size': int(float(f[4])),
                'cpg_n': int(float(f[5])),
                'gc_n': int(float(f[6])),
                'pct_cpg': float(f[7]),
                'pct_gc': float(f[8]),
                'oe': float(f[9]),
            }


def load_cpgi_fasta_dict(path: Optional[str]) -> Dict[str, str]:
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(path)
    seqs: Dict[str, List[str]] = {}
    cur = None
    with open_text(path, 'rt') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                cur = line[1:].split()[0]
                seqs[cur] = []
            else:
                if cur is None:
                    raise ValueError('FASTA sequence before header')
                seqs[cur].append(line.strip().upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def write_fasta_record(fh, name: str, seq: str, width: int = 60) -> None:
    fh.write(f'>{name}\n')
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + '\n')


def normalize_sv_contig_name(contig: str) -> str:
    return svh.normalize_sv_contig_name(contig)


def parse_sv_ins_tsv(path: Optional[str], *, contig_coordinate_base: int = 0) -> Dict[str, Dict[str, object]]:
    path = str(path or '').strip()
    if not path or not os.path.exists(path):
        raise FileNotFoundError(path)
    cache_key = f'{path}::base{int(contig_coordinate_base)}'
    cached = _SV_INS_CACHE.get(cache_key)
    if cached is not None:
        return cached
    out = svh.build_sv_index(path, contig_coordinate_base=int(contig_coordinate_base))
    _SV_INS_CACHE[cache_key] = out
    return out


def annotate_feature_sv_ins(contig: str, start0: int, end0: int, sv_index: Dict[str, Dict[str, object]]) -> Dict[str, object]:
    return svh.annotate_feature_overlaps(contig, start0, end0, sv_index)


HAL_EMPTY_FIELDS = {
    'hal_primary_chr': '',
    'hal_primary_start0': '',
    'hal_primary_end0': '',
    'hal_primary_start1': '',
    'hal_primary_end1': '',
    'hal_strand': '',
    'hal_query_cov': '',
    'hal_identity': '',
    'hal_aligned_bp': '',
    'hal_target_span_bp': '',
    'hal_block_count': '',
    'hal_multimap_n': 0,
    'hal_multimap_label': 'no_hal',
    'hal_admit_rule': '',
    'hal_primary_rank': '',
    'hal_all_intervals0': '',
    'hal_all_intervals': '',
    'hal_all_scores': '',
}


def _fmt_score(val: object, ndigits: int = 6) -> str:
    if val in ('', None, '.', 'nan', 'None'):
        return ''
    try:
        return f"{float(val):.{ndigits}f}"
    except Exception:
        return str(val)


def hal_empty_fields(label: str = 'no_hal') -> Dict[str, object]:
    out = dict(HAL_EMPTY_FIELDS)
    out['hal_multimap_label'] = label
    return out


def fid_len_from_string(fid: str) -> int:
    p = base.parse_fid(str(fid))
    return max(0, int(p['end0']) - int(p['start0']))


def parse_hal_psl_file(
    path: str,
    *,
    label: str,
    kind: str,
    min_coverage: float,
    min_identity: float,
    ambig_identity_delta: float,
    ambig_coverage_delta: float,
    ambig_aligned_bp_delta: int,
) -> Dict[str, Dict[str, object]]:


    if not path or not os.path.isfile(path):
        raise FileNotFoundError(path)
    grouped: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    with open_text(path, 'rt') as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith('#'):
                continue
            f = raw.rstrip('\n').split('\t')
            if len(f) != 22:
                raise ValueError(f'{path}: expected exactly 22 PSL-with-name fields')
            raw_fid = f[0]
            fid = str(raw_fid).strip()
            try:
                matches = int(f[1])
                mismatches = int(f[2])
                rep_matches = int(f[3])
                strand = str(f[9])
                q_start0 = int(f[12])
                q_end0 = int(f[13])
                from pancgi_contract import storage_id
                t_name = storage_id('C', str(f[14]))
                t_start0 = int(f[16])
                t_end0 = int(f[17])
                block_count = int(f[18])
                block_sizes = [int(x) for x in str(f[19]).rstrip(',').split(',') if x != '']
            except (ValueError, IndexError) as error:
                raise ValueError(f'{path}: malformed PSL record for {fid}') from error
            fid_len = fid_len_from_string(fid)
            q_span = max(1, int(q_end0) - int(q_start0))
            denom_len = max(1, int(fid_len or q_span))
            aligned_bp = sum(block_sizes) if block_sizes else max(0, matches + mismatches + rep_matches)
            identity_denom = max(1, matches + mismatches + rep_matches)
            identity = float(matches + rep_matches) / float(identity_denom)
            query_cov = float(aligned_bp) / float(denom_len)
            rec = {
                'fid': fid,
                'raw_fid': raw_fid,
                'chrom': t_name,
                'start0': int(t_start0),
                'end0': int(t_end0),
                'start1': int(t_start0) + 1,
                'end1': int(t_end0),
                'strand': strand,
                'query_cov': query_cov,
                'identity': identity,
                'aligned_bp': int(aligned_bp),
                'target_span_bp': max(1, int(t_end0) - int(t_start0)),
                'block_count': int(block_count),
            }
            grouped[fid].append(rec)

    resolved: Dict[str, Dict[str, object]] = {}
    for fid, hits in grouped.items():
        hits.sort(key=lambda r: (
            float(r['identity']),
            float(r['query_cov']),
            int(r['aligned_bp']),
            -int(r['block_count']),
            -int(r['target_span_bp']),
            str(r['chrom']),
            -int(r['start1']),
            -int(r['end1']),
        ), reverse=True)
        best = hits[0]
        label_out = 'unique' if len(hits) == 1 else 'multi_clear_primary'
        if float(best['query_cov']) < float(min_coverage) or float(best['identity']) < float(min_identity):
            label_out = 'low_quality'
        if len(hits) > 1:
            second = hits[1]
            near = (
                abs(float(best['identity']) - float(second['identity'])) <= float(ambig_identity_delta)
                and abs(float(best['query_cov']) - float(second['query_cov'])) <= float(ambig_coverage_delta)
                and abs(int(best['aligned_bp']) - int(second['aligned_bp'])) <= int(ambig_aligned_bp_delta)
            )
            if near:
                label_out = 'ambiguous_tie'
        all_intervals = []
        all_intervals0 = []
        all_scores = []
        for h in hits:
            all_intervals0.append(f"{h['chrom']}:{h['start0']}-{h['end0']}:{h['strand']}")
            all_intervals.append(f"{h['chrom']}:{h['start1']}-{h['end1']}:{h['strand']}")
            all_scores.append(f"{h['chrom']}:{h['start1']}-{h['end1']}|cov={_fmt_score(h['query_cov'])}|id={_fmt_score(h['identity'])}|bp={h['aligned_bp']}")
        resolved[fid] = {
            'hal_primary_chr': best['chrom'] if label_out not in ('low_quality', 'ambiguous_tie') else '',
            'hal_primary_start0': best['start0'] if label_out not in ('low_quality', 'ambiguous_tie') else '',
            'hal_primary_end0': best['end0'] if label_out not in ('low_quality', 'ambiguous_tie') else '',
            'hal_primary_start1': best['start1'] if label_out not in ('low_quality', 'ambiguous_tie') else '',
            'hal_primary_end1': best['end1'] if label_out not in ('low_quality', 'ambiguous_tie') else '',
            'hal_strand': best['strand'] if label_out not in ('low_quality', 'ambiguous_tie') else '',
            'hal_query_cov': _fmt_score(best['query_cov']),
            'hal_identity': _fmt_score(best['identity']),
            'hal_aligned_bp': int(best['aligned_bp']),
            'hal_target_span_bp': int(best['target_span_bp']),
            'hal_block_count': int(best['block_count']),
            'hal_multimap_n': len(hits),
            'hal_multimap_label': label_out,
            'hal_primary_rank': 1,
            'hal_all_intervals0': ';'.join(all_intervals0),
            'hal_all_intervals': ';'.join(all_intervals),
            'hal_all_scores': ';'.join(all_scores),
        }
    return resolved


def hal_fields_for_feature(fid: str, hal_idx: Dict[str, Dict[str, object]], *, kind: str) -> Dict[str, object]:
    if str(kind) == 'reference_primary':
        try:
            p = base.parse_fid(fid)
            out = dict(HAL_EMPTY_FIELDS)
            out.update({
                'hal_primary_chr': str(p['contig']),
                'hal_primary_start0': int(p['start0']),
                'hal_primary_end0': int(p['end0']),
                'hal_primary_start1': int(p['start0']) + 1,
                'hal_primary_end1': int(p['end0']),
                'hal_strand': '+',
                'hal_query_cov': '1.000000',
                'hal_identity': '1.000000',
                'hal_aligned_bp': int(p['end0']) - int(p['start0']),
                'hal_target_span_bp': int(p['end0']) - int(p['start0']),
                'hal_block_count': 1,
                'hal_multimap_n': 1,
                'hal_multimap_label': 'self_primary',
                'hal_admit_rule': 'self_primary',
                'hal_primary_rank': 1,
                'hal_all_intervals0': f"{p['contig']}:{int(p['start0'])}-{int(p['end0'])}:+",
                'hal_all_intervals': f"{p['contig']}:{int(p['start0']) + 1}-{int(p['end0'])}:+",
                'hal_all_scores': f"{p['contig']}:{int(p['start0']) + 1}-{int(p['end0'])}|cov=1.000000|id=1.000000|bp={int(p['end0']) - int(p['start0'])}",
            })
            return out
        except Exception as exc:
            raise ValueError('Invalid primary-reference feature identity') from exc
    return dict(hal_idx.get(str(fid)) or hal_empty_fields('no_hal'))


def is_hal_usable_for_ref(feat: Dict[str, object]) -> bool:
    label = str(feat.get('hal_multimap_label') or '')
    if label not in ('unique', 'self_primary'):
        return False
    return bool(str(feat.get('hal_primary_chr') or '') and feat.get('hal_primary_start1') not in ('', None) and feat.get('hal_primary_end1') not in ('', None))


def is_hal_usable_for_nonref(feat: Dict[str, object], args: argparse.Namespace) -> bool:
    if str(feat.get('hal_multimap_label') or '') != 'unique':
        return False
    if not is_hal_usable_for_ref(feat):
        return False
    try:
        return float(feat.get('hal_query_cov') or 0.0) >= float(getattr(args, 'hal_nonref_min_coverage', getattr(args, 'hal_min_coverage', 0.5)))
    except Exception:
        return False


def feature_hal_member_fields(feat: Dict[str, object]) -> Dict[str, object]:
    return {k: feat.get(k, v) for k, v in HAL_EMPTY_FIELDS.items()}


def hal_sv_prefilter_decision(h: Dict[str, object], sv_ann: Dict[str, object]) -> Tuple[bool, str, str]:

    h_label = str(h.get('hal_multimap_label') or '')
    sv_types = {x.strip().upper() for x in str(sv_ann.get('sv_overlap_types') or '').split(';') if x.strip()}
    has_ins = 'INS' in sv_types

    if has_ins:
        if h_label == 'unique':
            return True, 'hal_unique_cov_ge_0.5_with_ins', ''
        if h_label == 'multi_clear_primary':
            return True, 'sv_ins_primary_hal_multimap', ''
        if h_label == 'ambiguous_tie':
            return True, 'sv_ins_primary_hal_ambiguous', ''
        if h_label == 'low_quality':
            return True, 'sv_ins_primary_hal_cov_lt_0.5', ''
        if h_label == 'no_hal':
            return True, 'sv_ins_primary_no_hal', ''
        return True, f'sv_ins_primary_hal_{h_label or "unknown"}', ''

    if h_label == 'unique':
        return True, 'hal_unique_cov_ge_0.5', ''
    if h_label == 'multi_clear_primary':
        return False, '', 'hal_multimap_no_ins'
    if h_label == 'ambiguous_tie':
        return False, '', 'hal_ambiguous_no_ins'
    if h_label == 'low_quality':
        return False, '', 'hal_cov_lt_0.5_no_ins'
    if h_label == 'no_hal':
        return False, '', 'no_hal_no_ins'
    return False, '', 'no_usable_hal_or_ins_evidence'


def sv_best_interval_for_type(feat: Dict[str, object], svtype: str) -> Optional[Dict[str, object]]:
    details = parse_sv_overlap_detail_json(feat)
    wanted = str(svtype).upper()
    rows = [d for d in details if str(d.get('svtype') or '').upper() == wanted and str(d.get('chrom') or '')]
    if not rows:
        return None
    rows.sort(key=lambda d: (
        int(d.get('overlap_bp') or 0),
        int(d.get('svlen_abs') or 0),
        str(d.get('id') or ''),
    ), reverse=True)
    d = rows[0]
    chrom = str(d.get('chrom') or '')
    pos1 = _safe_int(d.get('pos1'))
    end1 = _safe_int(d.get('vcf_end1'), pos1)
    if not chrom or pos1 is None:
        return None
    return {'chrom': chrom, 'start1': pos1, 'end1': end1 if end1 is not None else pos1, 'id': str(d.get('id') or ''), 'svtype': wanted}


def nonref_site_seed_from_row(row: Dict[str, object], args: argparse.Namespace) -> Optional[Dict[str, object]]:
    hal_ok = is_hal_usable_for_nonref(row, args)
    ins_seed = sv_best_interval_for_type(row, 'INS')
    del_seed = sv_best_interval_for_type(row, 'DEL')
    if ins_seed is not None:
        return {
            'source': 'sv+hal' if hal_ok else 'sv_only',
            'chrom': str(ins_seed['chrom']),
            'start1': int(ins_seed['start1']),
            'end1': int(ins_seed['end1']),
        }
    if del_seed is not None:


        if not hal_ok:
            return None
        try:
            return {'source': 'sv+hal', 'chrom': str(row.get('hal_primary_chr')), 'start1': int(row.get('hal_primary_start1')), 'end1': int(row.get('hal_primary_end1'))}
        except Exception:
            return None
    if hal_ok:
        try:
            return {'source': 'hal_only', 'chrom': str(row.get('hal_primary_chr')), 'start1': int(row.get('hal_primary_start1')), 'end1': int(row.get('hal_primary_end1'))}
        except Exception:
            return None
    return None


def feature_sv_member_fields(feat: Dict[str, object]) -> Dict[str, object]:
    return {
        'sv_overlap_n': feat.get('sv_overlap_n', 0),
        'sv_overlap_ids': feat.get('sv_overlap_ids', ''),
        'sv_overlap_types': feat.get('sv_overlap_types', ''),
        'sv_overlap_classes': feat.get('sv_overlap_classes', ''),
        'sv_overlap_primary_sites': feat.get('sv_overlap_primary_sites', ''),
        'sv_overlap_primary_intervals': feat.get('sv_overlap_primary_intervals', ''),
        'sv_overlap_asm_intervals': feat.get('sv_overlap_asm_intervals', ''),
        'sv_overlap_bp': feat.get('sv_overlap_bp', ''),
        'sv_overlap_pct_cpgi': feat.get('sv_overlap_pct_cpgi', ''),
        'sv_overlap_pct_sv': feat.get('sv_overlap_pct_sv', ''),
        'sv_overlap_TR': feat.get('sv_overlap_TR', ''),
        'sv_overlap_CONFORMATION': feat.get('sv_overlap_CONFORMATION', ''),
        'sv_overlap_SD': feat.get('sv_overlap_SD', ''),
        'sv_overlap_ITYPE_N': feat.get('sv_overlap_ITYPE_N', ''),
        'sv_overlap_DTYPE_N': feat.get('sv_overlap_DTYPE_N', ''),
        'sv_overlap_FAM_N': feat.get('sv_overlap_FAM_N', ''),
        'sv_overlap_source': 'provided_sv_callset' if int(feat.get('sv_overlap_n', 0) or 0) > 0 else '',
        'sv_overlap_detail_json': feat.get('sv_overlap_detail_json', '[]'),
        'sv_ins_overlap_class': feat.get('sv_ins_overlap_class', ''),
        'sv_ins_n': feat.get('sv_ins_n', 0),
        'sv_ins_ids': feat.get('sv_ins_ids', ''),
        'sv_ins_primary_sites': feat.get('sv_ins_primary_sites', ''),
        'sv_ins_primary_intervals': feat.get('sv_ins_primary_intervals', ''),
        'sv_ins_asm_intervals': feat.get('sv_ins_asm_intervals', ''),
        'sv_ins_source': feat.get('sv_ins_source', ''),
        'sv_ins_longest_id': feat.get('sv_ins_longest_id', ''),
        'sv_ins_longest_chrom': feat.get('sv_ins_longest_chrom', ''),
        'sv_ins_longest_site1': feat.get('sv_ins_longest_site1', ''),
        'sv_ins_longest_start1': feat.get('sv_ins_longest_start1', ''),
        'sv_ins_longest_end1': feat.get('sv_ins_longest_end1', ''),
        'sv_ins_longest_len': feat.get('sv_ins_longest_len', ''),
        'sv_ins_overlap_bp': feat.get('sv_ins_overlap_bp', ''),
        'sv_ins_overlap_pct_cpgi': feat.get('sv_ins_overlap_pct_cpgi', ''),
        'sv_ins_overlap_pct_sv': feat.get('sv_ins_overlap_pct_sv', ''),
        'sv_ins_detail_json': feat.get('sv_ins_detail_json', '[]'),
    }


def parse_sv_overlap_detail_json(rec: Dict[str, object]) -> List[Dict[str, object]]:
    obj = json.loads(rec['sv_overlap_detail_json'])
    if not isinstance(obj, list) or not all(isinstance(x, dict) for x in obj):
        raise ValueError('Invalid serialized SV overlap evidence')
    return obj


def _safe_int(val: object, default: Optional[int] = None) -> Optional[int]:
    try:
        if val in ('', None, '.', 'nan', 'None'):
            return default
        return int(val)
    except Exception:
        try:
            return int(float(str(val)))
        except Exception:
            return default


def choose_coordinate_ref_candidate(feat: Dict[str, object], ref_interval_index: Dict[str, List[Tuple[int, int, str]]]) -> Optional[Dict[str, object]]:
    chrom = str(feat.get('primary_chr') or '').strip()
    start1 = _safe_int(feat.get('primary_start1'))
    end1 = _safe_int(feat.get('primary_end1'))
    if not chrom or start1 is None or end1 is None:
        return None
    candidates = query_ref_intervals_overlap(ref_interval_index, chrom, start1, end1)
    if not candidates:
        return None
    feat_len = max(1, int(end1) - int(start1) + 1)
    scored = []
    for s, e, fid in candidates:
        ref_len = max(1, int(e) - int(s) + 1)
        ov = interval_overlap_bp1(start1, end1, int(s), int(e))
        ro = float(ov) / float(min(feat_len, ref_len)) if min(feat_len, ref_len) > 0 else 0.0
        size_ratio = float(min(feat_len, ref_len)) / float(max(feat_len, ref_len)) if max(feat_len, ref_len) > 0 else 0.0
        gap = interval_gap_bp1(start1, end1, int(s), int(e))
        scored.append((ov, ro, size_ratio, -gap, -ref_len, str(fid), (s, e, fid)))
    scored.sort(reverse=True)
    best = scored[0]
    second = scored[1] if len(scored) > 1 else None
    return {
        'best_ref_fid': str(best[-1][2]),
        'best_ref_iv': best[-1],
        'best_overlap_bp': int(best[0]),
        'best_ro': float(best[1]),
        'best_size_ratio': float(best[2]),
        'candidate_n': len(candidates),
        'second_ref_fid': '' if second is None else str(second[-1][2]),
        'second_overlap_bp': '' if second is None else int(second[0]),
        'second_ro': '' if second is None else float(second[1]),
        'second_size_ratio': '' if second is None else float(second[2]),
        'source': 'coordinate_overlap',
    }


def choose_hal_ref_candidate(feat: Dict[str, object], ref_interval_index: Dict[str, List[Tuple[int, int, str]]], primary_refs: Dict[str, Dict[str, object]]) -> Optional[Dict[str, object]]:
    if not is_hal_usable_for_ref(feat):
        return None
    try:
        q = {
            'primary_chr': feat.get('hal_primary_chr'),
            'primary_start1': feat.get('hal_primary_start1'),
            'primary_end1': feat.get('hal_primary_end1'),
        }
        cand = choose_coordinate_ref_candidate(q, ref_interval_index)
    except Exception as exc:
        raise ValueError('Invalid accepted HAL projection') from exc
    if cand is None:
        return None
    cand['source'] = 'hal_liftover'
    return cand


def choose_sv_site_ref_candidate(feat: Dict[str, object], ref_interval_index: Dict[str, List[Tuple[int, int, str]]], primary_refs: Dict[str, Dict[str, object]]) -> Optional[Dict[str, object]]:
    details = json.loads(feat['sv_ins_detail_json'])
    if not isinstance(details, list) or len(details) != int(feat['sv_ins_n']):
        raise ValueError('Insertion details do not match the feature inventory')
    matched = []
    all_ref_ids = set()
    for ins in details:
        chrom, site = str(ins['chrom']), int(ins['pos1'])
        candidates = query_ref_intervals_overlap(ref_interval_index, chrom, site, site)
        if candidates:
            matched.append((ins, candidates))
            all_ref_ids.update(fid for _, _, fid in candidates)
    if not matched:
        return None
    ins, candidates = max(matched, key=lambda item: (
        int(item[0]['svlen_abs']), int(item[0]['overlap_bp']),
        int(item[0]['span_bp']), str(item[0]['id'])))
    site1 = int(ins['pos1'])
    start1 = end1 = start_q = end_q = site1
    scored = []
    for s, e, fid in candidates:
        ref_len = max(1, int(e) - int(s) + 1)
        contains_site = 1 if int(s) <= int(site1) <= int(e) else 0
        ov = interval_overlap_bp1(start_q, end_q, int(s), int(e))
        margin = min(abs(int(site1) - int(s)), abs(int(e) - int(site1))) if contains_site else -1
        scored.append((contains_site, ov, margin, -ref_len, str(fid), (s, e, fid)))
    scored.sort(reverse=True)
    best = scored[0]
    second = scored[1] if len(scored) > 1 else None
    return {
        'best_ref_fid': str(best[-1][2]),
        'best_ref_iv': best[-1],
        'candidate_n': len(all_ref_ids),
        'second_ref_fid': '' if second is None else str(second[-1][2]),
        'insert_site1': site1,
        'insert_interval_start1': start1 if start1 is not None else site1,
        'insert_interval_end1': end1 if end1 is not None else site1,
        'source': 'provided_sv_callset',
        'reference_ins_id': str(ins['id']),
        'all_matched_ref_ids': sorted(all_ref_ids),
    }


def cmd_make_cpgi_fasta(args: argparse.Namespace) -> None:
    rows = list(iter_internal_catalog(args.catalog))
    if not rows:
        raise ValueError('Internal sample catalog contains no rows')
    halrt.inspect_image(args.docker_bin, args.docker_image)
    log_dir = os.path.abspath(args.log_dir)
    os.makedirs(log_dir, exist_ok=True)
    contig_rows = mapping.read_tsv(args.contigs, mapping.CONTIG_COLUMNS)
    contig_index = {}
    for contig_row in contig_rows:
        if contig_row['mapping_status'] != 'confirmed':
            raise ValueError(f"CpGI FASTA extraction requires mapping_status=confirmed: line {contig_row['_line_number']}")
        key = (contig_row['genome_id'], contig_row['contig_id'])
        if key in contig_index:
            raise ValueError(f'Duplicate HAL contig mapping: {key}')
        contig_index[key] = contig_row['hal_sequence']

    def worker(row: Dict[str, str]) -> Dict[str, object]:
        return halrt.extract_cpgi_fasta_from_hal(
            row,
            args.hal,
            args.docker_bin,
            args.docker_image,
            args.hal2fasta,
            os.path.join(log_dir, f"{row['label']}.hal2fasta.log"),
            contig_index,
        )

    results: List[Dict[str, object]] = []
    threads = max(1, int(args.threads))
    if threads == 1:
        for row in rows:
            results.append(worker(row))
    else:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            results.extend(executor.map(worker, rows))
    order = {row['label']: index for index, row in enumerate(rows)}
    results.sort(key=lambda row: order[str(row['label'])])
    missing_rows: List[Tuple[str, str, str, int, int]] = []
    for result in results:
        for fid, expected_bp, observed_bp in result.get('missing', []):
            reason = 'missing_contig_or_interval' if int(observed_bp) == 0 else 'sequence_length_mismatch'
            missing_rows.append((str(result['label']), str(fid), reason, int(expected_bp), int(observed_bp)))
    if args.missing_report and missing_rows:
        ensure_parent(args.missing_report)
        with open_text(args.missing_report, 'wt') as out:
            w = csv.writer(out, delimiter='\t', lineterminator='\n')
            w.writerow(['label', 'fid', 'reason', 'expected_bp', 'observed_bp'])
            for row in missing_rows:
                w.writerow(row)

    if missing_rows:
        raise RuntimeError(
            f'CpGI FASTA construction failed strict completeness: {len(missing_rows)} BED intervals '
            'could not be extracted from HAL'
        )
    print(json_pretty({
        'sample_fastas_written': len(results),
        'records_written': sum(int(result['records']) for result in results),
        'missing_regions': len(missing_rows),
        'sequence_source': 'hal2fasta_upper',
    }), file=sys.stderr)


def build_feature_genome(row, args, out_fh, exc_w):
    n_written = 0
    n_excluded = 0
    n_partial = 0
    label = row['label']
    kind = row['kind']
    bed_path = row['bed']
    path_bed_dir = row['path_bed_dir']
    cpgi_fa = row['cpgi_fa']
    sv_tsv = row['sv_tsv']
    hal_psl = os.path.join(args.hal_psl_dir, f'{label}.to_primary.psl')

    for input_path, input_name in [(bed_path, 'bed'), (cpgi_fa, 'cpgi_fa')]:
        if not os.path.isfile(input_path):
            raise FileNotFoundError(f'{label}: required {input_name} not found: {input_path}')
    if not os.path.isdir(path_bed_dir):
        raise NotADirectoryError(f'{label}: required path_bed_dir not found: {path_bed_dir}')
    if kind != 'reference_primary' and not os.path.isfile(hal_psl):
        raise FileNotFoundError(f'{label}: required HAL PSL not found: {hal_psl}')
    if kind == 'assembly' and not os.path.isfile(sv_tsv):
        raise FileNotFoundError(f'{label}: required SV table not found: {sv_tsv}')

    bed_rows = list(iter_cpgi_bed(bed_path, kind))
    bed_idx = {str(record['fid']): record for record in bed_rows}
    by_contig: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for record in bed_rows:
        by_contig[str(record['contig'])].append(
            (int(record['start0']), int(record['end0']), str(record['fid']))
        )

    seq_from_cpgi = load_cpgi_fasta_dict(cpgi_fa)
    missing_sequences = [fid for fid in bed_idx if fid not in seq_from_cpgi]
    if missing_sequences:
        raise RuntimeError(
            f'{label}: CpGI FASTA lacks {len(missing_sequences)} BED features; '
            f'first missing IDs: {missing_sequences[:5]}'
        )
    sv_idx = (
        parse_sv_ins_tsv(
            sv_tsv,
            contig_coordinate_base=int(getattr(args, 'sv_contig_coordinate_base', 0)),
        )
        if kind == 'assembly'
        else {}
    )
    hal_idx = (
        parse_hal_psl_file(
            hal_psl,
            label=label,
            kind=kind,
            min_coverage=float(getattr(args, 'hal_min_coverage', 0.5)),
            min_identity=float(getattr(args, 'hal_min_identity', 0.0)),
            ambig_identity_delta=float(getattr(args, 'hal_ambig_identity_delta', 0.001)),
            ambig_coverage_delta=float(getattr(args, 'hal_ambig_coverage_delta', 0.01)),
            ambig_aligned_bp_delta=int(getattr(args, 'hal_ambig_aligned_bp_delta', 10)),
        )
        if kind != 'reference_primary'
        else {}
    )

    mapped_all: Dict[str, Dict[str, object]] = {}
    for contig, features in by_contig.items():
        graph_sample = row['graph_sample']
        graph_hap = row['graph_hap']
        pathbed = pbh.resolve_pathbed_file(
            path_bed_dir,
            sample=graph_sample,
            hap=graph_hap,
            contig=contig,
        )
        if pathbed is None:
            raise FileNotFoundError(
                f'{label}: no pathBED for explicit graph key '
                f'{graph_sample}#{graph_hap} contig {contig}'
            )
        mapped_all.update(
            pbh.map_features_on_contig(
                contig,
                features,
                pathbed,
                flank_bp=int(getattr(args, 'flank_bp', 1000)),
                flank_max_steps=int(getattr(args, 'flank_max_steps', 32)),
            )
        )

    for fid in bed_idx:
        if fid not in mapped_all:
            continue
        mapping = mapped_all[fid]
        if not mapping['steps']:
            exc_w.writerow([label, kind, fid, 'no_path_overlap'])
            n_excluded += 1
            continue
        if float(mapping['graph_cov']) < args.min_graph_cov:
            exc_w.writerow([label, kind, fid, f'graph_cov_lt_{args.min_graph_cov}'])
            n_excluded += 1
            n_partial += 1
            continue

        core = base.parse_fid(fid)
        hal = hal_fields_for_feature(fid, hal_idx, kind=kind)
        sv_annotation = annotate_feature_sv_ins(
            str(core['contig']),
            int(core['start0']),
            int(core['end0']),
            sv_idx,
        )
        if kind != 'reference_primary':
            admit, admit_rule, reject_reason = hal_sv_prefilter_decision(hal, sv_annotation)
            if not admit:
                exc_w.writerow([label, kind, fid, reject_reason])
                n_excluded += 1
                continue
            hal['hal_admit_rule'] = admit_rule

        bed_record = bed_idx[fid]
        sequence = seq_from_cpgi[fid]
        rec = {
            'fid': fid,
            'label': label,
            'kind': kind,
            'sample': core['sample'],
            'hap': core['hap'],
            'contig': core['contig'],
            'asm_start0': core['start0'],
            'asm_end0': core['end0'],
            'asm_len': int(core['end0']) - int(core['start0']),
            'bed_name': bed_record['bed_name'],
            'bed_size': bed_record['bed_size'],
            'cpg_n': bed_record['cpg_n'],
            'gc_n': bed_record['gc_n'],
            'pct_cpg': bed_record['pct_cpg'],
            'pct_gc': bed_record['pct_gc'],
            'oe': bed_record['oe'],
            'graph_path': mapping['graph_path'],
            'graph_bp': mapping['graph_bp'],
            'graph_cov': mapping['graph_cov'],
            'graph_nsteps': mapping['graph_nsteps'],
            'steps': mapping['steps'],
            'nodeints': mapping['nodeints'],
            'left_flank_steps': mapping.get('left_flank_steps', []),
            'right_flank_steps': mapping.get('right_flank_steps', []),
            'left_flank_bp': mapping.get('left_flank_bp', 0),
            'right_flank_bp': mapping.get('right_flank_bp', 0),
            'left_flank_nsteps': mapping.get('left_flank_nsteps', 0),
            'right_flank_nsteps': mapping.get('right_flank_nsteps', 0),
            'source_seq': sequence,
            'source_seq_len': len(sequence),
            'source_seq_source': 'cpgi_fa',
            'primary_chr': hal.get('hal_primary_chr'),
            'primary_start1': hal.get('hal_primary_start1'),
            'primary_end1': hal.get('hal_primary_end1'),
            'nonref_only': 1 if (kind == 'assembly' and not hal.get('hal_primary_chr')) else 0,
            'has_primary_coordinate': 1 if hal.get('hal_primary_chr') else 0,
            'position_coordinate_backend': 'hal_sv',
        }
        rec.update(sv_annotation)
        rec.update(hal)
        out_fh.write(json_compact(rec) + '\n')
        n_written += 1

    return {'features_written': n_written, 'excluded': n_excluded,
            'partial_graph_coverage_excluded': n_partial, 'input_records': len(bed_rows)}


def cmd_build_features_prod(args: argparse.Namespace) -> None:
    from pancgi_feature_build import run
    run(args)


def iter_feature_jsonl(path: str) -> Iterator[Dict[str, object]]:
    if str(path).endswith('.sqlite'):
        from pancgi_features import iter_raw
        yield from iter_raw(path)
        return
    with open_text(path, 'rt') as fh:
        for raw in fh:
            if raw.strip():
                yield json.loads(raw)


def feature_bp(rec: Dict[str, object]) -> int:
    if rec.get('source_seq_len') is not None:
        return int(rec['source_seq_len'])
    total = 0
    for node_id, s, e in rec['nodeints']:
        total += int(e) - int(s)
    return total


def shingle_strings(tokens: Sequence[str], k: int, max_mid_anchors: int = 8) -> List[str]:
    if not tokens:
        return []
    kk = min(k, len(tokens))
    if len(tokens) <= kk:
        return ['|'.join(tokens)]
    arr = ['|'.join(tokens[i:i + kk]) for i in range(len(tokens) - kk + 1)]
    if len(arr) > max_mid_anchors:
        idxs = sorted(set(int(round(i * (len(arr) - 1) / (max_mid_anchors - 1))) for i in range(max_mid_anchors)))
        arr = [arr[i] for i in idxs]
    return arr


def prep_feature(rec: Dict[str, object], anchor_k: int, max_mid_anchors: int = 8) -> Dict[str, object]:
    steps = rec['steps']
    tokens = [f"{'<' if int(step[1]) else '>'}{int(step[0])}" for step in steps]
    kk = min(anchor_k, len(tokens)) if tokens else 0
    left = tuple(tokens[:kk]) if kk else tuple()
    right = tuple(tokens[-kk:]) if kk else tuple()
    mids = shingle_strings(tokens, anchor_k, max_mid_anchors=max_mid_anchors)
    anchors = []
    if left:
        anchors.append('L:' + '|'.join(left))
    if right:
        anchors.append('R:' + '|'.join(right))
    anchors.extend('M:' + m for m in mids)
    shingles = set('M:' + m for m in mids)
    node_map: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for node_id, s, e in rec['nodeints']:
        node_map[int(node_id)].append((int(s), int(e)))
    p = dict(rec)
    p['_tokens'] = tokens
    p['_left'] = left
    p['_right'] = right
    p['_anchors'] = anchors
    p['_shingles'] = shingles
    p['_node_map'] = node_map
    p['_bp'] = feature_bp(rec)
    return p


def collect_shingle_df(features_path: str, anchor_k: int, max_mid_anchors: int) -> Tuple[int, Dict[str, int]]:
    if str(features_path).endswith('.sqlite'):
        from pancgi_features import frequencies
        return frequencies(features_path, anchor_k, max_mid_anchors)
    n = 0
    df: Counter = Counter()
    for rec in iter_feature_jsonl(features_path):
        p = prep_feature(rec, anchor_k=anchor_k, max_mid_anchors=max_mid_anchors)
        df.update(set(p['_shingles']))
        n += 1
    return n, dict(df)


def shingle_weights(n_features: int, df: Dict[str, int]) -> Dict[str, float]:
    return {s: math.log((n_features + 1.0) / (d + 1.0)) + 1.0 for s, d in df.items()}


def weighted_jaccard(a: set, b: set, weights: Optional[Dict[str, float]]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    if not weights:
        inter = len(a & b)
        union = len(a | b)
        return inter / union if union else 0.0

    inter_terms: List[float] = []
    union_terms: List[float] = []
    for k in sorted(a | b):
        w = float(weights.get(k, 1.0))
        if k in a and k in b:
            inter_terms.append(w)
        union_terms.append(w)
    inter = math.fsum(inter_terms)
    union = math.fsum(union_terms)
    return inter / union if union else 0.0


def interval_intersection_len(a: List[Tuple[int, int]], b: List[Tuple[int, int]]) -> int:
    i = j = 0
    total = 0
    a2 = sorted(a)
    b2 = sorted(b)
    while i < len(a2) and j < len(b2):
        s1, e1 = a2[i]
        s2, e2 = b2[j]
        ov = min(e1, e2) - max(s1, s2)
        if ov > 0:
            total += ov
        if e1 <= e2:
            i += 1
        else:
            j += 1
    return total


def shared_bp(a: Dict[str, object], b: Dict[str, object]) -> int:
    total = 0
    common_nodes = set(a['_node_map'].keys()) & set(b['_node_map'].keys())
    for node_id in common_nodes:
        total += interval_intersection_len(a['_node_map'][node_id], b['_node_map'][node_id])
    return total


def kmer_set(seq: str, k: int) -> set:
    if seq is None:
        return set()
    if len(seq) < k:
        return {seq} if seq else set()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def kmer_jaccard(seq1: Optional[str], seq2: Optional[str], k: int = 9) -> Optional[float]:
    if not seq1 or not seq2:
        return None
    a = kmer_set(seq1, k)
    b = kmer_set(seq2, k)
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def external_long_identity(seq1: Optional[str], seq2: Optional[str], args: argparse.Namespace) -> Optional[float]:
    if not seq1 or not seq2:
        raise ValueError('Missing sequence for WFA alignment')
    cmd_tpl = str(getattr(args, 'very_long_external_template', '') or '').strip()
    if not cmd_tpl:
        raise ValueError('The declared WFA command is required')
    with tempfile.TemporaryDirectory(prefix='cpgi_longalign_') as td:
        fa1 = os.path.join(td, 'seq1.fa')
        fa2 = os.path.join(td, 'seq2.fa')
        with open(fa1, 'wt') as fh:
            write_fasta_record(fh, 'seq1', seq1)
        with open(fa2, 'wt') as fh:
            write_fasta_record(fh, 'seq2', seq2)
        cmd = cmd_tpl.format(seq1=fa1, seq2=fa2, outdir=td)
        try:
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        except Exception as e:
            raise RuntimeError(f'very-long external command failed: {e}')
        stdout = (proc.stdout or '').strip()
        if not stdout:
            raise RuntimeError('WFA produced no result')
        try:
            result = json.loads(stdout)
        except (ValueError, TypeError) as error:
            raise RuntimeError('WFA output must be one JSON object') from error
        if not isinstance(result, dict):
            raise RuntimeError('WFA output must be one JSON object')
        identity = result.get('identity')
        if isinstance(identity, bool) or not isinstance(identity, (int, float)):
            raise RuntimeError('WFA output requires a numeric identity')
        if not math.isfinite(identity) or not 0 <= identity <= 1:
            raise RuntimeError('WFA identity must be finite and between zero and one')
        return float(identity)


def long_seq_similarity(seq1: Optional[str], seq2: Optional[str], args: argparse.Namespace) -> Optional[float]:
    if not seq1 or not seq2:
        return None
    thr = int(getattr(args, 'very_long_threshold', 0) or 0)
    if thr <= 0 or max(len(seq1), len(seq2)) < thr:
        return None
    backend = str(args.very_long_backend)
    if backend != 'external':
        raise ValueError('Publication allele clustering requires the declared WFA backend for long sequences')
    sim = external_long_identity(seq1, seq2, args)
    if sim is None or not math.isfinite(sim) or not 0 <= sim <= 1:
        raise RuntimeError('Invalid WFA identity; refusing alternate alignment')
    return sim


def parasail_available() -> bool:
    return parasail is not None


def parasail_matrix(match: int, mismatch: int):
    if parasail is None:
        return None
    return parasail.matrix_create("ACGTN", int(match), int(mismatch))


def parasail_semiglobal_identity(
    seq1: Optional[str],
    seq2: Optional[str],
    *,
    mode: str = 'sg',
    match: int = 2,
    mismatch: int = -3,
    gap_open: int = 5,
    gap_extend: int = 2,
) -> Optional[float]:
    if parasail is None or not seq1 or not seq2:
        raise RuntimeError('Parasail and both sequences are required')
    matrix = parasail_matrix(match, mismatch)
    fn = getattr(parasail, f'{mode}_stats_scan_32', None)
    if fn is None:
        raise RuntimeError('The declared 32-bit Parasail implementation is unavailable')
    result = fn(seq1, seq2, int(gap_open), int(gap_extend), matrix)
    if result.saturated:
        raise ArithmeticError('Parasail alignment saturated; no alternate backend was used')
    matches, aln_len = result.matches, result.length
    if aln_len <= 0 or not 0 <= matches <= aln_len:
        raise RuntimeError('Invalid Parasail alignment statistics')
    return float(matches) / float(aln_len)


def resolve_seq_backend(args: argparse.Namespace) -> str:
    backend = str(getattr(args, 'seq_backend', 'parasail')).lower()
    if backend != 'parasail':
        raise ValueError('Publication allele clustering requires Parasail')
    return backend


def seq_similarity(
    seq1: Optional[str],
    seq2: Optional[str],
    *,
    args: argparse.Namespace,
) -> Optional[float]:
    if not seq1 or not seq2:
        raise ValueError('Missing allele sequence')
    long_sim = long_seq_similarity(seq1, seq2, args)
    if long_sim is not None:
        return long_sim
    backend = resolve_seq_backend(args)
    if backend == 'parasail':
        if not parasail_available():
            raise RuntimeError('seq-backend=parasail but parasail is not installed. Please `pip install parasail`.')
        return parasail_semiglobal_identity(
            seq1,
            seq2,
            mode=str(getattr(args, 'parasail_mode', 'sg')),
            match=int(getattr(args, 'parasail_match', 2)),
            mismatch=int(getattr(args, 'parasail_mismatch', -3)),
            gap_open=int(getattr(args, 'parasail_gap_open', 5)),
            gap_extend=int(getattr(args, 'parasail_gap_extend', 2)),
        )
    raise ValueError(f'Unknown seq backend: {backend}')


def allele_freq_class(asm_hap_n: int, total_assemblies: int, major_flag: int, has_ref: int = 0) -> str:
    if asm_hap_n <= 0:
        return 'ref_only' if int(has_ref) == 1 else 'absent'
    if asm_hap_n == 1:
        return 'singleton'
    if total_assemblies > 0 and asm_hap_n >= total_assemblies:
        return 'shared'
    if int(major_flag) == 1:
        return 'major'
    return 'polymorphic'


def compare_feature_to_rep(
    feat: Dict[str, object],
    rep: Dict[str, object],
    args: argparse.Namespace,
    sh_weights: Optional[Dict[str, float]],
) -> Dict[str, object]:
    sh_bp = shared_bp(feat, rep)
    len_a = int(feat['_bp'])
    len_b = int(rep['_bp'])
    ro = min(sh_bp / max(1, len_a), sh_bp / max(1, len_b))
    size_ratio = min(len_a, len_b) / max(1, max(len_a, len_b))
    ctx = weighted_jaccard(set(feat['_shingles']), set(rep['_shingles']), sh_weights)
    boundary_match = int(feat['_left'] == rep['_left'] or feat['_right'] == rep['_right'])
    weak_seq = None
    rule = None
    rank = 0

    if ro >= args.exact_ro and size_ratio >= args.exact_size and ctx >= args.exact_ctx:
        rule = 'exact_graph'
        rank = 3
    elif ro >= args.ro_min and ctx >= args.ro_ctx:
        rule = 'ro_graph'
        rank = 2
    elif ro >= args.szro_ro and size_ratio >= args.szro_size and (ctx >= args.szro_ctx or boundary_match):
        if args.weak_seq_gate > 0:
            weak_seq = kmer_jaccard(feat.get('source_seq'), rep.get('source_seq'), k=args.kmer)
            if weak_seq is not None and weak_seq < args.weak_seq_gate:
                rule = None
                rank = 0
            else:
                rule = 'szro_graph'
                rank = 1
        else:
            rule = 'szro_graph'
            rank = 1

    composite = (rank * 10.0) + (ro * 3.0) + (ctx * 2.0) + size_ratio + (0.5 * (weak_seq if weak_seq is not None else 0.0))
    return {
        'rule': rule,
        'rule_rank': rank,
        'shared_bp': sh_bp,
        'ro_graph': ro,
        'size_ratio': size_ratio,
        'ctx_graph': ctx,
        'boundary_match': boundary_match,
        'weak_seq': weak_seq,
        'composite': composite,
    }


def pairwise_graph_similarity(
    a: Dict[str, object],
    b: Dict[str, object],
    sh_weights: Optional[Dict[str, float]],
) -> Dict[str, float]:
    sh_bp = shared_bp(a, b)
    len_a = int(a['_bp'])
    len_b = int(b['_bp'])
    ro = min(sh_bp / max(1, len_a), sh_bp / max(1, len_b))
    size_ratio = min(len_a, len_b) / max(1, max(len_a, len_b))
    ctx = weighted_jaccard(set(a['_shingles']), set(b['_shingles']), sh_weights)
    boundary = 1.0 if (a['_left'] == b['_left'] or a['_right'] == b['_right']) else 0.0
    score = (ro * 4.0) + (ctx * 2.0) + size_ratio + (0.25 * boundary)
    return {
        'shared_bp': sh_bp,
        'ro_graph': ro,
        'size_ratio': size_ratio,
        'ctx_graph': ctx,
        'boundary_match': boundary,
        'score': score,
    }

def pick_display_candidate(fids: Sequence[str], feat_idx: Dict[str, Dict[str, object]], default_fid: str) -> str:
    if not fids:
        return default_fid
    best = None
    best_record = None
    for fid in fids:
        f = feat_idx[fid]
        item = (
            base.kind_priority(str(f['kind'])),
            1 if f.get('primary_chr') else 0,
            int(f.get('source_seq_len') or 0),
            fid,
        )
        if best is None or item > best[0]:
            best = (item, fid)
    return best[1] if best else default_fid


def is_reference_kind(kind: object) -> bool:
    return str(kind).startswith('reference_')


def is_assembly_kind(kind: object) -> bool:
    return not is_reference_kind(kind)


def build_label_kind_map(catalog_path: Optional[str], feat_idx: Optional[Dict[str, Dict[str, object]]] = None) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if catalog_path and os.path.exists(catalog_path):
        with open_text(catalog_path, 'rt') as fh:
            reader = csv.DictReader(fh, delimiter='	')
            for row in reader:
                lab = str(row.get('label') or '').strip()
                kind = str(row.get('kind') or '').strip()
                if lab and kind and lab not in out:
                    out[lab] = kind
    if feat_idx is not None:
        for rec in feat_idx.values():
            lab = str(rec.get('label') or '')
            kind = str(rec.get('kind') or '')
            if lab and kind and lab not in out:
                out[lab] = kind
    return out


def fid_mid0(fid: str) -> int:
    p = base.parse_fid(fid)
    return int(round((int(p['start0']) + int(p['end0'])) / 2.0))


def nodeints_compact(nodeints: Sequence[Sequence[int]]) -> str:
    parts = []
    for node_id, s, e in nodeints:
        parts.append(f"{int(node_id)}:{int(s)}-{int(e)}")
    return ','.join(parts)


def member_ids_compact(member_ids: Sequence[str]) -> str:
    return ';'.join(str(x) for x in member_ids)


def member_labels_compact(labels: Sequence[str]) -> str:
    return ';'.join(str(x) for x in labels)


def member_graph_positions_compact(member_ids: Sequence[str], feat_idx: Dict[str, Dict[str, object]]) -> str:
    parts = []
    for fid in member_ids:
        f = feat_idx[fid]
        parts.append(f"{fid}|{nodeints_compact(f.get('nodeints', []))}")
    return ';'.join(parts)


def member_asm_midpoints_compact(member_ids: Sequence[str]) -> str:
    return ';'.join(f"{fid}|{fid_mid0(fid)}" for fid in member_ids)


def member_primary_midpoints_compact(member_ids: Sequence[str], feat_idx: Dict[str, Dict[str, object]]) -> str:
    parts = []
    for fid in member_ids:
        f = feat_idx[fid]
        if f.get('primary_chr'):
            mid = int(round((int(f['primary_start1']) + int(f['primary_end1'])) / 2.0))
            parts.append(f"{fid}|{f['primary_chr']}:{mid}")
    return ';'.join(parts)


def ref_interval1_from_feat(feat: Dict[str, object]) -> Optional[Tuple[str, int, int]]:
    chrom = str(feat.get('primary_chr') or '').strip()
    start1 = feat.get('primary_start1')
    end1 = feat.get('primary_end1')
    try:
        if chrom and start1 not in ('', None, '.', 'nan') and end1 not in ('', None, '.', 'nan'):
            return chrom, int(start1), int(end1)
    except Exception:
        pass
    if str(feat.get('kind') or '') == 'reference_primary':
        try:
            pf = base.parse_fid(str(feat['fid']))
            return str(pf['contig']), int(pf['start0']) + 1, int(pf['end0'])
        except Exception:
            return None
    return None


def build_ref_interval_index(primary_refs: Dict[str, Dict[str, object]]) -> Dict[str, List[Tuple[int, int, str]]]:
    out: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for fid, feat in primary_refs.items():
        iv = ref_interval1_from_feat(feat)
        if iv is None:
            continue
        chrom, start1, end1 = iv
        out[str(chrom)].append((int(start1), int(end1), str(fid)))
    for chrom in list(out.keys()):
        out[chrom].sort(key=lambda x: (int(x[0]), int(x[1]), str(x[2])))
    return out


def interval_overlap_bp1(a_start1: int, a_end1: int, b_start1: int, b_end1: int) -> int:
    s = max(int(a_start1), int(b_start1))
    e = min(int(a_end1), int(b_end1))
    return max(0, e - s + 1)


def query_ref_intervals_overlap(ref_interval_index: Dict[str, List[Tuple[int, int, str]]], chrom: str, start1: int, end1: int) -> List[Tuple[int, int, str]]:
    entries = ref_interval_index.get(str(chrom), [])
    out: List[Tuple[int, int, str]] = []
    if not entries:
        return out
    starts = [x[0] for x in entries]
    i = max(0, bisect.bisect_right(starts, int(end1)))
    j = i - 1
    while j >= 0:
        s, e, fid = entries[j]
        if interval_overlap_bp1(start1, end1, s, e) > 0:
            out.append((s, e, fid))
        j -= 1
    out.reverse()
    return out


def interval_gap_bp1(a_start1: int, a_end1: int, b_start1: int, b_end1: int) -> int:
    if int(a_end1) < int(b_start1):
        return int(b_start1) - int(a_end1) - 1
    if int(b_end1) < int(a_start1):
        return int(a_start1) - int(b_end1) - 1
    return 0


def nonref_site_interval_from_row(row: Dict[str, object], args: argparse.Namespace) -> Optional[Tuple[str, int, int]]:
    seed = nonref_site_seed_from_row(row, args)
    if seed is None:
        return None
    return str(seed['chrom']), int(seed['start1']), int(seed['end1'])


def cluster_nonref_site_rows(assign_rows: List[Dict[str, object]], args: argparse.Namespace) -> None:
    site_window = int(getattr(args, 'nonref_site_window_bp', 100) or 0)
    by_chrom: Dict[str, List[Tuple[int, int, int, int]]] = defaultdict(list)
    for i, row in enumerate(assign_rows):
        if str(row.get('anchor_type')) != 'nonref':
            continue
        seed = nonref_site_seed_from_row(row, args)
        iv = None if seed is None else (str(seed['chrom']), int(seed['start1']), int(seed['end1']))
        if iv is None:
            continue
        chrom, s, e = iv
        by_chrom[str(chrom)].append((int(s), int(e), i, int(row.get('seed_locus_id') or 0) if str(row.get('seed_locus_id') or '').isdigit() else i))

    parent: Dict[int, int] = {}
    def find(x: int) -> int:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda x: (x[0], x[1], x[3], x[2]))
        for _, _, idx, _ in rows:
            parent.setdefault(idx, idx)
        n = len(rows)
        for a in range(n):
            s1, e1, i1, _k1 = rows[a]
            for b in range(a + 1, n):
                s2, e2, i2, _k2 = rows[b]
                if s2 > e1 + site_window:
                    break
                if interval_gap_bp1(s1, e1, s2, e2) <= site_window:
                    union(i1, i2)

    clusters: Dict[Tuple[str, int], List[int]] = defaultdict(list)
    for chrom, rows in by_chrom.items():
        for _s, _e, idx, _sk in rows:
            clusters[(chrom, find(idx))].append(idx)

    ordered_clusters = []
    for (chrom, _root), idxs in clusters.items():
        ivs = [nonref_site_interval_from_row(assign_rows[i], args) for i in idxs]
        ivs = [x for x in ivs if x is not None]
        if not ivs:
            continue
        low = min(int(x[1]) for x in ivs)
        high = max(int(x[2]) for x in ivs)
        ordered_clusters.append((str(chrom), low, high, sorted(idxs)))
    ordered_clusters.sort(key=lambda x: (x[0], x[1], x[2], x[3][0]))

    for cluster_ord, (chrom, low, high, idxs) in enumerate(ordered_clusters, start=1):
        anchor_key = f'NRSITE::{chrom}::{cluster_ord}'
        cluster_n = len(idxs)
        source_counts = Counter(str((nonref_site_seed_from_row(assign_rows[i], args) or {}).get('source') or '') for i in idxs)
        cluster_source = source_counts.most_common(1)[0][0] if source_counts else ''
        for i in idxs:
            row = assign_rows[i]
            row['anchor_key'] = anchor_key
            row['anchor_assign_rule'] = 'nonref_projected_site_cluster' if cluster_n > 1 else 'nonref_projected_site_single'
            row['anchor_confidence'] = 'nonref_projected_site'
            row['nonref_site_cluster_n'] = cluster_n
            row['nonref_site_cluster_ord'] = cluster_ord
            row['nonref_site_chrom'] = chrom
            row['nonref_site_start1'] = low
            row['nonref_site_end1'] = high
            row['nonref_site_source'] = cluster_source
            row['nonref_site_source'] = cluster_source


_ANCHOR_FEAT_IDX: Optional[Dict[str, Dict[str, object]]] = None
_ANCHOR_SEED_LOCUS_BY_FID: Optional[Dict[str, str]] = None
_ANCHOR_PRIMARY_REFS: Optional[Dict[str, Dict[str, object]]] = None
_ANCHOR_REF_NODE_INDEX: Optional[Dict[int, set]] = None
_ANCHOR_REF_INTERVAL_INDEX: Optional[Dict[str, List[Tuple[int, int, str]]]] = None
_ANCHOR_SH_WEIGHTS: Optional[Dict[str, float]] = None
_ANCHOR_ARGS: Optional[argparse.Namespace] = None


def _anchor_assign_one_fid(fid: str) -> Dict[str, object]:
    global _ANCHOR_FEAT_IDX, _ANCHOR_SEED_LOCUS_BY_FID, _ANCHOR_PRIMARY_REFS, _ANCHOR_REF_NODE_INDEX, _ANCHOR_REF_INTERVAL_INDEX, _ANCHOR_SH_WEIGHTS, _ANCHOR_ARGS
    if _ANCHOR_FEAT_IDX is None or _ANCHOR_SEED_LOCUS_BY_FID is None or _ANCHOR_PRIMARY_REFS is None or _ANCHOR_REF_NODE_INDEX is None or _ANCHOR_REF_INTERVAL_INDEX is None or _ANCHOR_ARGS is None:
        raise RuntimeError('anchor worker globals not initialized')
    feat = _ANCHOR_FEAT_IDX[fid]
    seed_locus_id = _ANCHOR_SEED_LOCUS_BY_FID.get(fid, '')
    sv_fields = feature_sv_member_fields(feat)
    hal_fields = feature_hal_member_fields(feat)

    def blank_insert_info() -> Dict[str, object]:
        return {
            'insert_anchor_flag': 0,
            'insert_anchor_proj_type': '',
            'insert_anchor_source': '',
            'insert_site1': '',
            'insert_interval_start1': '',
            'insert_interval_end1': '',
            'insert_left_boundary1': '',
            'insert_right_boundary1': '',
            'multi_insert_n': 0,
            'multi_insert_types': '',
            'multi_insert_sites': '',
            'multi_insert_intervals': '',
            'multi_insert_bp_est': '',
            'multi_insert_longest_bp': '',
        }

    def sv_insert_info_from_feat() -> Dict[str, object]:
        if int(feat.get('sv_ins_n', 0) or 0) <= 0 or not str(feat.get('sv_ins_longest_chrom') or ''):
            return blank_insert_info()
        return {
            'insert_anchor_flag': 1,
            'insert_anchor_proj_type': 'sv_site',
            'insert_anchor_source': 'provided_sv_callset',
            'insert_site1': feat.get('sv_ins_longest_site1', ''),
            'insert_interval_start1': feat.get('sv_ins_longest_start1', ''),
            'insert_interval_end1': feat.get('sv_ins_longest_end1', ''),
            'insert_left_boundary1': '',
            'insert_right_boundary1': '',
            'multi_insert_n': int(feat.get('sv_ins_n', 0) or 0),
            'multi_insert_types': 'INS' if int(feat.get('sv_ins_n', 0) or 0) > 0 else '',
            'multi_insert_sites': feat.get('sv_ins_primary_sites', ''),
            'multi_insert_intervals': feat.get('sv_ins_primary_intervals', ''),
            'multi_insert_bp_est': feat.get('sv_ins_overlap_bp', ''),
            'multi_insert_longest_bp': feat.get('sv_ins_longest_len', ''),
        }

    def blank_nr_proj_info() -> Dict[str, object]:
        return {
            'nr_proj_type': '',
            'nr_proj_source': '',
            'nr_proj_chrom': '',
            'nr_site1': '',
            'nr_interval_start1': '',
            'nr_interval_end1': '',
            'nr_left_boundary1': '',
            'nr_right_boundary1': '',
            'nonref_site_cluster_n': 0,
            'nonref_site_cluster_ord': '',
            'nonref_site_chrom': '',
            'nonref_site_start1': '',
            'nonref_site_end1': '',
            'nonref_site_source': '',
        }

    if str(feat['kind']) == 'reference_primary':
        return {
            'fid': fid,
            'label': feat['label'],
            'kind': feat['kind'],
            'seed_locus_id': seed_locus_id,
            'anchor_key': f'REF::{fid}',
            'anchor_type': 'primary_ref',
            'anchor_ref_fid': fid,
            'anchor_assign_rule': 'self_primary_ref',
            'anchor_confidence': 'self_ref',
            'ref_overlap_n': 1,
            'best_ref_fid': fid,
            'best_ref_shared_bp': int(feat['_bp']),
            'best_ref_ro_graph': 1.0,
            'best_ref_size_ratio': 1.0,
            'best_ref_ctx_graph': 1.0,
            'second_ref_fid': '',
            'second_ref_shared_bp': '',
            'second_ref_ro_graph': '',
            'second_ref_size_ratio': '',
            'second_ref_ctx_graph': '',
            'is_anchor_ref': 1,
            **blank_insert_info(),
            **blank_nr_proj_info(),
            **sv_fields,
            **hal_fields,
        }

    hal_ref_cand = choose_hal_ref_candidate(
        feat, _ANCHOR_REF_INTERVAL_INDEX, _ANCHOR_PRIMARY_REFS
    )
    sv_ref_cand = choose_sv_site_ref_candidate(
        feat, _ANCHOR_REF_INTERVAL_INDEX, _ANCHOR_PRIMARY_REFS
    )

    anchor_key = ''
    anchor_type = 'nonref'
    anchor_ref_fid = ''
    rule = 'no_primary_ref_overlap'
    anchor_confidence = 'nonref'
    ref_overlap_n = 0
    best_ref_fid = ''
    best_ref_shared_bp = ''
    best_ref_ro_graph = ''
    best_ref_size_ratio = ''
    best_ref_ctx_graph = ''
    second_ref_fid = ''
    second_ref_shared_bp = ''
    second_ref_ro_graph = ''
    second_ref_size_ratio = ''
    second_ref_ctx_graph = ''
    insert_info = blank_insert_info()
    nr_proj_info = blank_nr_proj_info()

    if sv_ref_cand is not None:
        best_ref_fid = str(sv_ref_cand['best_ref_fid'])
        second_ref_fid = str(sv_ref_cand.get('second_ref_fid', ''))
        anchor_key = f'REF::{best_ref_fid}'
        anchor_type = 'primary_ref'
        anchor_ref_fid = best_ref_fid
        ref_overlap_n = int(sv_ref_cand.get('candidate_n', 0) or 0)
        best_ref_shared_bp = 1
        rule = 'sv_insert_ref_overlap' if ref_overlap_n == 1 else 'sv_insert_multi_ref_best'
        anchor_confidence = 'sv_insert_ref' if ref_overlap_n == 1 else 'sv_insert_ref_ambiguous'
        insert_info = {
            'insert_anchor_flag': 1,
            'insert_anchor_proj_type': 'sv_site',
            'insert_anchor_source': 'provided_sv_callset',
            'insert_site1': sv_ref_cand.get('insert_site1', ''),
            'insert_interval_start1': sv_ref_cand.get('insert_interval_start1', ''),
            'insert_interval_end1': sv_ref_cand.get('insert_interval_end1', ''),
            'insert_left_boundary1': '',
            'insert_right_boundary1': '',
            'multi_insert_n': int(feat.get('sv_ins_n', 0) or 0),
            'multi_insert_types': 'INS',
            'multi_insert_sites': feat.get('sv_ins_primary_sites', ''),
            'multi_insert_intervals': feat.get('sv_ins_primary_intervals', ''),
            'multi_insert_bp_est': feat.get('sv_ins_overlap_bp', ''),
            'multi_insert_longest_bp': feat.get('sv_ins_longest_len', ''),
        }
    elif hal_ref_cand is not None:
        best_ref_fid = str(hal_ref_cand['best_ref_fid'])
        second_ref_fid = str(hal_ref_cand.get('second_ref_fid', ''))
        anchor_key = f'REF::{best_ref_fid}'
        anchor_type = 'primary_ref'
        anchor_ref_fid = best_ref_fid
        ref_overlap_n = int(hal_ref_cand.get('candidate_n', 0) or 0)
        best_ref_shared_bp = int(hal_ref_cand.get('best_overlap_bp', 0) or 0)
        best_ref_ro_graph = f"{float(hal_ref_cand.get('best_ro', 0.0) or 0.0):.6f}"
        best_ref_size_ratio = f"{float(hal_ref_cand.get('best_size_ratio', 0.0) or 0.0):.6f}"
        if second_ref_fid:
            second_ref_shared_bp = hal_ref_cand.get('second_overlap_bp', '')
            second_ref_ro_graph = (
                '' if hal_ref_cand.get('second_ro', '') in ('', None)
                else f"{float(hal_ref_cand.get('second_ro', 0.0) or 0.0):.6f}"
            )
            second_ref_size_ratio = (
                '' if hal_ref_cand.get('second_size_ratio', '') in ('', None)
                else f"{float(hal_ref_cand.get('second_size_ratio', 0.0) or 0.0):.6f}"
            )
        rule = 'hal_ref_overlap' if ref_overlap_n == 1 else 'hal_multi_ref_best'
        anchor_confidence = 'hal_ref' if ref_overlap_n == 1 else 'hal_ref_ambiguous'
        if int(feat.get('sv_ins_n', 0) or 0) > 0:
            insert_info = sv_insert_info_from_feat()
    else:
        anchor_key = f'NR::{seed_locus_id or fid}'
        if int(feat.get('sv_ins_n', 0) or 0) > 0:
            insert_info = sv_insert_info_from_feat()

    return {
        'fid': fid,
        'label': feat['label'],
        'kind': feat['kind'],
        'seed_locus_id': seed_locus_id,
        'anchor_key': anchor_key,
        'anchor_type': anchor_type,
        'anchor_ref_fid': anchor_ref_fid,
        'anchor_assign_rule': rule,
        'reference_ins_id': sv_ref_cand['reference_ins_id'] if sv_ref_cand is not None else '',
        'all_ins_matched_ref_ids_json': json.dumps(sv_ref_cand['all_matched_ref_ids']) if sv_ref_cand is not None else '[]',
        'ref_overlap_n': ref_overlap_n,
        'best_ref_fid': best_ref_fid,
        'best_ref_shared_bp': best_ref_shared_bp,
        'best_ref_ro_graph': best_ref_ro_graph,
        'best_ref_size_ratio': best_ref_size_ratio,
        'best_ref_ctx_graph': best_ref_ctx_graph,
        'second_ref_fid': second_ref_fid,
        'second_ref_shared_bp': second_ref_shared_bp,
        'second_ref_ro_graph': second_ref_ro_graph,
        'second_ref_size_ratio': second_ref_size_ratio,
        'second_ref_ctx_graph': second_ref_ctx_graph,
        'anchor_confidence': anchor_confidence,
        'is_anchor_ref': 0,
        **insert_info,
        **nr_proj_info,
        **sv_fields,
        **hal_fields,
    }


def cmd_cluster_loci_prod(args: argparse.Namespace) -> None:
    n_features, df = collect_shingle_df(args.features, args.anchor_k, args.max_mid_anchors)
    sh_weights = shingle_weights(n_features, df)

    mem_fh, mem_w = write_tsv_header(args.out_members, [
        'locus_id', 'fid', 'label', 'kind', 'is_lead', 'is_display', 'ambiguous', 'second_best_locus', 'score_delta',
        'match_rule', 'shared_bp', 'ro_graph', 'size_ratio', 'ctx_graph', 'boundary_match', 'weak_seq',
        'primary_chr', 'primary_start1', 'primary_end1', 'nonref_only'
    ])

    loci: Dict[int, Dict[str, object]] = {}
    anchor_index: Dict[str, set] = defaultdict(set)
    high_freq_anchors: set = set()
    next_locus_id = 1

    def index_rep(locus_id: int, rep: Dict[str, object]) -> None:
        for a in rep['_anchors']:
            if a in high_freq_anchors:
                continue
            bucket = anchor_index[a]
            bucket.add(locus_id)
            if len(bucket) > args.max_anchor_bucket:
                high_freq_anchors.add(a)
                anchor_index.pop(a, None)

    def make_locus_from_feat(feat: Dict[str, object]) -> Dict[str, object]:
        return {
            'lead_fid': feat['fid'],
            'rep': feat,
            'display_fid': feat['fid'],
            'member_n': 0,
            'ambiguous_member_n': 0,
            'sample_labels': set(),
            'has_primary_ref_member': 1 if feat['kind'] == 'reference_primary' else 0,
            'has_comparison_ref_member': 1 if feat['kind'] == 'reference_comparison' else 0,
            'has_any_primary_annot': 1 if feat.get('primary_chr') else 0,
        }

    for raw_rec in iter_feature_jsonl(args.features):
        feat = prep_feature(raw_rec, anchor_k=args.anchor_k, max_mid_anchors=args.max_mid_anchors)
        candidate_ids: set = set()
        for a in feat['_anchors']:
            if a in high_freq_anchors:
                continue
            bucket = anchor_index.get(a)
            if bucket:
                candidate_ids.update(bucket)

        best = None
        second = None
        for lid in candidate_ids:
            locus = loci[lid]
            if feat['label'] in locus['sample_labels']:
                continue
            cmp = compare_feature_to_rep(feat, locus['rep'], args=args, sh_weights=sh_weights)
            if cmp['rule'] is None:
                continue
            item = (cmp['composite'], lid, cmp)
            if best is None or item[0] > best[0]:
                second = best
                best = item
            elif second is None or item[0] > second[0]:
                second = item

        ambiguous = 0
        second_best_locus = ''
        score_delta = ''
        if best is None:
            lid = next_locus_id
            next_locus_id += 1
            loci[lid] = make_locus_from_feat(feat)
            index_rep(lid, feat)
            cmp = {
                'rule': 'seed',
                'shared_bp': feat['_bp'],
                'ro_graph': 1.0,
                'size_ratio': 1.0,
                'ctx_graph': 1.0,
                'boundary_match': 1,
                'weak_seq': 1.0,
            }
        else:
            lid = best[1]
            cmp = best[2]
            if second is not None:
                delta = float(best[0]) - float(second[0])
                score_delta = f'{delta:.6f}'
                if delta < args.ambiguity_margin:
                    ambiguous = 1
                    second_best_locus = str(second[1])
                    loci[lid]['ambiguous_member_n'] += 1
            locus = loci[lid]
            if base.kind_priority(str(feat['kind'])) > base.kind_priority(str(locus['rep']['kind'])):
                locus['display_fid'] = feat['fid']
            elif (feat.get('primary_chr') and not locus['rep'].get('primary_chr')):
                locus['display_fid'] = feat['fid']
            if feat['kind'] == 'reference_primary':
                locus['has_primary_ref_member'] = 1
            if feat['kind'] == 'reference_comparison':
                locus['has_comparison_ref_member'] = 1
            if feat.get('primary_chr'):
                locus['has_any_primary_annot'] = 1

        loci[lid]['member_n'] += 1
        loci[lid]['sample_labels'].add(feat['label'])
        mem_w.writerow([
            lid,
            feat['fid'],
            feat['label'],
            feat['kind'],
            1 if feat['fid'] == loci[lid]['lead_fid'] else 0,
            1 if feat['fid'] == loci[lid]['display_fid'] else 0,
            ambiguous,
            second_best_locus,
            score_delta,
            cmp['rule'],
            cmp['shared_bp'],
            f"{cmp['ro_graph']:.6f}",
            f"{cmp['size_ratio']:.6f}",
            f"{cmp['ctx_graph']:.6f}",
            cmp['boundary_match'],
            '' if cmp['weak_seq'] is None else f"{cmp['weak_seq']:.6f}",
            feat.get('primary_chr', ''),
            feat.get('primary_start1', ''),
            feat.get('primary_end1', ''),
            feat.get('nonref_only', ''),
        ])
    mem_fh.close()

    loc_fh, loc_w = write_tsv_header(args.out_locus, [
        'locus_id', 'lead_fid', 'rep_fid', 'display_fid', 'member_n', 'hap_n', 'ambiguous_member_n',
        'has_primary_ref_member', 'has_comparison_ref_member', 'has_any_primary_annot', 'nonref_only_locus',
        'rep_len', 'rep_left_anchor', 'rep_right_anchor', 'rep_nodeints_json'
    ])
    for lid in sorted(loci):
        locus = loci[lid]
        rep = locus['rep']
        loc_w.writerow([
            lid,
            locus['lead_fid'],
            rep['fid'],
            locus['display_fid'],
            locus['member_n'],
            len(locus['sample_labels']),
            locus['ambiguous_member_n'],
            locus['has_primary_ref_member'],
            locus['has_comparison_ref_member'],
            locus['has_any_primary_annot'],
            1 if locus['has_any_primary_annot'] == 0 else 0,
            rep['_bp'],
            '|'.join(rep['_left']),
            '|'.join(rep['_right']),
            json_compact(rep['nodeints']),
        ])
    loc_fh.close()

    print(json_pretty({'n_loci': len(loci), 'n_features': n_features}), file=sys.stderr)


def load_prepped_feature_index(features_path: str, anchor_k: int, max_mid_anchors: int) -> Dict[str, Dict[str, object]]:
    idx: Dict[str, Dict[str, object]] = {}
    for rec in iter_feature_jsonl(features_path):
        p = prep_feature(rec, anchor_k=anchor_k, max_mid_anchors=max_mid_anchors)
        idx[str(p['fid'])] = p
    return idx


def choose_medoid_fid(member_fids: List[str], feat_idx: Dict[str, Dict[str, object]], sh_weights: Optional[Dict[str, float]], max_pairwise: int) -> str:
    fids = list(member_fids)
    if len(fids) == 1:
        return fids[0]

    if len(fids) > max_pairwise:
        refs = [fid for fid in fids if str(feat_idx[fid]['kind']).startswith('reference_')]
        others = [fid for fid in fids if fid not in refs]
        others.sort(key=lambda x: int(feat_idx[x].get('source_seq_len') or 0), reverse=True)
        keep = refs + others[:max(0, max_pairwise - len(refs))]
        keep = list(dict.fromkeys(keep))
        fids = keep
    best = None
    total_eps = 1e-12
    for fid in fids:
        score_terms = [10.0]
        for other in fids:
            if other == fid:
                continue
            sim = pairwise_graph_similarity(feat_idx[fid], feat_idx[other], sh_weights)
            score_terms.append(float(sim['score']))
        total = math.fsum(score_terms)
        tie_item = (
            base.kind_priority(str(feat_idx[fid]['kind'])),
            1 if feat_idx[fid].get('primary_chr') else 0,
            int(feat_idx[fid].get('source_seq_len') or 0),
            fid,
        )
        if best is None:
            best = (total, tie_item, fid)
            continue
        best_total, best_tie_item, _best_fid = best
        if total > best_total + total_eps:
            best = (total, tie_item, fid)
            continue
        if math.isclose(total, best_total, rel_tol=0.0, abs_tol=total_eps) and tie_item > best_tie_item:
            best = (total, tie_item, fid)
    return best[2]


_POLISH_FEAT_IDX: Optional[Dict[str, Dict[str, object]]] = None
_POLISH_SH_WEIGHTS: Optional[Dict[str, float]] = None
_POLISH_MAX_PAIRWISE: int = 256


_POLISH_STORE_READER: Optional["PolishSQLiteFeatureStore"] = None
_POLISH_STORE_READ_CHUNK: int = 900


def _strip_suffix_once(path: str, suffix: str) -> str:
    return path[:-len(suffix)] if path.endswith(suffix) else path


def default_polish_store_paths(features_path: str, db_override: str = '') -> Tuple[str, str, str]:
    if db_override:
        db_path = os.path.abspath(db_override)
        base = _strip_suffix_once(db_path, '.sqlite')
    else:
        base = os.path.abspath(features_path)
        for suffix in ('.gz', '.jsonl', '.tsv'):
            base = _strip_suffix_once(base, suffix)
        base = base + '.step04_polish'
        db_path = base + '.sqlite'
    shdf_path = base + '.shdf.gz'
    meta_path = base + '.meta.json'
    return db_path, shdf_path, meta_path


def _polish_store_expected_meta(features_path: str, anchor_k: int, max_mid_anchors: int) -> Dict[str, object]:
    st = os.stat(features_path)
    return {
        'format': 3,
        'features_path': os.path.abspath(features_path),
        'features_size': int(st.st_size),
        'features_mtime': float(st.st_mtime),
        'anchor_k': int(anchor_k),
        'max_mid_anchors': int(max_mid_anchors),
    }


def _polish_store_meta_matches(meta_path: str, expected: Dict[str, object]) -> bool:
    try:
        with open(meta_path, 'rt') as fh:
            meta = json.load(fh)
    except Exception:
        return False
    for k, v in expected.items():
        if meta.get(k) != v:
            return False
    return True


def _polish_store_pack_record(p: Dict[str, object]) -> bytes:
    payload = (
        p.get('label'),
        p.get('kind'),
        p.get('source_seq_len'),
        p.get('primary_chr'),
        p.get('primary_start1'),
        p.get('primary_end1'),
        p.get('nonref_only'),
        p.get('source_seq'),
        p.get('nodeints'),
        tuple(p.get('_left', ()) or ()),
        tuple(p.get('_right', ()) or ()),
        tuple(p.get('_anchors', ()) or ()),
        tuple(p.get('_shingles', ()) or ()),
        int(p.get('_bp') or 0),
    )
    return marshal.dumps(payload)


def _polish_store_unpack_record(fid: str, blob: bytes) -> Dict[str, object]:
    (
        label,
        kind,
        source_seq_len,
        primary_chr,
        primary_start1,
        primary_end1,
        nonref_only,
        source_seq,
        nodeints,
        left,
        right,
        anchors,
        shingles,
        bp,
    ) = marshal.loads(blob)
    node_map: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for node_id, s, e in nodeints:
        node_map[int(node_id)].append((int(s), int(e)))
    return {
        'fid': str(fid),
        'label': label,
        'kind': kind,
        'source_seq_len': source_seq_len,
        'primary_chr': primary_chr,
        'primary_start1': primary_start1,
        'primary_end1': primary_end1,
        'nonref_only': nonref_only,
        'source_seq': source_seq,
        'nodeints': nodeints,
        '_left': tuple(left),
        '_right': tuple(right),
        '_anchors': list(anchors),
        '_shingles': set(shingles),
        '_node_map': node_map,
        '_bp': int(bp),
    }


class PolishSQLiteFeatureStore:
    def __init__(self, db_path: str, *, readonly: bool = True):
        self.db_path = os.path.abspath(db_path)
        self.conn = sqlite3.connect(self.db_path, timeout=300, isolation_level=None, check_same_thread=False)
        self.conn.execute('PRAGMA temp_store=MEMORY')
        self.conn.execute('PRAGMA mmap_size=268435456')
        self.conn.execute('PRAGMA cache_size=-262144')
        if readonly:
            self.conn.execute('PRAGMA query_only=ON')

    def close(self) -> None:
        try:
            self.conn.close()
        except Exception:
            pass

    def __enter__(self) -> "PolishSQLiteFeatureStore":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def get(self, fid: str) -> Dict[str, object]:
        row = self.conn.execute('SELECT blob FROM polish_feature_store WHERE fid=?', (str(fid),)).fetchone()
        if row is None:
            raise KeyError(fid)
        return _polish_store_unpack_record(str(fid), row[0])

    def get_many(self, fids: Sequence[str], *, chunk_size: int = 900) -> Dict[str, Dict[str, object]]:
        out: Dict[str, Dict[str, object]] = {}
        if not fids:
            return out
        step = max(1, int(chunk_size or 900))
        for i in range(0, len(fids), step):
            part = [str(x) for x in fids[i:i + step]]
            sql = 'SELECT fid, blob FROM polish_feature_store WHERE fid IN (%s)' % (','.join('?' for _ in part))
            for fid, blob in self.conn.execute(sql, part):
                fid_s = str(fid)
                out[fid_s] = _polish_store_unpack_record(fid_s, blob)
        if len(out) != len(fids):
            missing = [str(fid) for fid in fids if str(fid) not in out]
            raise KeyError(f'Could not fetch {len(missing)} feature(s) from polish feature store; first missing: {missing[:3]}')
        return out

    def iter_sorted(self) -> Iterator[Dict[str, object]]:
        for fid, blob in self.conn.execute('SELECT fid, blob FROM polish_feature_store ORDER BY fid'):
            yield _polish_store_unpack_record(str(fid), blob)


def build_polish_feature_store(
    features_path: str,
    db_path: str,
    shdf_path: str,
    meta_path: str,
    *,
    anchor_k: int,
    max_mid_anchors: int,
    batch_size: int,
) -> Tuple[int, Dict[str, int]]:
    ensure_parent(db_path)
    for path in (db_path, shdf_path, meta_path):
        if path and os.path.exists(path):
            os.remove(path)
    conn = sqlite3.connect(db_path, timeout=300, isolation_level=None)
    try:
        conn.execute('PRAGMA journal_mode=OFF')
        conn.execute('PRAGMA synchronous=OFF')
        conn.execute('PRAGMA locking_mode=EXCLUSIVE')
        conn.execute('PRAGMA temp_store=MEMORY')
        conn.execute('PRAGMA mmap_size=268435456')
        conn.execute('PRAGMA cache_size=-262144')
        conn.execute('CREATE TABLE polish_feature_store (fid TEXT PRIMARY KEY, blob BLOB) WITHOUT ROWID')
        cur = conn.cursor()
        conn.execute('BEGIN')
        df: Counter = Counter()
        batch: List[Tuple[str, bytes]] = []
        n_features = 0
        for rec in iter_feature_jsonl(features_path):
            p = prep_feature(rec, anchor_k=anchor_k, max_mid_anchors=max_mid_anchors)
            df.update(set(p['_shingles']))
            batch.append((str(p['fid']), sqlite3.Binary(_polish_store_pack_record(p))))
            if len(batch) >= max(1, int(batch_size)):
                cur.executemany('INSERT INTO polish_feature_store(fid, blob) VALUES (?, ?)', batch)
                batch.clear()
            n_features += 1
            if n_features % 1000000 == 0:
                print(f'[polish-store] cached {n_features} features', file=sys.stderr)
        if batch:
            cur.executemany('INSERT INTO polish_feature_store(fid, blob) VALUES (?, ?)', batch)
        conn.commit()
        with gzip.open(shdf_path, 'wb') as fh:
            marshal.dump((int(n_features), dict(df)), fh)
        meta = _polish_store_expected_meta(features_path, anchor_k=anchor_k, max_mid_anchors=max_mid_anchors)
        meta['record_count'] = int(n_features)
        with open(meta_path, 'wt') as fh:
            json.dump(meta, fh, sort_keys=True)
        return int(n_features), dict(df)
    except Exception:
        try:
            conn.close()
        finally:
            for path in (db_path, shdf_path, meta_path):
                if path and os.path.exists(path):
                    try:
                        os.remove(path)
                    except Exception:
                        pass
        raise
    finally:
        try:
            conn.close()
        except Exception:
            pass


def load_or_build_polish_feature_store(args: argparse.Namespace) -> Tuple[str, int, Dict[str, int]]:
    if str(args.features).endswith('.sqlite'):
        n, df = collect_shingle_df(args.features, args.anchor_k, args.max_mid_anchors)
        return args.features, n, df
    db_path, shdf_path, meta_path = default_polish_store_paths(args.features, getattr(args, 'feature_store_db', '') or '')
    expected = _polish_store_expected_meta(args.features, anchor_k=args.anchor_k, max_mid_anchors=args.max_mid_anchors)
    reuse_ok = (not bool(getattr(args, 'feature_store_rebuild', False))) and os.path.exists(db_path) and os.path.exists(shdf_path) and os.path.exists(meta_path) and _polish_store_meta_matches(meta_path, expected)
    if reuse_ok:
        print(f'[polish-store] reuse {db_path}', file=sys.stderr)
        with gzip.open(shdf_path, 'rb') as fh:
            n_features, df = marshal.load(fh)
        return db_path, int(n_features), dict(df)
    print(f'[polish-store] build {db_path}', file=sys.stderr)
    n_features, df = build_polish_feature_store(
        args.features,
        db_path,
        shdf_path,
        meta_path,
        anchor_k=int(args.anchor_k),
        max_mid_anchors=int(args.max_mid_anchors),
        batch_size=int(getattr(args, 'feature_store_batch_size', 4096) or 4096),
    )
    return db_path, int(n_features), dict(df)


def read_members_by_locus_stream(path: str) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = defaultdict(list)
    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='	')
        for row in reader:
            out[str(row['locus_id'])].append(str(row['fid']))
    return out


def _init_polish_store_worker(db_path: str, read_chunk: int) -> None:
    global _POLISH_STORE_READER, _POLISH_STORE_READ_CHUNK
    _POLISH_STORE_READ_CHUNK = max(1, int(read_chunk or 900))
    _POLISH_STORE_READER = PolishSQLiteFeatureStore(db_path, readonly=True)


def _choose_medoid_store_task(task: Tuple[str, List[str]]) -> Tuple[str, str, str, int, int]:
    global _POLISH_STORE_READER, _POLISH_SH_WEIGHTS, _POLISH_MAX_PAIRWISE, _POLISH_STORE_READ_CHUNK
    if _POLISH_STORE_READER is None:
        raise RuntimeError('polish sqlite worker globals not initialized')
    locus_id, fids = task
    feat_map = _POLISH_STORE_READER.get_many(fids, chunk_size=_POLISH_STORE_READ_CHUNK)
    medoid = choose_medoid_fid(fids, feat_map, _POLISH_SH_WEIGHTS, max_pairwise=_POLISH_MAX_PAIRWISE)
    display = pick_display_candidate(fids, feat_map, medoid)
    display_feat = feat_map[display]
    return locus_id, medoid, display, base.kind_priority(str(display_feat['kind'])), 1 if display_feat.get('primary_chr') else 0


def cmd_polish_loci_sqlite_exact(args: argparse.Namespace) -> None:
    db_path, n_features, df = load_or_build_polish_feature_store(args)
    sh_weights = shingle_weights(n_features, df)
    del df
    gc.collect()

    prov_loc = read_tsv(args.in_locus)
    prov_loc['locus_id'] = prov_loc['locus_id'].astype(str)
    lead_by_locus = dict(zip(prov_loc['locus_id'], prov_loc['lead_fid']))
    del prov_loc
    gc.collect()

    members_by_locus = read_members_by_locus_stream(args.in_members)
    medoid_tasks = [(locus_id, fids) for locus_id, fids in members_by_locus.items()]

    medoid_info: Dict[str, Dict[str, object]] = {}
    global _POLISH_SH_WEIGHTS, _POLISH_MAX_PAIRWISE
    _POLISH_SH_WEIGHTS = sh_weights
    _POLISH_MAX_PAIRWISE = int(args.max_medoid_members)
    read_chunk = max(1, int(getattr(args, 'feature_store_read_chunk', 900) or 900))

    with PolishSQLiteFeatureStore(db_path, readonly=True) as store:
        if int(getattr(args, 'threads', 1)) > 1:
            ctx = mp.get_context('fork')
            chunksize = max(1, int(getattr(args, 'chunksize', 1)))
            with ctx.Pool(
                processes=int(args.threads),
                maxtasksperchild=int(getattr(args, 'maxtasksperchild', 0) or 0) or None,
                initializer=_init_polish_store_worker,
                initargs=(db_path, read_chunk),
            ) as pool:
                for locus_id, medoid, display, display_kind_priority, display_has_primary in pool.imap(_choose_medoid_store_task, medoid_tasks, chunksize=chunksize):
                    medoid_info[locus_id] = {
                        'medoid_fid': medoid,
                        'orig_lead_fid': lead_by_locus.get(locus_id, medoid),
                        'display_fid': display,
                        'display_kind_priority': int(display_kind_priority),
                        'display_has_primary': int(display_has_primary),
                    }
        else:
            for locus_id, fids in medoid_tasks:
                feat_map = store.get_many(fids, chunk_size=read_chunk)
                medoid = choose_medoid_fid(fids, feat_map, sh_weights, max_pairwise=args.max_medoid_members)
                display = pick_display_candidate(fids, feat_map, medoid)
                display_feat = feat_map[display]
                medoid_info[locus_id] = {
                    'medoid_fid': medoid,
                    'orig_lead_fid': lead_by_locus.get(locus_id, medoid),
                    'display_fid': display,
                    'display_kind_priority': int(base.kind_priority(str(display_feat['kind']))),
                    'display_has_primary': 1 if display_feat.get('primary_chr') else 0,
                }

        rep_loci: Dict[str, Dict[str, object]] = {}
        for locus_id, info in medoid_info.items():
            medoid_fid = str(info['medoid_fid'])
            rep = store.get(medoid_fid)
            rep_loci[locus_id] = {
                'seed_locus_id': locus_id,
                'lead_fid': str(info['orig_lead_fid']),
                'medoid_fid': medoid_fid,
                'rep': rep,
                'display_fid': str(info['display_fid']),
                '_display_kind_priority': int(info['display_kind_priority']),
                '_display_has_primary': int(info['display_has_primary']),
                'member_n': 0,
                'ambiguous_member_n': 0,
                'sample_labels': set(),
                'has_primary_ref_member': 1 if rep['kind'] == 'reference_primary' else 0,
                'has_comparison_ref_member': 1 if rep['kind'] == 'reference_comparison' else 0,
                'has_any_primary_annot': 1 if rep.get('primary_chr') else 0,
            }

        del members_by_locus
        del medoid_tasks
        del medoid_info
        gc.collect()

        anchor_index: Dict[str, set] = defaultdict(set)
        high_freq_anchors: set = set()
        for locus_id, locus in rep_loci.items():
            rep = locus['rep']
            for a in rep['_anchors']:
                if a in high_freq_anchors:
                    continue
                bucket = anchor_index[a]
                bucket.add(locus_id)
                if len(bucket) > args.max_anchor_bucket:
                    high_freq_anchors.add(a)
                    anchor_index.pop(a, None)

        out_mem_fh, out_mem_w = write_tsv_header(args.out_members, [
            'locus_id', 'seed_locus_id', 'fid', 'label', 'kind', 'is_lead', 'is_medoid', 'is_display', 'ambiguous', 'second_best_locus', 'score_delta',
            'match_rule', 'shared_bp', 'ro_graph', 'size_ratio', 'ctx_graph', 'boundary_match', 'weak_seq',
            'primary_chr', 'primary_start1', 'primary_end1', 'nonref_only'
        ])

        used_ids = {int(x) for x in rep_loci.keys()}
        next_locus_id = max(used_ids) + 1 if used_ids else 1

        def new_singleton_locus(feat: Dict[str, object]) -> Tuple[str, Dict[str, object]]:
            nonlocal next_locus_id
            lid = str(next_locus_id)
            next_locus_id += 1
            locus = {
                'seed_locus_id': lid,
                'lead_fid': feat['fid'],
                'medoid_fid': feat['fid'],
                'rep': feat,
                'display_fid': feat['fid'],
                '_display_kind_priority': int(base.kind_priority(str(feat['kind']))),
                '_display_has_primary': 1 if feat.get('primary_chr') else 0,
                'member_n': 0,
                'ambiguous_member_n': 0,
                'sample_labels': set(),
                'has_primary_ref_member': 1 if feat['kind'] == 'reference_primary' else 0,
                'has_comparison_ref_member': 1 if feat['kind'] == 'reference_comparison' else 0,
                'has_any_primary_annot': 1 if feat.get('primary_chr') else 0,
            }
            rep_loci[lid] = locus
            for a in feat['_anchors']:
                if a in high_freq_anchors:
                    continue
                bucket = anchor_index[a]
                bucket.add(lid)
                if len(bucket) > args.max_anchor_bucket:
                    high_freq_anchors.add(a)
                    anchor_index.pop(a, None)
            return lid, locus

        for feat in store.iter_sorted():
            candidate_ids: set = set()
            for a in feat['_anchors']:
                if a in high_freq_anchors:
                    continue
                bucket = anchor_index.get(a)
                if bucket:
                    candidate_ids.update(bucket)

            best = None
            second = None
            for lid in candidate_ids:
                locus = rep_loci[lid]
                if feat['label'] in locus['sample_labels']:
                    continue
                cmp = compare_feature_to_rep(feat, locus['rep'], args=args, sh_weights=sh_weights)
                if cmp['rule'] is None:
                    continue
                item = (cmp['composite'], lid, cmp)
                if best is None or item[0] > best[0]:
                    second = best
                    best = item
                elif second is None or item[0] > second[0]:
                    second = item

            ambiguous = 0
            second_best_locus = ''
            score_delta = ''
            if best is None:
                lid, locus = new_singleton_locus(feat)
                cmp = {
                    'rule': 'medoid_seed',
                    'shared_bp': feat['_bp'],
                    'ro_graph': 1.0,
                    'size_ratio': 1.0,
                    'ctx_graph': 1.0,
                    'boundary_match': 1,
                    'weak_seq': 1.0,
                }
            else:
                lid = best[1]
                locus = rep_loci[lid]
                cmp = best[2]
                if second is not None:
                    delta = float(best[0]) - float(second[0])
                    score_delta = f'{delta:.6f}'
                    if delta < args.ambiguity_margin:
                        ambiguous = 1
                        second_best_locus = str(second[1])
                        locus['ambiguous_member_n'] += 1
                feat_kind_priority = int(base.kind_priority(str(feat['kind'])))
                if feat_kind_priority > int(locus['_display_kind_priority']):
                    locus['display_fid'] = feat['fid']
                    locus['_display_kind_priority'] = feat_kind_priority
                    locus['_display_has_primary'] = 1 if feat.get('primary_chr') else 0
                elif feat.get('primary_chr') and not int(locus['_display_has_primary']):
                    locus['display_fid'] = feat['fid']
                    locus['_display_kind_priority'] = feat_kind_priority
                    locus['_display_has_primary'] = 1
                if feat['kind'] == 'reference_primary':
                    locus['has_primary_ref_member'] = 1
                if feat['kind'] == 'reference_comparison':
                    locus['has_comparison_ref_member'] = 1
                if feat.get('primary_chr'):
                    locus['has_any_primary_annot'] = 1

            locus['member_n'] += 1
            locus['sample_labels'].add(feat['label'])
            out_mem_w.writerow([
                lid,
                locus['seed_locus_id'],
                feat['fid'],
                feat['label'],
                feat['kind'],
                1 if feat['fid'] == locus['lead_fid'] else 0,
                1 if feat['fid'] == locus['medoid_fid'] else 0,
                1 if feat['fid'] == locus['display_fid'] else 0,
                ambiguous,
                second_best_locus,
                score_delta,
                cmp['rule'],
                cmp['shared_bp'],
                f"{cmp['ro_graph']:.6f}",
                f"{cmp['size_ratio']:.6f}",
                f"{cmp['ctx_graph']:.6f}",
                cmp['boundary_match'],
                '' if cmp['weak_seq'] is None else f"{cmp['weak_seq']:.6f}",
                feat.get('primary_chr', ''),
                feat.get('primary_start1', ''),
                feat.get('primary_end1', ''),
                feat.get('nonref_only', ''),
            ])
        out_mem_fh.close()

        out_loc_fh, out_loc_w = write_tsv_header(args.out_locus, [
            'locus_id', 'seed_locus_id', 'lead_fid', 'medoid_fid', 'display_fid',
            'member_n', 'hap_n', 'ambiguous_member_n',
            'has_primary_ref_member', 'has_comparison_ref_member', 'has_any_primary_annot', 'nonref_only_locus',
            'medoid_len', 'medoid_left_anchor', 'medoid_right_anchor', 'medoid_nodeints_json'
        ])
        for lid in sorted(rep_loci, key=lambda x: int(x)):
            locus = rep_loci[lid]
            rep = locus['rep']
            out_loc_w.writerow([
                lid,
                locus['seed_locus_id'],
                locus['lead_fid'],
                locus['medoid_fid'],
                locus['display_fid'],
                locus['member_n'],
                len(locus['sample_labels']),
                locus['ambiguous_member_n'],
                locus['has_primary_ref_member'],
                locus['has_comparison_ref_member'],
                locus['has_any_primary_annot'],
                1 if locus['has_any_primary_annot'] == 0 else 0,
                rep['_bp'],
                '|'.join(rep['_left']),
                '|'.join(rep['_right']),
                json_compact(rep['nodeints']),
            ])
        out_loc_fh.close()

    print(json_pretty({'n_polished_loci': len(rep_loci)}), file=sys.stderr)


def cmd_polish_loci(args: argparse.Namespace) -> None:
    return cmd_polish_loci_sqlite_exact(args)


def cmd_anchor_loci_primary(args: argparse.Namespace) -> None:
    from pancgi_anchor_bounded import run
    from pancgi_anchor_disk import execute
    run(args, execute)


def tokens_to_shingles(tokens: Sequence[str], k: int, max_mid_anchors: int = 8) -> set:
    return set('M:' + s for s in shingle_strings(tokens, k, max_mid_anchors=max_mid_anchors))


def determine_orientation_to_medoid(feat: Dict[str, object], medoid: Dict[str, object], anchor_k: int, sh_weights: Optional[Dict[str, float]]) -> Tuple[int, float, float]:
    fwd = weighted_jaccard(tokens_to_shingles(feat['_tokens'], anchor_k), set(medoid['_shingles']), sh_weights)
    rev = weighted_jaccard(tokens_to_shingles(reverse_tokens(feat['_tokens']), anchor_k), set(medoid['_shingles']), sh_weights)
    if rev > fwd:
        return -1, fwd, rev
    return 1, fwd, rev


def choose_seq_medoid(member_fids: List[str], seq_map: Dict[str, Optional[str]], feat_idx: Dict[str, Dict[str, object]], args: argparse.Namespace) -> str:
    if len(member_fids) == 1:
        return member_fids[0]
    fids = list(member_fids)
    max_n = int(getattr(args, 'max_allele_medoid_members', 0) or 0)
    if max_n > 0 and len(fids) > max_n:
        refs = [fid for fid in fids if str(feat_idx[fid]['kind']).startswith('reference_')]
        others = [fid for fid in fids if fid not in refs]
        others.sort(key=lambda x: (int(feat_idx[x].get('source_seq_len') or 0), base.kind_priority(str(feat_idx[x]['kind'])), x), reverse=True)
        keep = refs + others[:max(0, max_n - len(refs))]
        fids = list(dict.fromkeys(keep))
    best = None
    for fid in fids:
        total = 0.0
        for other in fids:
            if other == fid:
                total += 10.0
                continue
            sim = seq_similarity(seq_map[fid], seq_map[other], args=args)
            total += (sim if sim is not None else 0.0)
        item = (
            total,
            base.kind_priority(str(feat_idx[fid]['kind'])),
            1 if feat_idx[fid].get('primary_chr') else 0,
            int(feat_idx[fid].get('source_seq_len') or 0),
            fid,
        )
        if best is None or item > best[0]:
            best = (item, fid)
    return best[1]


def choose_allele_rep(member_fids: List[str], seq_map: Dict[str, Optional[str]], feat_idx: Dict[str, Dict[str, object]], args: argparse.Namespace) -> str:
    strategy = str(args.allele_rep_strategy)
    asm_only = [fid for fid in member_fids if is_assembly_kind(feat_idx[fid]['kind'])]
    if strategy == 'asm_seq_medoid':
        if asm_only:
            return choose_seq_medoid(asm_only, seq_map, feat_idx, args)
        return choose_seq_medoid(member_fids, seq_map, feat_idx, args)
    if strategy == 'seq_medoid':
        return choose_seq_medoid(member_fids, seq_map, feat_idx, args)
    if strategy == 'ref_first':
        refs = [fid for fid in member_fids if is_reference_kind(feat_idx[fid]['kind'])]
        if refs:
            return choose_seq_medoid(refs, seq_map, feat_idx, args)
        return choose_seq_medoid(member_fids, seq_map, feat_idx, args)
    if strategy == 'asm_longest':
        pool = asm_only if asm_only else member_fids
        return max(pool, key=lambda x: (int(feat_idx[x].get('source_seq_len') or 0), x))
    if strategy == 'longest':
        return max(member_fids, key=lambda x: (int(feat_idx[x].get('source_seq_len') or 0), x))
    return member_fids[0]


def pair_len_ratio(fid1: str, fid2: str, feat_idx: Dict[str, Dict[str, object]]) -> float:
    len1 = int(feat_idx[fid1].get('source_seq_len') or 0)
    len2 = int(feat_idx[fid2].get('source_seq_len') or 0)
    return min(len1, len2) / max(1, max(len1, len2))


def seq_pair_similarity_cached(fid1: str, fid2: str, seq_map: Dict[str, Optional[str]], args: argparse.Namespace, cache: Dict[Tuple[str, str], Optional[float]]) -> Optional[float]:
    key = (fid1, fid2) if fid1 <= fid2 else (fid2, fid1)
    if key not in cache:
        cache[key] = seq_similarity(seq_map.get(fid1), seq_map.get(fid2), args=args)
    return cache[key]


def pair_passes_allele_threshold(fid1: str, fid2: str, seq_map: Dict[str, Optional[str]], feat_idx: Dict[str, Dict[str, object]], args: argparse.Namespace, cache: Dict[Tuple[str, str], Optional[float]]) -> Tuple[bool, Optional[float], float]:
    lr = pair_len_ratio(fid1, fid2, feat_idx)
    if lr < float(args.min_len_ratio):
        return False, None, lr
    sim = seq_pair_similarity_cached(fid1, fid2, seq_map, args, cache)
    if sim is None:
        return False, None, lr
    return bool(sim >= float(args.identity)), sim, lr


def greedy_rep_seed_clusters(members_sorted: Sequence[str], norm_seq: Dict[str, Optional[str]], feat_idx: Dict[str, Dict[str, object]], args: argparse.Namespace, cache: Dict[Tuple[str, str], Optional[float]]) -> List[List[str]]:
    provisional: List[Dict[str, object]] = []
    for fid in members_sorted:
        best = None
        for idx, a in enumerate(provisional):
            rep_fid = str(a['seed_rep_fid'])
            ok, sim, len_ratio = pair_passes_allele_threshold(fid, rep_fid, norm_seq, feat_idx, args, cache)
            if not ok:
                continue
            item = (float(sim), idx, float(len_ratio))
            if best is None or item > best:
                best = item
        if best is None:
            provisional.append({'seed_rep_fid': fid, 'members': [fid]})
        else:
            _, idx, _ = best
            provisional[idx]['members'].append(fid)
    return [list(x['members']) for x in provisional]


def clique_refine_clusters(seed_clusters: Sequence[Sequence[str]], global_order: Sequence[str], norm_seq: Dict[str, Optional[str]], feat_idx: Dict[str, Dict[str, object]], args: argparse.Namespace, cache: Dict[Tuple[str, str], Optional[float]]) -> List[List[str]]:
    order_rank = {fid: i for i, fid in enumerate(global_order)}
    refined: List[List[str]] = []
    for members in seed_clusters:
        ordered = sorted(list(members), key=lambda fid: order_rank.get(fid, 10**12))
        subclusters: List[List[str]] = []
        for fid in ordered:
            best = None
            for idx, cluster in enumerate(subclusters):
                sims = []
                ok_all = True
                min_len_ratio = 1.0
                for other in cluster:
                    ok, sim, len_ratio = pair_passes_allele_threshold(fid, other, norm_seq, feat_idx, args, cache)
                    min_len_ratio = min(min_len_ratio, float(len_ratio))
                    if not ok:
                        ok_all = False
                        break
                    sims.append(float(sim))
                if not ok_all:
                    continue
                score = (statistics.mean(sims) if sims else 10.0, min_len_ratio, -idx)
                if best is None or score > best[0]:
                    best = (score, idx)
            if best is None:
                subclusters.append([fid])
            else:
                _, idx = best
                subclusters[idx].append(fid)
        refined.extend(subclusters)
    return refined


def cluster_members_with_mode(members_sorted: Sequence[str], norm_seq: Dict[str, Optional[str]], feat_idx: Dict[str, Dict[str, object]], args: argparse.Namespace) -> Tuple[List[List[str]], Dict[Tuple[str, str], Optional[float]]]:
    cache: Dict[Tuple[str, str], Optional[float]] = {}
    mode = str(getattr(args, 'allele_cluster_mode', 'allpairs_clique'))
    if mode == 'allpairs_clique':
        seed_clusters = [list(members_sorted)]
        clusters = clique_refine_clusters(seed_clusters, members_sorted, norm_seq, feat_idx, args, cache)
    elif mode == 'greedy':
        clusters = greedy_rep_seed_clusters(members_sorted, norm_seq, feat_idx, args, cache)
    else:
        seed_clusters = greedy_rep_seed_clusters(members_sorted, norm_seq, feat_idx, args, cache)
        clusters = clique_refine_clusters(seed_clusters, members_sorted, norm_seq, feat_idx, args, cache)
    return clusters, cache


def build_center_star_alignment(records: Sequence[Tuple[str, str]], args: argparse.Namespace) -> Tuple[List[Tuple[str, str]], Dict[str, object]]:
    if len(records) < 2:
        return list(records), {'n_seq': len(records), 'aligned_len': len(records[0][1]) if records else 0}
    pairwise2 = load_pairwise2()
    center_idx = max(range(len(records)), key=lambda i: (len(records[i][1]), -i))
    center_head, center_seq = records[center_idx]
    profiles: List[Tuple[str, List[str], List[str]]] = []
    max_gaps = [0] * (len(center_seq) + 1)
    for head, seq in records:
        if seq == center_seq:
            ins = [''] * (len(center_seq) + 1)
            bases = list(seq)
        else:
            aln = pairwise2.align.globalms(
                center_seq,
                seq,
                float(getattr(args, 'msa_match', 2.0)),
                float(getattr(args, 'msa_mismatch', -3.0)),
                float(getattr(args, 'msa_gap_open', -5.0)),
                float(getattr(args, 'msa_gap_extend', -2.0)),
                one_alignment_only=True,
            )
            if not aln:
                raise RuntimeError('center-star pairwise2 returned no alignment')
            center_aln, seq_aln, _score, _start, _end = aln[0]
            ins = [''] * (len(center_seq) + 1)
            bases = []
            raw_idx = 0
            for ca, sa in zip(center_aln, seq_aln):
                if ca == '-':
                    ins[raw_idx] += sa
                    continue
                bases.append(sa)
                raw_idx += 1
            if len(bases) != len(center_seq):
                raise RuntimeError(f'center-star aligned base count mismatch: {len(bases)} vs {len(center_seq)}')
        profiles.append((head, ins, bases))
        for pos, val in enumerate(ins):
            if len(val) > max_gaps[pos]:
                max_gaps[pos] = len(val)
    center_parts: List[str] = []
    for pos in range(len(center_seq) + 1):
        if max_gaps[pos] > 0:
            center_parts.append('-' * max_gaps[pos])
        if pos < len(center_seq):
            center_parts.append(center_seq[pos])
    center_aln = ''.join(center_parts)
    out: List[Tuple[str, str]] = []
    for head, ins, bases in profiles:
        parts: List[str] = []
        for pos in range(len(center_seq) + 1):
            cur = ins[pos]
            if cur:
                parts.append(cur)
            if len(cur) < max_gaps[pos]:
                parts.append('-' * (max_gaps[pos] - len(cur)))
            if pos < len(center_seq):
                parts.append(bases[pos])
        seq_aln = ''.join(parts)
        if len(seq_aln) != len(center_aln):
            raise RuntimeError(f'center-star aligned length mismatch for {head}: {len(seq_aln)} vs {len(center_aln)}')
        out.append((head, seq_aln))
    meta = {'n_seq': len(records), 'center_header': center_head, 'center_index': center_idx, 'center_len': len(center_seq), 'aligned_len': len(center_aln)}
    return out, meta


def run_mafft_alignment(records: Sequence[Tuple[str, str]], args: argparse.Namespace) -> List[Tuple[str, str]]:
    mafft_bin = str(getattr(args, 'msa_mafft_bin', 'mafft') or 'mafft')
    with tempfile.TemporaryDirectory(prefix='cpgi_mafft_') as td:
        in_fa = os.path.join(td, 'input.fa')
        out_fa = os.path.join(td, 'output.fa')
        with open(in_fa, 'wt') as fh:
            for head, seq in records:
                write_fasta_record(fh, head, seq)
        cmd = [mafft_bin, '--anysymbol', '--quiet', '--thread', '1', in_fa]
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        with open(out_fa, 'wt') as fh:
            fh.write(proc.stdout)
        return read_simple_fasta(out_fa)


def read_simple_fasta(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    head = None
    buf: List[str] = []
    with open(path, 'rt') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if head is not None:
                    records.append((head, ''.join(buf)))
                head = line[1:]
                buf = []
            else:
                buf.append(line)
    if head is not None:
        records.append((head, ''.join(buf)))
    return records


def maybe_write_msa(locus_dir: str, records: Sequence[Tuple[str, str]], args: argparse.Namespace) -> Optional[str]:
    backend = str(getattr(args, 'dump_msa_backend', 'none') or 'none').lower()
    if backend == 'none':
        return None
    min_members = int(getattr(args, 'dump_msa_min_members', 2) or 2)
    max_members = int(getattr(args, 'dump_msa_max_members', 0) or 0)
    if len(records) < min_members:
        return None
    if max_members > 0 and len(records) > max_members:
        return None
    out_name = str(getattr(args, 'dump_msa_name', 'center_star_msa.fa') or 'center_star_msa.fa')
    out_path = os.path.join(locus_dir, out_name)
    if backend == 'center_star':
        aligned, _meta = build_center_star_alignment(records, args)
    elif backend == 'mafft':
        aligned = run_mafft_alignment(records, args)
    else:
        raise ValueError(f'Unknown dump_msa_backend: {backend}')
    with open(out_path, 'wt') as fh:
        for head, seq in aligned:
            write_fasta_record(fh, head, seq)
    return out_path


def dump_locus_cluster_artifacts(
    locus_id: str,
    members: Sequence[str],
    norm_seq: Dict[str, Optional[str]],
    feat_idx: Dict[str, Dict[str, object]],
    allele_rows: Sequence[Dict[str, object]],
    allele_member_rows: Sequence[Dict[str, object]],
    args: argparse.Namespace,
) -> None:
    outdir = str(getattr(args, 'dump_locus_dir', '') or '')
    if not outdir:
        return
    min_members = int(getattr(args, 'dump_min_members', 0) or 0)
    if min_members > 0 and len(members) < min_members:
        return
    locus_dir = os.path.join(outdir, f'locus_{locus_id}')
    os.makedirs(locus_dir, exist_ok=True)
    msa_records: List[Tuple[str, str]] = []
    with open(os.path.join(locus_dir, 'normalized.fa'), 'wt') as fh:
        for fid in members:
            seq = norm_seq.get(fid)
            if not seq:
                continue
            head = f"{fid}|label={feat_idx[fid]['label']}|kind={feat_idx[fid]['kind']}"
            msa_records.append((head, str(seq)))
            write_fasta_record(fh, head, str(seq))
    with open(os.path.join(locus_dir, 'allele_clusters.tsv'), 'wt') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['allele_id', 'fid', 'label', 'kind', 'is_allele_rep', 'orientation_norm', 'seq_similarity', 'len_ratio'])
        for mr in allele_member_rows:
            w.writerow([
                mr['allele_id'], mr['fid'], mr['label'], mr['kind'], mr['is_allele_rep'], mr['orientation_norm'],
                '' if mr['seq_similarity'] is None else f"{float(mr['seq_similarity']):.6f}",
                f"{float(mr['len_ratio']):.6f}",
            ])
    if bool(getattr(args, 'dump_pairwise_matrix', False)):
        with gzip.open(os.path.join(locus_dir, 'pairwise_similarity.tsv.gz'), 'wt') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['fid1', 'fid2', 'seq_similarity'])
            mm = list(members)
            for i, fid1 in enumerate(mm):
                for j in range(i, len(mm)):
                    fid2 = mm[j]
                    sim = seq_similarity(norm_seq.get(fid1), norm_seq.get(fid2), args=args)
                    w.writerow([fid1, fid2, '' if sim is None else f"{float(sim):.6f}"])


ALLELE_CATALOG_HEADER = [
    'allele_id', 'locus_id', 'allele_rep_fid', 'display_fid', 'display_kind', 'display_label',
    'major_allele', 'member_n', 'hap_n', 'asm_hap_n', 'ref_hap_n', 'asm_member_n', 'ref_member_n',
    'has_primary_ref_member', 'has_comparison_ref_member', 'primary_ref_member_fids', 'comparison_ref_member_fids',
    'rep_seq_len', 'rep_len', 'rep_cpg_n', 'rep_gc_n', 'rep_pct_gc', 'rep_oe', 'allele_rep_strategy', 'allele_cluster_mode',
    'allele_member_ids', 'allele_member_labels'
]
ALLELE_MEMBER_HEADER = [
    'allele_id', 'locus_id', 'fid', 'label', 'kind', 'is_allele_rep', 'orientation_norm', 'ctx_fwd', 'ctx_rev',
    'seq_similarity', 'len_ratio'
]


def locus_sort_key(locus_id: str):
    s = str(locus_id)
    return (0, int(s)) if s.isdigit() else (1, s)


def read_locus_catalog_medoid(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'locus_id' not in (reader.fieldnames or []):
            raise ValueError(f'locus catalog missing locus_id column: {path}')
        for row in reader:
            locus_id = str(row.get('locus_id') or '')
            if not locus_id:
                continue
            medoid = str(row.get('medoid_fid') or row.get('lead_fid') or row.get('display_fid') or '')
            if not medoid:
                raise ValueError(f'locus catalog row missing medoid_fid/lead_fid for locus {locus_id}')
            out[locus_id] = medoid
    return out


def count_members_by_locus_stream(path: str, selected_loci: Optional[set] = None) -> Dict[str, int]:
    out: Dict[str, int] = defaultdict(int)
    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'locus_id' not in (reader.fieldnames or []):
            raise ValueError(f'locus members table missing locus_id column: {path}')
        for row in reader:
            locus_id = str(row.get('locus_id') or '')
            if not locus_id:
                continue
            if selected_loci is not None and locus_id not in selected_loci:
                continue
            out[locus_id] += 1
    return out


def read_locus_ids_file(path: str) -> List[str]:
    ids: List[str] = []
    seen = set()
    with open_text(path, 'rt') as fh:
        first = True
        locus_idx: Optional[int] = None
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if first:
                first = False
                if 'locus_id' in fields:
                    locus_idx = fields.index('locus_id')
                    continue
            val = fields[locus_idx if locus_idx is not None and locus_idx < len(fields) else 0].strip()
            if not val or val == 'locus_id':
                continue
            if val not in seen:
                ids.append(val)
                seen.add(val)
    return ids


def parse_locus_ids_arg(value: str) -> List[str]:
    out: List[str] = []
    seen = set()
    for part in str(value or '').replace(',', '\n').splitlines():
        x = part.strip()
        if x and x not in seen:
            out.append(x)
            seen.add(x)
    return out


def select_locus_ids_from_args(all_locus_ids: Sequence[str], args: argparse.Namespace) -> Tuple[List[str], Dict[str, object]]:
    selected = list(all_locus_ids)
    info: Dict[str, object] = {'total_loci': len(all_locus_ids)}
    explicit: List[str] = []
    if getattr(args, 'locus_ids_file', ''):
        explicit.extend(read_locus_ids_file(str(args.locus_ids_file)))
    if getattr(args, 'locus_ids', ''):
        explicit.extend(parse_locus_ids_arg(str(args.locus_ids)))
    if explicit:
        allowed = set(explicit)
        selected = [x for x in selected if x in allowed]
        info['explicit_loci'] = len(allowed)

    start = getattr(args, 'locus_start_index', None)
    end = getattr(args, 'locus_end_index', None)
    if start is not None or end is not None:
        s0 = 0 if start is None else int(start)
        e0 = len(selected) if end is None else int(end)
        if s0 < 0 or e0 < s0:
            raise ValueError(f'Invalid locus index range: start={start}, end={end}')
        selected = selected[s0:e0]
        info['range_start_index'] = s0
        info['range_end_index'] = e0

    shard_count = int(getattr(args, 'locus_shard_count', 0) or 0)
    if shard_count > 0:
        shard_index = int(getattr(args, 'locus_shard_index', 0) or 0)
        if shard_index < 0 or shard_index >= shard_count:
            raise ValueError(f'--locus-shard-index must be 0 <= index < --locus-shard-count; got {shard_index}/{shard_count}')
        mode = str(getattr(args, 'locus_shard_mode', 'contiguous') or 'contiguous')
        if mode == 'modulo':
            selected = [lid for i, lid in enumerate(selected) if (i % shard_count) == shard_index]
        elif mode == 'contiguous':
            n = len(selected)
            lo = (n * shard_index) // shard_count
            hi = (n * (shard_index + 1)) // shard_count
            selected = selected[lo:hi]
            info['contiguous_shard_start'] = lo
            info['contiguous_shard_end'] = hi
        else:
            raise ValueError(f'Unknown --locus-shard-mode: {mode}')
        info['shard_index'] = shard_index
        info['shard_count'] = shard_count
        info['shard_mode'] = mode
    info['selected_loci'] = len(selected)
    return selected, info


def read_members_for_loci_stream(path: str, selected_loci: Optional[set] = None) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = defaultdict(list)
    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'locus_id' not in (reader.fieldnames or []) or 'fid' not in (reader.fieldnames or []):
            raise ValueError(f'locus members table must contain locus_id and fid columns: {path}')
        for row in reader:
            locus_id = str(row.get('locus_id') or '')
            if not locus_id:
                continue
            if selected_loci is not None and locus_id not in selected_loci:
                continue
            out[locus_id].append(str(row['fid']))
    return out


def build_allele_tasks_from_args(args: argparse.Namespace) -> Tuple[List[Tuple[str, List[str], str]], Dict[str, object]]:
    medoid_by_locus = read_locus_catalog_medoid(args.locus_catalog)
    all_locus_ids = sorted(medoid_by_locus, key=locus_sort_key)
    selected_locus_ids, info = select_locus_ids_from_args(all_locus_ids, args)
    selected_set = set(selected_locus_ids)
    read_all_members = len(selected_locus_ids) == len(all_locus_ids) and not cluster_allele_uses_locus_filter(args)
    members_by_locus = read_members_for_loci_stream(args.locus_members, None if read_all_members else selected_set)
    tasks: List[Tuple[str, List[str], str]] = []
    missing_medoid = []
    empty_loci = []
    for locus_id in selected_locus_ids:
        fids = members_by_locus.get(locus_id, [])
        if not fids:
            empty_loci.append(locus_id)
            continue
        medoid = str(medoid_by_locus.get(locus_id) or '')
        if not medoid:
            missing_medoid.append(locus_id)
            continue
        tasks.append((locus_id, fids, medoid))
    if missing_medoid:
        raise ValueError(f'{len(missing_medoid)} selected loci missing medoid_fid; first: {missing_medoid[:3]}')
    info['task_loci'] = len(tasks)
    info['empty_selected_loci'] = len(empty_loci)
    return tasks, info


def write_locus_task_list(tasks: Sequence[Tuple[str, List[str], str]], path: str) -> None:
    fh, w = write_tsv_header(path, ['global_selected_index', 'locus_id', 'medoid_fid', 'member_n'])
    try:
        for i, (locus_id, members, medoid) in enumerate(tasks):
            w.writerow([i, locus_id, medoid, len(members)])
    finally:
        fh.close()


def load_prepped_feature_index_selected(features_path: str, needed_fids: set, anchor_k: int, max_mid_anchors: int) -> Dict[str, Dict[str, object]]:
    if str(features_path).endswith('.sqlite'):
        from pancgi_features import selected_raw
        collect_shingle_df(features_path, anchor_k, max_mid_anchors)
        return {fid: prep_feature(rec, anchor_k, max_mid_anchors) for fid, rec in selected_raw(features_path, needed_fids).items()}
    needed = set(map(str, needed_fids))
    idx: Dict[str, Dict[str, object]] = {}
    if not needed:
        return idx
    for rec in iter_feature_jsonl(features_path):
        fid = str(rec.get('fid') or '')
        if fid not in needed:
            continue
        idx[fid] = prep_feature(rec, anchor_k=anchor_k, max_mid_anchors=max_mid_anchors)
        if len(idx) == len(needed):
            break
    if len(idx) != len(needed):
        missing = sorted(needed - set(idx))
        raise KeyError(f'Could not load {len(missing)} selected feature(s) from {features_path}; first missing: {missing[:5]}')
    return idx


def _strip_known_feature_suffixes(path: str) -> str:
    base_path = os.path.abspath(path)
    for suffix in ('.gz', '.jsonl', '.tsv'):
        base_path = _strip_suffix_once(base_path, suffix)
    return base_path


def default_allele_store_paths(features_path: str, db_override: str = '') -> Tuple[str, str, str]:
    if db_override:
        db_path = os.path.abspath(db_override)
        base_path = _strip_suffix_once(db_path, '.sqlite')
    else:
        base_path = _strip_known_feature_suffixes(features_path) + '.step07_allele'
        db_path = base_path + '.sqlite'
    return db_path, base_path + '.shdf.gz', base_path + '.meta.json'


def _allele_store_expected_meta(features_path: str, anchor_k: int, max_mid_anchors: int) -> Dict[str, object]:
    st = os.stat(features_path)
    return {
        'format': 2,
        'purpose': 'step07_allele_feature_store',
        'features_path': os.path.abspath(features_path),
        'features_size': int(st.st_size),
        'features_mtime': float(st.st_mtime),
        'anchor_k': int(anchor_k),
        'max_mid_anchors': int(max_mid_anchors),
    }


def _allele_store_meta_matches(meta_path: str, expected: Dict[str, object]) -> bool:
    try:
        with open(meta_path, 'rt') as fh:
            meta = json.load(fh)
    except Exception:
        return False
    return all(meta.get(k) == v for k, v in expected.items())


def _allele_store_pack_record(p: Dict[str, object]) -> bytes:
    payload = (
        p.get('label'), p.get('kind'), p.get('source_seq_len'), p.get('asm_len'),
        p.get('cpg_n'), p.get('gc_n'), p.get('pct_gc'), p.get('oe'),
        p.get('primary_chr'), p.get('primary_start1'), p.get('primary_end1'),
        p.get('nonref_only'),
        p.get('source_seq'), p.get('nodeints'),
        tuple(p.get('_tokens', ()) or ()), tuple(p.get('_left', ()) or ()), tuple(p.get('_right', ()) or ()),
        tuple(p.get('_anchors', ()) or ()), tuple(p.get('_shingles', ()) or ()), int(p.get('_bp') or 0),
    )
    return marshal.dumps(payload)


def _allele_store_unpack_record(fid: str, blob: bytes) -> Dict[str, object]:
    (
        label, kind, source_seq_len, asm_len,
        cpg_n, gc_n, pct_gc, oe,
        primary_chr, primary_start1, primary_end1,
        nonref_only,
        source_seq, nodeints,
        tokens, left, right, anchors, shingles, bp,
    ) = marshal.loads(blob)
    node_map: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for node_id, s, e in (nodeints or []):
        node_map[int(node_id)].append((int(s), int(e)))
    return {
        'fid': str(fid),
        'label': label,
        'kind': kind,
        'source_seq_len': source_seq_len,
        'asm_len': asm_len,
        'cpg_n': cpg_n,
        'gc_n': gc_n,
        'pct_gc': pct_gc,
        'oe': oe,
        'primary_chr': primary_chr,
        'primary_start1': primary_start1,
        'primary_end1': primary_end1,
        'nonref_only': nonref_only,
        'source_seq': source_seq,
        'nodeints': nodeints,
        '_tokens': list(tokens),
        '_left': tuple(left),
        '_right': tuple(right),
        '_anchors': list(anchors),
        '_shingles': set(shingles),
        '_node_map': node_map,
        '_bp': int(bp),
    }


class AlleleSQLiteFeatureStore:
    def __init__(self, db_path: str, *, readonly: bool = True):
        self.db_path = os.path.abspath(db_path)
        self.conn = sqlite3.connect(self.db_path, timeout=300, isolation_level=None, check_same_thread=False)
        self.conn.execute('PRAGMA temp_store=MEMORY')
        self.conn.execute('PRAGMA mmap_size=268435456')
        self.conn.execute('PRAGMA cache_size=-262144')
        if readonly:
            self.conn.execute('PRAGMA query_only=ON')

    def close(self) -> None:
        try:
            self.conn.close()
        except Exception:
            pass

    def __enter__(self) -> "AlleleSQLiteFeatureStore":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def get_many(self, fids: Sequence[str], *, chunk_size: int = 900) -> Dict[str, Dict[str, object]]:
        unique_fids = list(dict.fromkeys(str(x) for x in fids))
        out: Dict[str, Dict[str, object]] = {}
        if not unique_fids:
            return out
        step = max(1, int(chunk_size or 900))
        for i in range(0, len(unique_fids), step):
            part = unique_fids[i:i + step]
            sql = 'SELECT fid, blob FROM allele_feature_store WHERE fid IN (%s)' % ','.join('?' for _ in part)
            for fid, blob in self.conn.execute(sql, part):
                fid_s = str(fid)
                out[fid_s] = _allele_store_unpack_record(fid_s, blob)
        if len(out) != len(unique_fids):
            missing = [fid for fid in unique_fids if fid not in out]
            raise KeyError(f'Could not fetch {len(missing)} allele feature(s) from {self.db_path}; first missing: {missing[:5]}')
        return out


def build_allele_feature_store(
    features_path: str,
    db_path: str,
    shdf_path: str,
    meta_path: str,
    *,
    anchor_k: int,
    max_mid_anchors: int,
    batch_size: int,
) -> Tuple[int, Dict[str, int]]:
    ensure_parent(db_path)
    for path in (db_path, shdf_path, meta_path):
        if path and os.path.exists(path):
            os.remove(path)
    conn = sqlite3.connect(db_path, timeout=300, isolation_level=None)
    try:
        conn.execute('PRAGMA journal_mode=OFF')
        conn.execute('PRAGMA synchronous=OFF')
        conn.execute('PRAGMA locking_mode=EXCLUSIVE')
        conn.execute('PRAGMA temp_store=MEMORY')
        conn.execute('PRAGMA mmap_size=268435456')
        conn.execute('PRAGMA cache_size=-262144')
        conn.execute('CREATE TABLE allele_feature_store (fid TEXT PRIMARY KEY, blob BLOB) WITHOUT ROWID')
        cur = conn.cursor()
        conn.execute('BEGIN')
        df: Counter = Counter()
        batch: List[Tuple[str, bytes]] = []
        n_features = 0
        for rec in iter_feature_jsonl(features_path):
            p = prep_feature(rec, anchor_k=anchor_k, max_mid_anchors=max_mid_anchors)
            df.update(set(p['_shingles']))
            batch.append((str(p['fid']), sqlite3.Binary(_allele_store_pack_record(p))))
            if len(batch) >= max(1, int(batch_size)):
                cur.executemany('INSERT INTO allele_feature_store(fid, blob) VALUES (?, ?)', batch)
                batch.clear()
            n_features += 1
            if n_features % 1000000 == 0:
                print(f'[allele-store] cached {n_features} features', file=sys.stderr)
        if batch:
            cur.executemany('INSERT INTO allele_feature_store(fid, blob) VALUES (?, ?)', batch)
        conn.commit()
        with gzip.open(shdf_path, 'wb') as fh:
            marshal.dump((int(n_features), dict(df)), fh)
        meta = _allele_store_expected_meta(features_path, anchor_k=anchor_k, max_mid_anchors=max_mid_anchors)
        meta['record_count'] = int(n_features)
        with open(meta_path, 'wt') as fh:
            json.dump(meta, fh, sort_keys=True)
        return int(n_features), dict(df)
    except Exception:
        try:
            conn.close()
        finally:
            for path in (db_path, shdf_path, meta_path):
                if path and os.path.exists(path):
                    try:
                        os.remove(path)
                    except Exception:
                        pass
        raise
    finally:
        try:
            conn.close()
        except Exception:
            pass


def load_or_build_allele_feature_store(args: argparse.Namespace) -> Tuple[str, int, Dict[str, int]]:
    db_path, shdf_path, meta_path = default_allele_store_paths(args.features, getattr(args, 'feature_store_db', '') or '')
    expected = _allele_store_expected_meta(args.features, anchor_k=args.anchor_k, max_mid_anchors=args.max_mid_anchors)
    reuse_ok = (not bool(getattr(args, 'feature_store_rebuild', False))) and os.path.exists(db_path) and os.path.exists(shdf_path) and os.path.exists(meta_path) and _allele_store_meta_matches(meta_path, expected)
    if reuse_ok:
        print(f'[allele-store] reuse {db_path}', file=sys.stderr)
        with gzip.open(shdf_path, 'rb') as fh:
            n_features, df = marshal.load(fh)
        return db_path, int(n_features), dict(df)
    print(f'[allele-store] build {db_path}', file=sys.stderr)
    n_features, df = build_allele_feature_store(
        args.features, db_path, shdf_path, meta_path,
        anchor_k=int(args.anchor_k), max_mid_anchors=int(args.max_mid_anchors),
        batch_size=int(getattr(args, 'feature_store_batch_size', 4096) or 4096),
    )
    return db_path, int(n_features), dict(df)


def cluster_allele_uses_locus_filter(args: argparse.Namespace) -> bool:
    return bool(getattr(args, 'locus_ids_file', '') or getattr(args, 'locus_ids', '') or int(getattr(args, 'locus_shard_count', 0) or 0) > 0 or getattr(args, 'locus_start_index', None) is not None or getattr(args, 'locus_end_index', None) is not None)


def load_allele_feature_context(args: argparse.Namespace, tasks: Sequence[Tuple[str, List[str], str]], info: Dict[str, object]) -> Tuple[Dict[str, Dict[str, object]], Dict[str, float], Dict[str, object]]:
    needed_fids: List[str] = []
    for _locus_id, members, medoid in tasks:
        needed_fids.extend(members)
        needed_fids.append(str(medoid))
    needed_unique = list(dict.fromkeys(needed_fids))
    mode = str(getattr(args, 'feature_load_mode', 'auto') or 'auto')
    if mode not in ('auto', 'jsonl', 'sqlite'):
        raise ValueError(f'Unknown --feature-load-mode: {mode}')
    use_sqlite = mode == 'sqlite' or (mode == 'auto' and (cluster_allele_uses_locus_filter(args) or bool(getattr(args, 'feature_store_db', '') or '')))
    out_info = dict(info)
    out_info['needed_features'] = len(needed_unique)
    if str(args.features).endswith('.sqlite'):
        n, df = collect_shingle_df(args.features, args.anchor_k, args.max_mid_anchors)
        feats = load_prepped_feature_index_selected(args.features, set(needed_unique), args.anchor_k, args.max_mid_anchors)
        out_info.update(feature_load_mode_resolved='shared_sqlite', features_in_store=n, feature_store_db=args.features)
        return feats, shingle_weights(n, df), out_info
    out_info['feature_load_mode_resolved'] = 'sqlite' if use_sqlite else 'jsonl'
    if use_sqlite:
        db_path, n_features, df = load_or_build_allele_feature_store(args)
        sh_weights = shingle_weights(n_features, df)
        with AlleleSQLiteFeatureStore(db_path, readonly=True) as store:
            feat_idx = store.get_many(needed_unique, chunk_size=int(getattr(args, 'feature_store_read_chunk', 900) or 900))
        out_info['feature_store_db'] = db_path
        out_info['features_in_store'] = int(n_features)
        return feat_idx, sh_weights, out_info
    n_features, df = collect_shingle_df(args.features, args.anchor_k, args.max_mid_anchors)
    sh_weights = shingle_weights(n_features, df)
    if cluster_allele_uses_locus_filter(args):
        feat_idx = load_prepped_feature_index_selected(args.features, set(needed_unique), anchor_k=args.anchor_k, max_mid_anchors=args.max_mid_anchors)
    else:
        feat_idx = load_prepped_feature_index(args.features, anchor_k=args.anchor_k, max_mid_anchors=args.max_mid_anchors)
    out_info['features_scanned_for_df'] = int(n_features)
    return feat_idx, sh_weights, out_info


def iter_tsv_groups_by_locus(path: str, expected_header: List[str]) -> Iterator[Tuple[str, List[List[str]]]]:
    with open_text(path, 'rt') as fh:
        reader = csv.reader(fh, delimiter='\t')
        try:
            header = next(reader)
        except StopIteration:
            return
        if header != expected_header:
            raise ValueError(f'Header mismatch in shard {path}')
        locus_idx = header.index('locus_id')
        cur_locus = None
        rows: List[List[str]] = []
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f'Malformed row in shard {path}: expected {len(header)} fields, got {len(row)}')
            lid = row[locus_idx]
            if cur_locus is None:
                cur_locus = lid
            if lid != cur_locus:
                yield str(cur_locus), rows
                cur_locus = lid
                rows = []
            rows.append(row)
        if cur_locus is not None:
            yield str(cur_locus), rows


def merge_sharded_tsv_by_locus(shard_paths: Sequence[str], locus_order: Sequence[str], out_path: str, header: List[str], *, require_all_loci: bool = True) -> Dict[str, object]:
    import heapq
    rank = {str(lid): i for i, lid in enumerate(locus_order)}
    out_tmp = temp_output_path(out_path)
    out_fh, out_w = write_tsv_header(out_tmp, header)
    heap: List[Tuple[int, int, str, List[List[str]]]] = []
    iters = []
    seen_loci = set()
    try:
        for idx, path in enumerate(shard_paths):
            it = iter_tsv_groups_by_locus(path, header)
            iters.append(it)
            try:
                lid, rows = next(it)
            except StopIteration:
                continue
            if lid not in rank:
                raise ValueError(f'Shard {path} contains locus not in locus list: {lid}')
            heapq.heappush(heap, (rank[lid], idx, lid, rows))
        n_groups = 0
        n_rows = 0
        last_rank = -1
        while heap:
            r, idx, lid, rows = heapq.heappop(heap)
            if lid in seen_loci:
                raise ValueError(f'Duplicate locus in allele shards: {lid}')
            if r < last_rank:
                raise ValueError(f'Shard merge order violation around locus {lid}')
            last_rank = r
            seen_loci.add(lid)
            for row in rows:
                out_w.writerow(row)
                n_rows += 1
            n_groups += 1
            try:
                next_lid, next_rows = next(iters[idx])
            except StopIteration:
                continue
            if next_lid not in rank:
                raise ValueError(f'Shard {shard_paths[idx]} contains locus not in locus list: {next_lid}')
            heapq.heappush(heap, (rank[next_lid], idx, next_lid, next_rows))
        out_fh.close()
        out_fh = None
        if require_all_loci and len(seen_loci) != len(locus_order):
            missing = [lid for lid in locus_order if lid not in seen_loci][:10]
            raise ValueError(f'Merged shards cover {len(seen_loci)}/{len(locus_order)} loci; first missing: {missing}')
        finalize_output_path(out_tmp, out_path)
        return {'groups_written': n_groups, 'rows_written': n_rows, 'missing_loci': max(0, len(locus_order) - len(seen_loci))}
    except Exception:
        try:
            if out_fh is not None:
                out_fh.close()
        finally:
            if os.path.exists(out_tmp):
                try:
                    os.remove(out_tmp)
                except Exception:
                    pass
        raise

_ALLELE_WORKER_FEAT_IDX: Optional[Dict[str, Dict[str, object]]] = None
_ALLELE_WORKER_SH_WEIGHTS: Optional[Dict[str, float]] = None
_ALLELE_WORKER_ARGS: Optional[argparse.Namespace] = None


def _process_one_locus_alleles(task: Tuple[str, List[str], str]) -> Tuple[str, List[Dict[str, object]], List[Dict[str, object]]]:
    global _ALLELE_WORKER_FEAT_IDX, _ALLELE_WORKER_SH_WEIGHTS, _ALLELE_WORKER_ARGS
    if _ALLELE_WORKER_FEAT_IDX is None or _ALLELE_WORKER_ARGS is None:
        raise RuntimeError('allele worker globals not initialized')
    feat_idx = _ALLELE_WORKER_FEAT_IDX
    sh_weights = _ALLELE_WORKER_SH_WEIGHTS
    args = _ALLELE_WORKER_ARGS
    locus_id, members, medoid_fid = task
    medoid = feat_idx[medoid_fid]

    norm_seq: Dict[str, Optional[str]] = {}
    orient: Dict[str, int] = {}
    ctx_pair: Dict[str, Tuple[float, float]] = {}
    for fid in members:
        f = feat_idx[fid]
        sign, fwd, rev = determine_orientation_to_medoid(f, medoid, anchor_k=args.anchor_k, sh_weights=sh_weights)
        orient[fid] = sign
        ctx_pair[fid] = (fwd, rev)
        seq = f.get('source_seq')
        norm_seq[fid] = reverse_complement(seq) if sign == -1 else seq

    members_sorted = sorted(
        members,
        key=lambda fid: (
            0 if fid == medoid_fid else 1,
            0 if is_assembly_kind(feat_idx[fid]['kind']) else 1,
            -(int(feat_idx[fid].get('source_seq_len') or 0)),
            fid,
        )
    )
    clustered_members, pair_cache = cluster_members_with_mode(members_sorted, norm_seq, feat_idx, args)

    allele_rows: List[Dict[str, object]] = []
    allele_member_rows: List[Dict[str, object]] = []
    for afids in clustered_members:
        afids = list(afids)
        rep_fid = choose_allele_rep(afids, norm_seq, feat_idx, args)
        display_fid = rep_fid
        allele_id = rep_fid
        rep_feat = feat_idx[rep_fid]
        asm_labels = sorted({str(feat_idx[fid]['label']) for fid in afids if is_assembly_kind(feat_idx[fid]['kind'])})
        ref_labels = sorted({str(feat_idx[fid]['label']) for fid in afids if is_reference_kind(feat_idx[fid]['kind'])})
        primary_ref_fids = [fid for fid in afids if str(feat_idx[fid]['kind']) == 'reference_primary']
        comparison_ref_fids = [fid for fid in afids if str(feat_idx[fid]['kind']) == 'reference_comparison']
        asm_member_n = sum(1 for fid in afids if is_assembly_kind(feat_idx[fid]['kind']))
        ref_member_n = len(afids) - asm_member_n
        for fid in afids:
            len_ratio = min(int(feat_idx[fid].get('source_seq_len') or 0), int(rep_feat.get('source_seq_len') or 0)) / max(1, max(int(feat_idx[fid].get('source_seq_len') or 0), int(rep_feat.get('source_seq_len') or 0)))
            sim = seq_pair_similarity_cached(fid, rep_fid, norm_seq, args, pair_cache)
            allele_member_rows.append({
                'allele_id': allele_id,
                'locus_id': locus_id,
                'fid': fid,
                'label': feat_idx[fid]['label'],
                'kind': feat_idx[fid]['kind'],
                'is_allele_rep': 1 if fid == rep_fid else 0,
                'orientation_norm': orient[fid],
                'ctx_fwd': ctx_pair[fid][0],
                'ctx_rev': ctx_pair[fid][1],
                'seq_similarity': sim,
                'len_ratio': len_ratio,
            })
        allele_rows.append({
            'allele_id': allele_id,
            'locus_id': locus_id,
            'allele_rep_fid': rep_fid,
            'display_fid': display_fid,
            'display_kind': feat_idx[display_fid]['kind'],
            'display_label': feat_idx[display_fid]['label'],
            'member_n': len(afids),
            'hap_n': len({feat_idx[fid]['label'] for fid in afids}),
            'asm_hap_n': len(asm_labels),
            'ref_hap_n': len(ref_labels),
            'asm_member_n': asm_member_n,
            'ref_member_n': ref_member_n,
            'has_primary_ref_member': 1 if primary_ref_fids else 0,
            'has_comparison_ref_member': 1 if comparison_ref_fids else 0,
            'primary_ref_member_fids': ';'.join(primary_ref_fids),
            'comparison_ref_member_fids': ';'.join(comparison_ref_fids),
            'rep_seq_len': rep_feat.get('source_seq_len') or rep_feat.get('asm_len') or rep_feat['_bp'],
            'rep_len': rep_feat.get('asm_len') or rep_feat['_bp'],
            'rep_cpg_n': rep_feat.get('cpg_n'),
            'rep_gc_n': rep_feat.get('gc_n'),
            'rep_pct_gc': rep_feat.get('pct_gc'),
            'rep_oe': rep_feat.get('oe'),
            'member_fids': afids,
            'member_labels': sorted({str(feat_idx[fid]['label']) for fid in afids}),
            'allele_cluster_mode': str(getattr(args, 'allele_cluster_mode', 'allpairs_clique')),
        })
    if allele_rows:
        major = max(
            allele_rows,
            key=lambda x: (
                int(x['asm_hap_n']),
                int(x['asm_member_n']),
                int(x['hap_n']),
                int(x['member_n']),
                int(x['rep_seq_len'] or 0),
                str(x['allele_rep_fid']),
            ),
        )
        major_id = str(major['allele_id'])
    else:
        major_id = ''
    for a in allele_rows:
        a['major_allele'] = 1 if a['allele_id'] == major_id else 0
        a['allele_rep_strategy'] = args.allele_rep_strategy
    dump_locus_cluster_artifacts(locus_id, members, norm_seq, feat_idx, allele_rows, allele_member_rows, args)
    return locus_id, allele_rows, allele_member_rows

def cmd_cluster_alleles_prod(args: argparse.Namespace) -> None:
    tasks, task_info = build_allele_tasks_from_args(args)
    if getattr(args, 'write_locus_list', ''):
        write_locus_task_list(tasks, str(args.write_locus_list))
    if bool(getattr(args, 'locus_list_only', False)):
        print(json_pretty({'mode': 'locus_list_only', **task_info, 'write_locus_list': str(getattr(args, 'write_locus_list', '') or '')}), file=sys.stderr)
        return

    out_allele_tmp = temp_output_path(args.out_allele)
    out_allele_members_tmp = temp_output_path(args.out_allele_members)
    allele_fh, allele_w = write_tsv_header(out_allele_tmp, ALLELE_CATALOG_HEADER)
    allele_mem_fh, allele_mem_w = write_tsv_header(out_allele_members_tmp, ALLELE_MEMBER_HEADER)

    feat_idx: Dict[str, Dict[str, object]] = {}
    sh_weights: Dict[str, float] = {}
    feature_info: Dict[str, object] = dict(task_info)
    if tasks:
        feat_idx, sh_weights, feature_info = load_allele_feature_context(args, tasks, task_info)

    global _ALLELE_WORKER_FEAT_IDX, _ALLELE_WORKER_SH_WEIGHTS, _ALLELE_WORKER_ARGS
    _ALLELE_WORKER_FEAT_IDX = feat_idx
    _ALLELE_WORKER_SH_WEIGHTS = sh_weights
    _ALLELE_WORKER_ARGS = args

    n_alleles = 0
    n_loci_new = 0
    flush_every = int(getattr(args, 'flush_every_loci', 0) or 0)

    def maybe_flush() -> None:
        if flush_every > 0 and (n_loci_new % flush_every) == 0:
            try:
                allele_fh.flush()
                allele_mem_fh.flush()
            except Exception:
                pass

    def write_result(result):
        nonlocal n_alleles, n_loci_new
        _locus_id, allele_rows, allele_member_rows = result
        for mr in allele_member_rows:
            allele_mem_w.writerow([
                mr['allele_id'], mr['locus_id'], mr['fid'], mr['label'], mr['kind'], mr['is_allele_rep'], mr['orientation_norm'],
                f"{float(mr['ctx_fwd']):.6f}", f"{float(mr['ctx_rev']):.6f}",
                '' if mr['seq_similarity'] is None else f"{float(mr['seq_similarity']):.6f}",
                f"{float(mr['len_ratio']):.6f}",
            ])
        for a in allele_rows:
            allele_w.writerow([
                a['allele_id'], a['locus_id'], a['allele_rep_fid'], a['display_fid'], a['display_kind'], a['display_label'],
                a['major_allele'], a['member_n'], a['hap_n'], a['asm_hap_n'], a['ref_hap_n'], a['asm_member_n'], a['ref_member_n'],
                a['has_primary_ref_member'], a['has_comparison_ref_member'], a['primary_ref_member_fids'], a['comparison_ref_member_fids'],
                a['rep_seq_len'], a['rep_len'], a['rep_cpg_n'], a['rep_gc_n'], a['rep_pct_gc'], a['rep_oe'], a['allele_rep_strategy'], a.get('allele_cluster_mode', ''),
                ';'.join(a['member_fids']), ';'.join(a['member_labels'])
            ])
            n_alleles += 1
        n_loci_new += 1
        maybe_flush()

    n_threads = int(getattr(args, 'threads', 1))
    if tasks:
        if n_threads > 1:
            backend = resolve_parallel_backend(args, default='process')
            chunksize = max(1, int(getattr(args, 'chunksize', 1)))
            if backend == 'thread':
                with ThreadPoolExecutor(max_workers=n_threads) as ex:
                    for result in ex.map(_process_one_locus_alleles, tasks):
                        write_result(result)
            else:
                ctx = mp.get_context(resolve_mp_start_method(args, default='fork'))
                with ctx.Pool(processes=n_threads, maxtasksperchild=int(getattr(args, 'maxtasksperchild', 0) or 0) or None) as pool:
                    for result in pool.imap(_process_one_locus_alleles, tasks, chunksize=chunksize):
                        write_result(result)
        else:
            for task in tasks:
                write_result(_process_one_locus_alleles(task))

    allele_fh.close()
    allele_mem_fh.close()
    finalize_output_path(out_allele_tmp, args.out_allele)
    finalize_output_path(out_allele_members_tmp, args.out_allele_members)
    summary = {
        'n_alleles': n_alleles,
        'loci_clustered_this_run': n_loci_new,
        'threads': int(getattr(args, 'threads', 1)),
        'parallel_backend': resolve_parallel_backend(args, default='process'),
    }
    summary.update(feature_info)
    print(json_pretty(summary), file=sys.stderr)

def cmd_build_allele_feature_store(args: argparse.Namespace) -> None:
    db_path, shdf_path, meta_path = default_allele_store_paths(args.features, getattr(args, 'out_db', '') or '')
    expected = _allele_store_expected_meta(args.features, anchor_k=args.anchor_k, max_mid_anchors=args.max_mid_anchors)
    reuse_ok = (not bool(getattr(args, 'rebuild', False))) and os.path.exists(db_path) and os.path.exists(shdf_path) and os.path.exists(meta_path) and _allele_store_meta_matches(meta_path, expected)
    if reuse_ok:
        with gzip.open(shdf_path, 'rb') as fh:
            n_features, df = marshal.load(fh)
        print(json_pretty({'reused': True, 'db': db_path, 'shdf': shdf_path, 'meta': meta_path, 'features': int(n_features), 'shingles': len(df)}), file=sys.stderr)
        return
    n_features, df = build_allele_feature_store(
        args.features, db_path, shdf_path, meta_path,
        anchor_k=int(args.anchor_k), max_mid_anchors=int(args.max_mid_anchors),
        batch_size=int(getattr(args, 'batch_size', 4096) or 4096),
    )
    print(json_pretty({'reused': False, 'db': db_path, 'shdf': shdf_path, 'meta': meta_path, 'features': int(n_features), 'shingles': len(df)}), file=sys.stderr)


def _write_locus_list_rows(path: str, rows: Sequence[Tuple[int, str, str, int]]) -> None:
    fh, w = write_tsv_header(path, ['global_index', 'locus_id', 'medoid_fid', 'member_n'])
    try:
        for row in rows:
            w.writerow(list(row))
    finally:
        fh.close()


def cmd_export_allele_loci(args: argparse.Namespace) -> None:
    medoid_by_locus = read_locus_catalog_medoid(args.locus_catalog)
    counts = count_members_by_locus_stream(args.locus_members)
    locus_ids = [lid for lid in sorted(medoid_by_locus, key=locus_sort_key) if counts.get(lid, 0) > 0]
    rows = [(i, lid, medoid_by_locus[lid], int(counts.get(lid, 0))) for i, lid in enumerate(locus_ids)]
    if getattr(args, 'out', ''):
        _write_locus_list_rows(str(args.out), rows)
    shard_count = int(getattr(args, 'shard_count', 0) or 0)
    shard_paths: List[str] = []
    if shard_count > 0:
        out_dir = str(getattr(args, 'shard_dir', '') or '')
        if not out_dir:
            raise ValueError('--shard-dir is required when --shard-count > 0')
        os.makedirs(out_dir, exist_ok=True)
        mode = str(getattr(args, 'shard_mode', 'contiguous') or 'contiguous')
        prefix = str(getattr(args, 'shard_prefix', 'allele_loci') or 'allele_loci')
        n = len(rows)
        for idx in range(shard_count):
            if mode == 'modulo':
                part = [row for j, row in enumerate(rows) if (j % shard_count) == idx]
            elif mode == 'contiguous':
                lo = (n * idx) // shard_count
                hi = (n * (idx + 1)) // shard_count
                part = rows[lo:hi]
            else:
                raise ValueError(f'Unknown --shard-mode: {mode}')
            path = os.path.join(out_dir, f'{prefix}.shard_{idx:03d}_of_{shard_count:03d}.tsv')
            _write_locus_list_rows(path, part)
            shard_paths.append(path)
    print(json_pretty({'loci': len(rows), 'out': str(getattr(args, 'out', '') or ''), 'shard_count': shard_count, 'shard_paths': shard_paths}), file=sys.stderr)


def iter_tsv_dict_rows(path: str) -> Iterator[Dict[str, str]]:
    with open_text(path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            yield {str(k): ('' if v is None else str(v)) for k, v in row.items()}


def tsv_header(path: str) -> List[str]:
    with open_text(path, 'rt') as fh:
        raw = fh.readline()
    return raw.rstrip('\n').split('\t') if raw else []


class TsvWriterPool:
    def __init__(self, header_by_path: Dict[str, List[str]], max_open: int = 64):
        self.header_by_path = header_by_path
        self.max_open = max(1, int(max_open))
        self.handles: "OrderedDict[str, Tuple[object, object]]" = OrderedDict()
        self.initialized: set = set()

    def writer(self, path: str):
        if path in self.handles:
            fh, w = self.handles.pop(path)
            self.handles[path] = (fh, w)
            return w
        if len(self.handles) >= self.max_open:
            old_path, (old_fh, _old_w) = self.handles.popitem(last=False)
            old_fh.close()
        if path in self.initialized:
            ensure_parent(path)
            fh = open_text(path, 'at')
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        else:
            fh, w = write_tsv_header(path, self.header_by_path[path])
            self.initialized.add(path)
        self.handles[path] = (fh, w)
        return w

    def close(self) -> None:
        while self.handles:
            _path, (fh, _w) = self.handles.popitem(last=False)
            fh.close()


class JsonlGzWriterPool:
    def __init__(self, paths: Sequence[str], max_open: int = 64):
        self.path_set = set(str(p) for p in paths)
        self.max_open = max(1, int(max_open))
        self.handles: "OrderedDict[str, object]" = OrderedDict()

    def write(self, path: str, line: str) -> None:
        if path not in self.path_set:
            raise KeyError(path)
        if path in self.handles:
            fh = self.handles.pop(path)
            self.handles[path] = fh
        else:
            if len(self.handles) >= self.max_open:
                _old_path, old_fh = self.handles.popitem(last=False)
                old_fh.close()
            ensure_parent(path)
            fh = gzip.open(path, 'at') if path.endswith('.gz') else open(path, 'at')
            self.handles[path] = fh
        self.handles[path].write(line)

    def close(self) -> None:
        while self.handles:
            _path, fh = self.handles.popitem(last=False)
            fh.close()


def cmd_export_locus_chunks(args: argparse.Namespace) -> None:
    os.makedirs(args.out_dir, exist_ok=True)
    loc_header = tsv_header(args.locus_catalog)
    mem_header = tsv_header(args.locus_members)
    if 'locus_id' not in loc_header or 'locus_id' not in mem_header or 'fid' not in mem_header:
        raise ValueError('locus catalog/members must contain locus_id and members must contain fid')
    counts = count_members_by_locus_stream(args.locus_members)
    loc_rows: Dict[str, Dict[str, str]] = {}
    for row in iter_tsv_dict_rows(args.locus_catalog):
        lid = str(row.get('locus_id') or '')
        if lid and counts.get(lid, 0) > 0:
            loc_rows[lid] = row
    locus_ids = sorted(loc_rows, key=locus_sort_key)
    chunk_count = int(getattr(args, 'chunk_count', 0) or 0)
    if chunk_count <= 0:
        chunk_size = max(1, int(getattr(args, 'chunk_size', 0) or 0))
        chunk_count = max(1, int(math.ceil(len(locus_ids) / float(chunk_size))))
    mode = str(getattr(args, 'chunk_mode', 'contiguous') or 'contiguous')
    locus_rank = {lid: i for i, lid in enumerate(locus_ids)}
    chunks: List[List[str]] = []
    for idx in range(chunk_count):
        if mode == 'modulo':
            part = [lid for j, lid in enumerate(locus_ids) if (j % chunk_count) == idx]
        elif mode == 'contiguous':
            lo = (len(locus_ids) * idx) // chunk_count
            hi = (len(locus_ids) * (idx + 1)) // chunk_count
            part = locus_ids[lo:hi]
        else:
            raise ValueError(f'Unknown --chunk-mode: {mode}')
        chunks.append(part)
    locus_to_chunk: Dict[str, int] = {}
    for idx, lids in enumerate(chunks):
        for lid in lids:
            locus_to_chunk[lid] = idx

    prefix = str(getattr(args, 'prefix', 'chunk') or 'chunk')
    chunk_meta: List[Dict[str, object]] = []
    loc_paths: Dict[int, str] = {}
    mem_paths: Dict[int, str] = {}
    ids_paths: Dict[int, str] = {}
    feat_paths: Dict[int, str] = {}
    for idx, lids in enumerate(chunks):
        tag = f'{prefix}_{idx:04d}_of_{chunk_count:04d}'
        ids_path = os.path.join(args.out_dir, tag + '.locus_ids.tsv')
        loc_path = os.path.join(args.out_dir, tag + '.locus_catalog.tsv.gz')
        mem_path = os.path.join(args.out_dir, tag + '.locus_members.tsv.gz')
        feat_path = os.path.join(args.out_dir, tag + '.features.jsonl.gz')
        ids_paths[idx] = ids_path
        loc_paths[idx] = loc_path
        mem_paths[idx] = mem_path
        feat_paths[idx] = feat_path
        with open_text(ids_path, 'wt') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['global_index', 'chunk_index', 'locus_id', 'member_n'])
            for lid in lids:
                w.writerow([locus_rank.get(lid, ''), idx, lid, int(counts.get(lid, 0))])
        loc_fh, loc_w = write_tsv_header(loc_path, loc_header)
        try:
            for lid in lids:
                row = loc_rows[lid]
                loc_w.writerow([row.get(c, '') for c in loc_header])
        finally:
            loc_fh.close()

    mem_pool = TsvWriterPool({p: mem_header for p in mem_paths.values()}, max_open=int(getattr(args, 'max_open_files', 64) or 64))
    fid_to_chunk: Dict[str, int] = {}
    member_counts = Counter()
    try:
        for row in iter_tsv_dict_rows(args.locus_members):
            lid = str(row.get('locus_id') or '')
            if lid not in locus_to_chunk:
                continue
            idx = locus_to_chunk[lid]
            mem_pool.writer(mem_paths[idx]).writerow([row.get(c, '') for c in mem_header])
            fid = str(row.get('fid') or '')
            if fid:
                fid_to_chunk[fid] = idx
            member_counts[idx] += 1
    finally:
        mem_pool.close()

    feature_counts = Counter()
    if bool(getattr(args, 'write_feature_subsets', False)):
        if not str(getattr(args, 'features', '') or ''):
            raise ValueError('--features is required with --write-feature-subsets')

        for pth in feat_paths.values():
            if os.path.exists(pth):
                os.remove(pth)
        feat_pool = JsonlGzWriterPool(list(feat_paths.values()), max_open=int(getattr(args, 'max_open_files', 64) or 64))
        try:
            with open_text(args.features, 'rt') as fh:
                for raw in fh:
                    if not raw.strip():
                        continue
                    try:
                        rec = json.loads(raw)
                    except Exception:
                        continue
                    fid = str(rec.get('fid') or '')
                    idx = fid_to_chunk.get(fid)
                    if idx is None:
                        continue
                    feat_pool.write(feat_paths[idx], raw)
                    feature_counts[idx] += 1
        finally:
            feat_pool.close()

    catalog_path = os.path.join(args.out_dir, prefix + '.chunk_catalog.tsv')
    mf_fh, mf_w = write_tsv_header(catalog_path, ['chunk_index', 'chunk_count', 'locus_n', 'member_n', 'feature_n', 'locus_ids', 'locus_catalog', 'locus_members', 'features'])
    try:
        for idx, lids in enumerate(chunks):
            mf_w.writerow([
                idx, chunk_count, len(lids), int(member_counts.get(idx, 0)), int(feature_counts.get(idx, 0)) if bool(getattr(args, 'write_feature_subsets', False)) else '',
                ids_paths[idx], loc_paths[idx], mem_paths[idx], feat_paths[idx] if bool(getattr(args, 'write_feature_subsets', False)) else '',
            ])
            chunk_meta.append({'chunk_index': idx, 'locus_n': len(lids), 'member_n': int(member_counts.get(idx, 0)), 'feature_n': int(feature_counts.get(idx, 0)) if bool(getattr(args, 'write_feature_subsets', False)) else None})
    finally:
        mf_fh.close()
    print(json_pretty({'loci': len(locus_ids), 'chunk_count': chunk_count, 'catalog': os.path.abspath(catalog_path), 'write_feature_subsets': bool(getattr(args, 'write_feature_subsets', False)), 'chunks_preview': chunk_meta[:5]}), file=sys.stderr)


def expand_input_paths(paths: Sequence[str]) -> List[str]:
    out: List[str] = []
    for path in paths:
        matches = sorted(glob.glob(path))
        if matches:
            out.extend(matches)
        else:
            out.append(path)
    return out


def cmd_merge_allele_shards(args: argparse.Namespace) -> None:
    locus_order = read_locus_ids_file(args.locus_list)
    if not locus_order:
        raise ValueError(f'No locus IDs read from --locus-list {args.locus_list}')
    allele_shards = expand_input_paths(args.allele_shards)
    member_shards = expand_input_paths(args.allele_member_shards)
    missing = [p for p in allele_shards + member_shards if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f'Missing shard file(s): {missing[:5]}')
    require_all = not bool(getattr(args, 'allow_missing_loci', False))
    allele_info = merge_sharded_tsv_by_locus(allele_shards, locus_order, args.out_allele, ALLELE_CATALOG_HEADER, require_all_loci=require_all)
    member_info = merge_sharded_tsv_by_locus(member_shards, locus_order, args.out_allele_members, ALLELE_MEMBER_HEADER, require_all_loci=require_all)
    print(json_pretty({'locus_list_n': len(locus_order), 'allele_shards': len(allele_shards), 'allele_member_shards': len(member_shards), 'allele_catalog': allele_info, 'allele_members': member_info}), file=sys.stderr)


def group_members(df: pd.DataFrame, key: str) -> Dict[str, List[Dict[str, object]]]:
    out: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for _, row in df.iterrows():
        rec = {c: row[c] for c in df.columns}
        out[str(row[key])].append(rec)
    return out


def locus_annotation_summary(member_fids: Sequence[str], feat_idx: Dict[str, Dict[str, object]]) -> Dict[str, object]:
    anns = []
    for fid in member_fids:
        f = feat_idx[fid]
        if f.get('primary_chr'):
            anns.append((str(f['primary_chr']), int(f['primary_start1']), int(f['primary_end1'])))
    if not anns:
        return {'source': 'none'}
    chrom = dominant_value([a[0] for a in anns])
    anns2 = [a for a in anns if a[0] == chrom]
    starts = [a[1] for a in anns2]
    ends = [a[2] for a in anns2]
    mids = [int(round((a[1] + a[2]) / 2.0)) for a in anns2]
    return {
        'source': 'member_median',
        'chrom': chrom,
        'start_median1': median_int(starts),
        'mid_median1': median_int(mids),
        'end_median1': median_int(ends),
        'n_annotated_members': len(anns2),
    }


def emit_merged_fasta(path: str, allele_df: pd.DataFrame, feat_idx: Dict[str, Dict[str, object]]) -> None:
    ensure_parent(path)
    with open_text(path, 'wt') as out:
        for _, row in allele_df.iterrows():
            fid = str(row['allele_rep_fid'])
            seq = feat_idx[fid].get('source_seq')
            if not seq:
                continue
            write_fasta_record(out, fid, str(seq))


def cmd_strict_genotype_prod(args: argparse.Namespace) -> None:
    from pancgi_genotyping.pipeline import run
    run(args)


def validate_strict_matrix_compat(df: pd.DataFrame, *, label_order: List[str], expected_ids: Sequence[str], name: str) -> pd.DataFrame:
    if list(df.columns) != label_order or df.index.has_duplicates or list(df.index) != [str(x) for x in expected_ids]:
        raise ValueError(f'{name}: identifiers or column order differ from the catalog')
    if not set(df.to_numpy().ravel()) <= {'0', '1', 'NA'}:
        raise ValueError(f'{name}: invalid or missing genotype cell')
    missing_cols = [lab for lab in label_order if lab not in df.columns]
    if missing_cols:
        raise RuntimeError(f"{name} is missing label columns required by catalog: {missing_cols[:10]}{'...' if len(missing_cols) > 10 else ''}")
    out = df.loc[:, label_order].copy()
    expected = [str(x) for x in expected_ids]
    missing_ids = [x for x in expected if x not in out.index]
    if missing_ids:
        raise RuntimeError(f"{name} is missing ids required by catalogs: {missing_ids[:10]}{'...' if len(missing_ids) > 10 else ''}")
    return out


def cmd_emit_merged(args: argparse.Namespace) -> None:
    feat_idx: Dict[str, Dict[str, object]] = {}
    feature_label_order: List[str] = []
    seen_labels = set()
    for rec in iter_feature_jsonl(args.features):
        feat_idx[str(rec['fid'])] = rec
        lab = str(rec['label'])
        if lab not in seen_labels:
            feature_label_order.append(lab)
            seen_labels.add(lab)

    catalog_labels = iter_catalog_labels(args.catalog)
    strict_locus_gt = pd.read_csv(args.locus_genotype_tsv, sep='	', dtype=str, index_col=0, keep_default_na=False)
    strict_allele_gt = pd.read_csv(args.allele_genotype_tsv, sep='	', dtype=str, index_col=0, keep_default_na=False)
    strict_locus_gt.index = strict_locus_gt.index.map(str)
    strict_allele_gt.index = strict_allele_gt.index.map(str)
    label_kind = build_label_kind_map(args.catalog, feat_idx)
    label_order = list(catalog_labels)
    assembly_labels = [lab for lab in label_order if is_assembly_kind(label_kind.get(lab, 'assembly'))]
    total_assemblies = len(assembly_labels)

    loc = read_tsv(args.locus_catalog)
    loc_mem = read_tsv(args.locus_members)
    al = read_tsv(args.allele_catalog)
    al_mem = read_tsv(args.allele_members)
    for df in (loc, loc_mem, al, al_mem):
        for col in ('locus_id', 'allele_id', 'fid', 'allele_rep_fid', 'medoid_fid', 'lead_fid', 'display_fid', 'anchor_ref_fid'):
            if col in df.columns:
                df[col] = df[col].astype(str)

    loc_by_id = {str(r['locus_id']): {c: r[c] for c in loc.columns} for _, r in loc.iterrows()}
    locus_members = group_members(loc_mem, 'locus_id')
    allele_members = group_members(al_mem, 'allele_id')
    strict_locus_gt = validate_strict_matrix_compat(strict_locus_gt, label_order=label_order, expected_ids=loc['locus_id'].astype(str).tolist(), name='strict locus genotype matrix')
    strict_allele_gt = validate_strict_matrix_compat(strict_allele_gt, label_order=label_order, expected_ids=al['allele_id'].astype(str).tolist(), name='strict allele genotype matrix')

    locus_primary: Dict[str, Dict[str, object]] = {}
    for locus_id, rows in locus_members.items():
        fids = [str(r['fid']) for r in rows]
        locus_primary[locus_id] = locus_annotation_summary(fids, feat_idx)

    base_header = [
        'merged_id', 'allele_id', 'locus_id', 'locus_type', 'anchor_ref_fid', 'locus_insert_anchor_member_n',
        'rep_anchor_assign_rule', 'rep_anchor_confidence', 'rep_insert_anchor_flag', 'rep_insert_anchor_proj_type',
        'major_allele', 'allele_rep_fid', 'allele_rep_strategy', 'allele_cluster_mode', 'display_fid',
        'lead_fid', 'medoid_fid', 'locus_display_fid',
        'has_primary_ref', 'has_comparison_ref', 'primary_ref_fids', 'comparison_ref_fids',
        'allele_member_n_total', 'allele_member_n_asm', 'allele_member_n_ref',
        'allele_hap_n_total', 'allele_hap_n_asm', 'allele_hap_n_ref',
        'locus_member_n_total', 'locus_member_n_asm', 'locus_member_n_ref',
        'locus_hap_n_total', 'locus_hap_n_asm', 'locus_hap_n_ref',
        'n_total_assemblies', 'allele_asm_freq', 'allele_freq_class', 'ref_only_flag',
        'rep_fid', 'rep_label', 'rep_kind', 'rep_is_assembly', 'rep_seq_len_bp', 'rep_cpgi_len_bp', 'rep_len', 'rep_cpg_n', 'rep_gc_n', 'rep_pct_gc', 'rep_oe',
        'rep_asm_mid0', 'rep_graph_nodeints',
        'locus_primary_source', 'locus_primary_candidate_fid', 'locus_primary_chr', 'locus_primary_start_median1', 'locus_primary_mid_median1', 'locus_primary_end_median1',
        'locus_nonref_site_source', 'locus_nonref_site_chrom', 'locus_nonref_site_start1', 'locus_nonref_site_end1',
        'rep_primary_source', 'rep_primary_candidate_fid', 'rep_primary_type', 'rep_primary_chr', 'rep_primary_start1', 'rep_primary_mid1', 'rep_primary_end1', 'rep_primary_anchor1',
        'rep_sv_overlap_class', 'rep_sv_ins_n', 'rep_sv_ids', 'rep_sv_primary_sites', 'rep_sv_primary_intervals', 'rep_sv_asm_intervals', 'rep_sv_source', 'rep_sv_longest_id', 'rep_sv_longest_chrom', 'rep_sv_longest_site1', 'rep_sv_longest_start1', 'rep_sv_longest_end1', 'rep_sv_longest_len',
        'allele_member_ids', 'allele_member_labels', 'allele_member_graph_positions', 'allele_member_asm_midpoints', 'allele_member_primary_midpoints',
        'locus_member_ids', 'locus_member_labels', 'locus_member_graph_positions', 'locus_member_asm_midpoints', 'locus_member_primary_midpoints',
        'strict_genotype_mode', 'allele_callable_asm_n', 'allele_callable_total_n', 'locus_callable_asm_n', 'locus_callable_total_n', 'allele_gt_1_asm_n', 'allele_gt_0_asm_n', 'allele_gt_na_asm_n', 'allele_gt_1_total_n', 'allele_gt_0_total_n', 'allele_gt_na_total_n', 'locus_gt_1_asm_n', 'locus_gt_0_asm_n', 'locus_gt_na_asm_n', 'locus_gt_1_total_n', 'locus_gt_0_total_n', 'locus_gt_na_total_n'
    ]
    extra_header = [
        'locus_primary_cpgi_id', 'merged_cpgi_n', 'merged_cpgi_ids', 'locus_cpgi_n', 'locus_cpgi_ids',
        'contains_primary_member', 'contains_comparison_member', 'allele_freq',
        'rep_primary_primary_source', 'position_annotation_source', 'rep_primary_primary_chr', 'rep_primary_primary_start1', 'rep_primary_primary_end1', 'rep_primary_primary_site1',
        'rep_hal_multimap_label', 'rep_hal_admit_rule', 'rep_hal_multimap_n', 'rep_hal_primary_chr', 'rep_hal_primary_start0', 'rep_hal_primary_end0', 'rep_hal_primary_start1', 'rep_hal_primary_end1', 'rep_hal_query_cov', 'rep_hal_identity', 'rep_hal_all_intervals0', 'rep_hal_all_intervals',
        'rep_sv_all_n', 'rep_sv_all_ids', 'rep_sv_all_types', 'rep_sv_all_classes', 'rep_sv_all_primary_sites', 'rep_sv_all_primary_intervals', 'rep_sv_all_asm_intervals',
        'rep_sv_all_overlap_bp', 'rep_sv_all_overlap_pct_cpgi', 'rep_sv_all_overlap_pct_sv',
        'rep_sv_all_TR', 'rep_sv_all_CONFORMATION', 'rep_sv_all_SD', 'rep_sv_all_ITYPE_N', 'rep_sv_all_DTYPE_N', 'rep_sv_all_FAM_N', 'rep_sv_all_source'
    ]
    out_header = base_header + extra_header + label_order
    out_fh, out_w = write_tsv_header(args.out_tsv, out_header)

    big_rows = []
    for _, arow in al.iterrows():
        locus_id = str(arow['locus_id'])
        allele_id = str(arow['allele_id'])
        rep_fid = str(arow['allele_rep_fid'])
        rep_feat = feat_idx[rep_fid]
        locus_row = loc_by_id[locus_id]
        locus_member_rows = locus_members.get(locus_id, [])
        allele_member_rows = allele_members.get(allele_id, [])
        rep_member_row = next((r for r in locus_member_rows if str(r.get('fid')) == rep_fid), {})

        locus_member_ids = [str(r['fid']) for r in locus_member_rows]
        allele_member_ids = [str(r['fid']) for r in allele_member_rows]
        locus_labels_all = sorted({str(r['label']) for r in locus_member_rows})
        allele_labels_all = sorted({str(r['label']) for r in allele_member_rows})
        locus_labels_asm = sorted({str(r['label']) for r in locus_member_rows if is_assembly_kind(r['kind'])})
        allele_labels_asm = sorted({str(r['label']) for r in allele_member_rows if is_assembly_kind(r['kind'])})
        locus_labels_ref = sorted({str(r['label']) for r in locus_member_rows if is_reference_kind(r['kind'])})
        allele_labels_ref = sorted({str(r['label']) for r in allele_member_rows if is_reference_kind(r['kind'])})
        primary_ref_fids = [fid for fid in locus_member_ids if str(feat_idx[fid]['kind']) == 'reference_primary']
        comparison_ref_fids = [fid for fid in locus_member_ids if str(feat_idx[fid]['kind']) == 'reference_comparison']

        rep_primary_source = 'unresolved'
        rep_primary_candidate_fid = rep_fid if rep_feat.get('primary_chr') else ''
        rep_primary_type = 'unresolved'
        rep_primary_chr = rep_feat.get('primary_chr')
        rep_primary_start1 = rep_feat.get('primary_start1')
        rep_primary_end1 = rep_feat.get('primary_end1')
        rep_primary_anchor1 = ''

        if rep_feat.get('hal_primary_chr'):
            rep_primary_source = 'hal_liftover'
            rep_primary_chr = rep_feat.get('hal_primary_chr')
            rep_primary_start1 = rep_feat.get('hal_primary_start1')
            rep_primary_end1 = rep_feat.get('hal_primary_end1')
            rep_primary_type = 'span'
        elif rep_primary_chr:
            rep_primary_source = 'reference_coordinate'
            rep_primary_type = 'span'
        if not rep_primary_chr:
            summary = locus_primary.get(locus_id, {'source': 'none'})
            if summary.get('chrom'):
                rep_primary_source = str(summary.get('source', 'formal_locus_members'))
                rep_primary_candidate_fid = ''
                rep_primary_type = 'span'
                rep_primary_chr = summary.get('chrom')
                rep_primary_start1 = summary.get('start_median1')
                rep_primary_end1 = summary.get('end_median1')
                rep_primary_anchor1 = summary.get('mid_median1')
            elif locus_row.get('nonref_site_chrom'):
                rep_primary_source = str(locus_row.get('nonref_site_source') or 'provided_sv_callset')
                rep_primary_type = 'site'
                rep_primary_chr = locus_row.get('nonref_site_chrom')
                rep_primary_start1 = locus_row.get('nonref_site_start1')
                rep_primary_end1 = locus_row.get('nonref_site_end1')

        rep_primary_mid1 = ''
        if rep_primary_start1 not in (None, '', '.') and rep_primary_end1 not in (None, '', '.'):
            rep_primary_mid1 = int(round((int(rep_primary_start1) + int(rep_primary_end1)) / 2.0))

        rep_sv_overlap_class = rep_feat.get('sv_ins_overlap_class', '')
        rep_sv_ins_n = rep_feat.get('sv_ins_n', 0)
        rep_sv_ids = rep_feat.get('sv_ins_ids', '')
        rep_sv_primary_sites = rep_feat.get('sv_ins_primary_sites', '')
        rep_sv_primary_intervals = rep_feat.get('sv_ins_primary_intervals', '')
        rep_sv_asm_intervals = rep_feat.get('sv_ins_asm_intervals', '')
        rep_sv_source = rep_feat.get('sv_ins_source', '')
        rep_sv_longest_id = rep_feat.get('sv_ins_longest_id', '')
        rep_sv_longest_chrom = rep_feat.get('sv_ins_longest_chrom', '')
        rep_sv_longest_site1 = rep_feat.get('sv_ins_longest_site1', '')
        rep_sv_longest_start1 = rep_feat.get('sv_ins_longest_start1', '')
        rep_sv_longest_end1 = rep_feat.get('sv_ins_longest_end1', '')
        rep_sv_longest_len = rep_feat.get('sv_ins_longest_len', '')
        rep_sv_all_n = rep_feat.get('sv_overlap_n', 0)
        rep_sv_all_ids = rep_feat.get('sv_overlap_ids', '')
        rep_sv_all_types = rep_feat.get('sv_overlap_types', '')
        rep_sv_all_classes = rep_feat.get('sv_overlap_classes', '')
        rep_sv_all_primary_sites = rep_feat.get('sv_overlap_primary_sites', '')
        rep_sv_all_primary_intervals = rep_feat.get('sv_overlap_primary_intervals', '')
        rep_sv_all_asm_intervals = rep_feat.get('sv_overlap_asm_intervals', '')
        rep_sv_all_overlap_bp = rep_feat.get('sv_overlap_bp', '')
        rep_sv_all_overlap_pct_cpgi = rep_feat.get('sv_overlap_pct_cpgi', '')
        rep_sv_all_overlap_pct_sv = rep_feat.get('sv_overlap_pct_sv', '')
        rep_sv_all_TR = rep_feat.get('sv_overlap_TR', '')
        rep_sv_all_CONFORMATION = rep_feat.get('sv_overlap_CONFORMATION', '')
        rep_sv_all_SD = rep_feat.get('sv_overlap_SD', '')
        rep_sv_all_ITYPE_N = rep_feat.get('sv_overlap_ITYPE_N', '')
        rep_sv_all_DTYPE_N = rep_feat.get('sv_overlap_DTYPE_N', '')
        rep_sv_all_FAM_N = rep_feat.get('sv_overlap_FAM_N', '')
        rep_sv_all_source = rep_feat.get('sv_overlap_source', '')
        allele_asm_hap_n = int(arow['asm_hap_n']) if 'asm_hap_n' in arow else len(allele_labels_asm)
        allele_ref_hap_n = int(arow['ref_hap_n']) if 'ref_hap_n' in arow else len(allele_labels_ref)
        allele_hap_n_total = int(arow['hap_n']) if 'hap_n' in arow else len(allele_labels_all)
        allele_member_n_total = int(arow['member_n']) if 'member_n' in arow else len(allele_member_ids)
        allele_member_n_asm = int(arow['asm_member_n']) if 'asm_member_n' in arow else sum(1 for fid in allele_member_ids if is_assembly_kind(feat_idx[fid]['kind']))
        allele_member_n_ref = int(arow['ref_member_n']) if 'ref_member_n' in arow else allele_member_n_total - allele_member_n_asm
        locus_hap_n_asm = int(locus_row['asm_hap_n']) if 'asm_hap_n' in locus_row else len(locus_labels_asm)
        locus_hap_n_ref = int(locus_row['ref_hap_n']) if 'ref_hap_n' in locus_row else len(locus_labels_ref)
        locus_hap_n_total = int(locus_row['hap_n']) if 'hap_n' in locus_row else len(locus_labels_all)
        locus_member_n_total = int(locus_row['member_n']) if 'member_n' in locus_row else len(locus_member_ids)
        locus_member_n_asm = sum(1 for fid in locus_member_ids if is_assembly_kind(feat_idx[fid]['kind']))
        locus_member_n_ref = locus_member_n_total - locus_member_n_asm

        try:
            major_flag = int(arow['major_allele']) if 'major_allele' in arow and str(arow['major_allele']) not in ('', 'nan', 'None') else 0
        except Exception:
            major_flag = 0
        allele_freq = (float(allele_asm_hap_n) / float(total_assemblies)) if total_assemblies else None
        has_ref = 1 if (len(primary_ref_fids) + len(comparison_ref_fids)) > 0 else 0
        allele_gt_row = strict_allele_gt.loc[allele_id]
        locus_gt_row = strict_locus_gt.loc[locus_id]
        allele_vals_asm = [str(allele_gt_row[lab]) for lab in assembly_labels]
        allele_vals_total = [str(allele_gt_row[lab]) for lab in label_order]
        allele_callable_asm_n = sum(1 for v in allele_vals_asm if v != 'NA')
        allele_callable_total_n = sum(1 for v in allele_vals_total if v != 'NA')
        allele_gt_1_asm_n = sum(1 for v in allele_vals_asm if v == '1')
        allele_gt_0_asm_n = sum(1 for v in allele_vals_asm if v == '0')
        allele_gt_na_asm_n = sum(1 for v in allele_vals_asm if v == 'NA')
        allele_gt_1_total_n = sum(1 for v in allele_vals_total if v == '1')
        allele_gt_0_total_n = sum(1 for v in allele_vals_total if v == '0')
        allele_gt_na_total_n = sum(1 for v in allele_vals_total if v == 'NA')
        locus_vals_asm = [str(locus_gt_row[lab]) for lab in assembly_labels]
        locus_vals_total = [str(locus_gt_row[lab]) for lab in label_order]
        locus_callable_asm_n = sum(1 for v in locus_vals_asm if v != 'NA')
        locus_callable_total_n = sum(1 for v in locus_vals_total if v != 'NA')
        locus_gt_1_asm_n = sum(1 for v in locus_vals_asm if v == '1')
        locus_gt_0_asm_n = sum(1 for v in locus_vals_asm if v == '0')
        locus_gt_na_asm_n = sum(1 for v in locus_vals_asm if v == 'NA')
        locus_gt_1_total_n = sum(1 for v in locus_vals_total if v == '1')
        locus_gt_0_total_n = sum(1 for v in locus_vals_total if v == '0')
        locus_gt_na_total_n = sum(1 for v in locus_vals_total if v == 'NA')

        rep_primary_primary_source = 'unresolved'
        rep_primary_primary_chr = ''
        rep_primary_primary_start1 = ''
        rep_primary_primary_end1 = ''
        rep_primary_primary_site1 = ''
        rep_sv_class_tokens = {x for x in str(rep_sv_overlap_class or '').split(';') if x}
        if 'inside_ins' in rep_sv_class_tokens and rep_sv_longest_chrom and (rep_sv_longest_site1 not in ('', None, '.', 'nan') or rep_sv_longest_start1 not in ('', None, '.', 'nan') or rep_sv_longest_end1 not in ('', None, '.', 'nan')):
            rep_primary_primary_source = 'provided_sv_callset'
            rep_primary_primary_chr = rep_sv_longest_chrom
            rep_primary_primary_site1 = rep_sv_longest_site1 if rep_sv_longest_site1 not in ('', None, '.', 'nan') else ''
            rep_primary_primary_start1 = rep_sv_longest_start1 if rep_sv_longest_start1 not in ('', None, '.', 'nan') else rep_primary_primary_site1
            rep_primary_primary_end1 = rep_sv_longest_end1 if rep_sv_longest_end1 not in ('', None, '.', 'nan') else rep_primary_primary_site1
        elif rep_feat.get('hal_primary_chr'):
            rep_primary_primary_source = 'hal_liftover'
            rep_primary_primary_chr = rep_feat.get('hal_primary_chr')
            rep_primary_primary_start1 = rep_feat.get('hal_primary_start1')
            rep_primary_primary_end1 = rep_feat.get('hal_primary_end1')
            if rep_primary_primary_start1 not in ('', None, '.', 'nan') and rep_primary_primary_end1 not in ('', None, '.', 'nan'):
                rep_primary_primary_site1 = int(round((int(rep_primary_primary_start1) + int(rep_primary_primary_end1)) / 2.0))
        elif rep_primary_chr:
            rep_primary_primary_source = rep_primary_source or 'formal_coordinate'
            rep_primary_primary_chr = rep_primary_chr
            rep_primary_primary_start1 = rep_primary_start1
            rep_primary_primary_end1 = rep_primary_end1
            if rep_primary_anchor1 not in ('', None, '.', 'nan'):
                rep_primary_primary_site1 = rep_primary_anchor1
            elif rep_primary_start1 not in ('', None, '.', 'nan') and rep_primary_end1 not in ('', None, '.', 'nan'):
                rep_primary_primary_site1 = int(round((int(rep_primary_start1) + int(rep_primary_end1)) / 2.0))
        elif str(locus_row.get('anchor_type', '')) == 'primary_ref':
            anchor_ref_fid = str(locus_row.get('anchor_ref_fid') or '').strip()
            if anchor_ref_fid:
                try:
                    anchor_feat = feat_idx.get(anchor_ref_fid, {})
                    anchor_chr = str(anchor_feat.get('primary_chr') or '')
                    anchor_start1 = anchor_feat.get('primary_start1', '')
                    anchor_end1 = anchor_feat.get('primary_end1', '')
                    if not anchor_chr or anchor_start1 in ('', None, '.', 'nan') or anchor_end1 in ('', None, '.', 'nan'):
                        core = base.parse_fid(anchor_ref_fid)
                        anchor_chr = str(core.get('contig') or '')
                        anchor_start1 = int(core.get('start0') or 0) + 1
                        anchor_end1 = int(core.get('end0') or 0)
                    if anchor_chr:
                        rep_primary_primary_source = 'locus_anchor_ref'
                        rep_primary_primary_chr = anchor_chr
                        rep_primary_primary_start1 = anchor_start1
                        rep_primary_primary_end1 = anchor_end1
                        if anchor_start1 not in ('', None, '.', 'nan') and anchor_end1 not in ('', None, '.', 'nan'):
                            rep_primary_primary_site1 = int(round((int(anchor_start1) + int(anchor_end1)) / 2.0))
                except Exception:
                    pass

        row = {
            'merged_id': rep_fid,
            'allele_id': allele_id,
            'locus_id': locus_id,
            'locus_type': locus_row.get('anchor_type', 'unknown'),
            'anchor_ref_fid': '' if str(locus_row.get('anchor_ref_fid', '')) in ('nan', 'None') else str(locus_row.get('anchor_ref_fid', '')),
            'locus_insert_anchor_member_n': locus_row.get('insert_anchor_member_n', ''),
            'rep_anchor_assign_rule': rep_member_row.get('anchor_assign_rule', ''),
            'rep_anchor_confidence': rep_member_row.get('anchor_confidence', ''),
            'rep_insert_anchor_flag': rep_member_row.get('insert_anchor_flag', ''),
            'rep_insert_anchor_proj_type': rep_member_row.get('insert_anchor_proj_type', ''),
            'major_allele': major_flag,
            'allele_rep_fid': rep_fid,
            'allele_rep_strategy': arow.get('allele_rep_strategy', ''),
            'allele_cluster_mode': arow.get('allele_cluster_mode', ''),
            'display_fid': rep_fid,
            'lead_fid': str(locus_row.get('lead_fid', '')),
            'medoid_fid': str(locus_row.get('medoid_fid', '')),
            'locus_display_fid': str(locus_row.get('display_fid', '')),
            'has_primary_ref': 1 if primary_ref_fids else 0,
            'has_comparison_ref': 1 if comparison_ref_fids else 0,
            'primary_ref_fids': ';'.join(primary_ref_fids),
            'comparison_ref_fids': ';'.join(comparison_ref_fids),
            'allele_member_n_total': allele_member_n_total,
            'allele_member_n_asm': allele_member_n_asm,
            'allele_member_n_ref': allele_member_n_ref,
            'allele_hap_n_total': allele_hap_n_total,
            'allele_hap_n_asm': allele_asm_hap_n,
            'allele_hap_n_ref': allele_ref_hap_n,
            'locus_member_n_total': locus_member_n_total,
            'locus_member_n_asm': locus_member_n_asm,
            'locus_member_n_ref': locus_member_n_ref,
            'locus_hap_n_total': locus_hap_n_total,
            'locus_hap_n_asm': locus_hap_n_asm,
            'locus_hap_n_ref': locus_hap_n_ref,
            'n_total_assemblies': total_assemblies,
            'allele_asm_freq': '' if allele_freq is None else f'{allele_freq:.6f}',
            'allele_freq_class': allele_freq_class(allele_asm_hap_n, total_assemblies, major_flag, has_ref=has_ref),
            'ref_only_flag': 1 if allele_asm_hap_n == 0 and has_ref else 0,
            'rep_fid': rep_fid,
            'rep_label': rep_feat.get('label'),
            'rep_kind': rep_feat.get('kind'),
            'rep_is_assembly': 1 if is_assembly_kind(rep_feat.get('kind')) else 0,
            'rep_seq_len_bp': int(rep_feat.get('source_seq_len') or rep_feat.get('asm_len') or 0),
            'rep_cpgi_len_bp': int(rep_feat.get('asm_len') or rep_feat.get('_bp') or 0),
            'rep_len': int(rep_feat.get('source_seq_len') or rep_feat.get('asm_len') or 0),
            'rep_cpg_n': rep_feat.get('cpg_n'),
            'rep_gc_n': rep_feat.get('gc_n'),
            'rep_pct_gc': rep_feat.get('pct_gc'),
            'rep_oe': rep_feat.get('oe'),
            'rep_asm_mid0': fid_mid0(rep_fid),
            'rep_graph_nodeints': nodeints_compact(rep_feat.get('nodeints', [])),
            'locus_primary_source': locus_primary.get(locus_id, {}).get('source', ''),
            'locus_primary_candidate_fid': locus_primary.get(locus_id, {}).get('candidate_fid', ''),
            'locus_primary_chr': locus_primary.get(locus_id, {}).get('chrom', ''),
            'locus_primary_start_median1': locus_primary.get(locus_id, {}).get('start_median1', ''),
            'locus_primary_mid_median1': locus_primary.get(locus_id, {}).get('mid_median1', ''),
            'locus_primary_end_median1': locus_primary.get(locus_id, {}).get('end_median1', ''),
            'locus_nonref_site_source': locus_row.get('nonref_site_source', ''),
            'locus_nonref_site_chrom': locus_row.get('nonref_site_chrom', ''),
            'locus_nonref_site_start1': locus_row.get('nonref_site_start1', ''),
            'locus_nonref_site_end1': locus_row.get('nonref_site_end1', ''),
            'rep_primary_source': rep_primary_source,
            'rep_primary_candidate_fid': rep_primary_candidate_fid,
            'rep_primary_type': rep_primary_type,
            'rep_primary_chr': rep_primary_chr,
            'rep_primary_start1': rep_primary_start1,
            'rep_primary_mid1': rep_primary_mid1,
            'rep_primary_end1': rep_primary_end1,
            'rep_primary_anchor1': rep_primary_anchor1,
            'rep_sv_overlap_class': rep_sv_overlap_class,
            'rep_sv_ins_n': rep_sv_ins_n,
            'rep_sv_ids': rep_sv_ids,
            'rep_sv_primary_sites': rep_sv_primary_sites,
            'rep_sv_primary_intervals': rep_sv_primary_intervals,
            'rep_sv_asm_intervals': rep_sv_asm_intervals,
            'rep_sv_source': rep_sv_source,
            'rep_sv_longest_id': rep_sv_longest_id,
            'rep_sv_longest_chrom': rep_sv_longest_chrom,
            'rep_sv_longest_site1': rep_sv_longest_site1,
            'rep_sv_longest_start1': rep_sv_longest_start1,
            'rep_sv_longest_end1': rep_sv_longest_end1,
            'rep_sv_longest_len': rep_sv_longest_len,
            'allele_member_ids': member_ids_compact(allele_member_ids),
            'allele_member_labels': member_labels_compact(allele_labels_all),
            'allele_member_graph_positions': member_graph_positions_compact(allele_member_ids, feat_idx),
            'allele_member_asm_midpoints': member_asm_midpoints_compact(allele_member_ids),
            'allele_member_primary_midpoints': member_primary_midpoints_compact(allele_member_ids, feat_idx),
            'locus_member_ids': member_ids_compact(locus_member_ids),
            'locus_member_labels': member_labels_compact(locus_labels_all),
            'locus_member_graph_positions': member_graph_positions_compact(locus_member_ids, feat_idx),
            'locus_member_asm_midpoints': member_asm_midpoints_compact(locus_member_ids),
            'locus_member_primary_midpoints': member_primary_midpoints_compact(locus_member_ids, feat_idx),
            'strict_genotype_mode': 'multi_scale_graph_genotyping',
            'allele_callable_asm_n': allele_callable_asm_n,
            'allele_callable_total_n': allele_callable_total_n,
            'locus_callable_asm_n': locus_callable_asm_n,
            'locus_callable_total_n': locus_callable_total_n,
            'allele_gt_1_asm_n': allele_gt_1_asm_n,
            'allele_gt_0_asm_n': allele_gt_0_asm_n,
            'allele_gt_na_asm_n': allele_gt_na_asm_n,
            'allele_gt_1_total_n': allele_gt_1_total_n,
            'allele_gt_0_total_n': allele_gt_0_total_n,
            'allele_gt_na_total_n': allele_gt_na_total_n,
            'locus_gt_1_asm_n': locus_gt_1_asm_n,
            'locus_gt_0_asm_n': locus_gt_0_asm_n,
            'locus_gt_na_asm_n': locus_gt_na_asm_n,
            'locus_gt_1_total_n': locus_gt_1_total_n,
            'locus_gt_0_total_n': locus_gt_0_total_n,
            'locus_gt_na_total_n': locus_gt_na_total_n,
        }
        row.update({
            'locus_primary_cpgi_id': '' if str(locus_row.get('anchor_ref_fid', '')) in ('nan', 'None') else str(locus_row.get('anchor_ref_fid', '')),
            'merged_cpgi_n': allele_member_n_total,
            'merged_cpgi_ids': member_ids_compact(allele_member_ids),
            'locus_cpgi_n': locus_member_n_total,
            'locus_cpgi_ids': member_ids_compact(locus_member_ids),
            'contains_primary_member': 1 if primary_ref_fids else 0,
            'contains_comparison_member': 1 if comparison_ref_fids else 0,
            'allele_freq': '' if allele_freq is None else f'{allele_freq:.6f}',
            'rep_primary_primary_source': rep_primary_primary_source,
            'position_annotation_source': rep_primary_primary_source,
            'rep_primary_primary_chr': rep_primary_primary_chr,
            'rep_primary_primary_start1': rep_primary_primary_start1,
            'rep_primary_primary_end1': rep_primary_primary_end1,
            'rep_primary_primary_site1': rep_primary_primary_site1,
            'rep_hal_multimap_label': rep_feat.get('hal_multimap_label', ''),
            'rep_hal_admit_rule': rep_feat.get('hal_admit_rule', ''),
            'rep_hal_multimap_n': rep_feat.get('hal_multimap_n', ''),
            'rep_hal_primary_chr': rep_feat.get('hal_primary_chr', ''),
            'rep_hal_primary_start0': rep_feat.get('hal_primary_start0', ''),
            'rep_hal_primary_end0': rep_feat.get('hal_primary_end0', ''),
            'rep_hal_primary_start1': rep_feat.get('hal_primary_start1', ''),
            'rep_hal_primary_end1': rep_feat.get('hal_primary_end1', ''),
            'rep_hal_query_cov': rep_feat.get('hal_query_cov', ''),
            'rep_hal_identity': rep_feat.get('hal_identity', ''),
            'rep_hal_all_intervals0': rep_feat.get('hal_all_intervals0', ''),
            'rep_hal_all_intervals': rep_feat.get('hal_all_intervals', ''),
            'rep_sv_all_n': rep_sv_all_n,
            'rep_sv_all_ids': rep_sv_all_ids,
            'rep_sv_all_types': rep_sv_all_types,
            'rep_sv_all_classes': rep_sv_all_classes,
            'rep_sv_all_primary_sites': rep_sv_all_primary_sites,
            'rep_sv_all_primary_intervals': rep_sv_all_primary_intervals,
            'rep_sv_all_asm_intervals': rep_sv_all_asm_intervals,
            'rep_sv_all_overlap_bp': rep_sv_all_overlap_bp,
            'rep_sv_all_overlap_pct_cpgi': rep_sv_all_overlap_pct_cpgi,
            'rep_sv_all_overlap_pct_sv': rep_sv_all_overlap_pct_sv,
            'rep_sv_all_TR': rep_sv_all_TR,
            'rep_sv_all_CONFORMATION': rep_sv_all_CONFORMATION,
            'rep_sv_all_SD': rep_sv_all_SD,
            'rep_sv_all_ITYPE_N': rep_sv_all_ITYPE_N,
            'rep_sv_all_DTYPE_N': rep_sv_all_DTYPE_N,
            'rep_sv_all_FAM_N': rep_sv_all_FAM_N,
            'rep_sv_all_source': rep_sv_all_source,
        })
        allele_label_set = set(allele_labels_all)
        for lab in label_order:
            row[lab] = str(allele_gt_row[lab])
        out_w.writerow([row.get(c, '') for c in out_header])
        big_rows.append(row)
    out_fh.close()

    from pancgi_merged_schema import merged_dataframe
    big_df = merged_dataframe(big_rows, out_header, label_order)
    if bool(getattr(args, 'write_parquet', False)) or str(getattr(args, 'out_parquet', '') or ''):
        pq_path = str(getattr(args, 'out_parquet', '') or '')
        if not pq_path:
            if args.out_tsv.endswith('.tsv.gz'):
                pq_path = args.out_tsv[:-7] + '.parquet'
            elif args.out_tsv.endswith('.tsv'):
                pq_path = args.out_tsv[:-4] + '.parquet'
            else:
                pq_path = args.out_tsv + '.parquet'
        write_optional_parquet(big_df, pq_path)

    if args.out_fasta:
        emit_merged_fasta(args.out_fasta, al, feat_idx)

    print(json_pretty({'merged_rows': len(big_rows), 'out_tsv': os.path.abspath(args.out_tsv)}), file=sys.stderr)


_ARTIFACT_FEAT_IDX: Optional[Dict[str, Dict[str, object]]] = None
_ARTIFACT_SH_WEIGHTS: Optional[Dict[str, float]] = None
_ARTIFACT_ARGS: Optional[argparse.Namespace] = None
_ARTIFACT_ALLELE_MEMBERS: Optional[Dict[str, List[Dict[str, object]]]] = None


def _export_one_locus_artifact(task: Tuple[str, List[str], str]) -> Dict[str, object]:
    global _ARTIFACT_FEAT_IDX, _ARTIFACT_SH_WEIGHTS, _ARTIFACT_ARGS, _ARTIFACT_ALLELE_MEMBERS
    if _ARTIFACT_FEAT_IDX is None or _ARTIFACT_ARGS is None:
        raise RuntimeError('artifact worker globals not initialized')
    feat_idx = _ARTIFACT_FEAT_IDX
    sh_weights = _ARTIFACT_SH_WEIGHTS
    args = _ARTIFACT_ARGS
    locus_id, members, medoid_fid = task
    if int(getattr(args, 'min_members', 0) or 0) and len(members) < int(getattr(args, 'min_members', 0)):
        return {'locus_id': locus_id, 'member_n': len(members), 'status': 'skipped_min_members'}
    if int(getattr(args, 'max_members', 0) or 0) and len(members) > int(getattr(args, 'max_members', 0)):
        return {'locus_id': locus_id, 'member_n': len(members), 'status': 'skipped_max_members'}
    medoid = feat_idx[medoid_fid]
    norm_seq: Dict[str, Optional[str]] = {}
    orient: Dict[str, int] = {}
    ctx_pair: Dict[str, Tuple[float, float]] = {}
    for fid in members:
        f = feat_idx[fid]
        sign, fwd, rev = determine_orientation_to_medoid(f, medoid, anchor_k=args.anchor_k, sh_weights=sh_weights)
        orient[fid] = sign
        ctx_pair[fid] = (fwd, rev)
        seq = f.get('source_seq')
        norm_seq[fid] = reverse_complement(seq) if sign == -1 else seq
    allele_member_rows = []
    for aid, rows in (_ARTIFACT_ALLELE_MEMBERS or {}).items():
        for r in rows:
            if str(r.get('locus_id')) == locus_id:
                fid = str(r['fid'])
                allele_member_rows.append({
                    'allele_id': str(r['allele_id']),
                    'fid': fid,
                    'label': feat_idx[fid]['label'],
                    'kind': feat_idx[fid]['kind'],
                    'is_allele_rep': int(r.get('is_allele_rep', 0)),
                    'orientation_norm': orient.get(fid, 1),
                    'seq_similarity': r.get('seq_similarity', ''),
                    'len_ratio': r.get('len_ratio', ''),
                })
    locus_dir = os.path.join(args.out_dir, f'locus_{locus_id}')
    os.makedirs(locus_dir, exist_ok=True)
    msa_records: List[Tuple[str, str]] = []
    with open(os.path.join(locus_dir, 'normalized.fa'), 'wt') as fh:
        for fid in members:
            seq = norm_seq.get(fid)
            if not seq:
                continue
            head = f"{fid}|label={feat_idx[fid]['label']}|kind={feat_idx[fid]['kind']}"
            msa_records.append((head, str(seq)))
            write_fasta_record(fh, head, str(seq))
    with open(os.path.join(locus_dir, 'allele_clusters.tsv'), 'wt') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['allele_id', 'fid', 'label', 'kind', 'is_allele_rep', 'orientation_norm', 'seq_similarity', 'len_ratio'])
        for mr in allele_member_rows:
            w.writerow([mr['allele_id'], mr['fid'], mr['label'], mr['kind'], mr['is_allele_rep'], mr['orientation_norm'], mr['seq_similarity'], mr['len_ratio']])
    pairwise_written = 0
    if bool(getattr(args, 'pairwise_matrix', False)):
        with gzip.open(os.path.join(locus_dir, 'pairwise_similarity.tsv.gz'), 'wt') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['fid1', 'fid2', 'seq_similarity'])
            for i, fid1 in enumerate(members):
                for j in range(i, len(members)):
                    fid2 = members[j]
                    sim = seq_similarity(norm_seq.get(fid1), norm_seq.get(fid2), args=args)
                    w.writerow([fid1, fid2, '' if sim is None else f"{float(sim):.6f}"])
                    pairwise_written += 1
    msa_path = maybe_write_msa(locus_dir, msa_records, args)
    return {'locus_id': locus_id, 'member_n': len(members), 'pairwise_n': pairwise_written, 'msa_written': 1 if msa_path else 0, 'status': 'ok'}


def cmd_export_locus_artifacts_prod(args: argparse.Namespace) -> None:
    feat_idx = load_prepped_feature_index(args.features, anchor_k=args.anchor_k, max_mid_anchors=args.max_mid_anchors)
    n_features, df = collect_shingle_df(args.features, args.anchor_k, args.max_mid_anchors)
    sh_weights = shingle_weights(n_features, df)
    loc = read_tsv(args.locus_catalog)
    mem = read_tsv(args.locus_members)
    al_mem = read_tsv(args.allele_members)
    loc['locus_id'] = loc['locus_id'].astype(str)
    mem['locus_id'] = mem['locus_id'].astype(str)
    al_mem['locus_id'] = al_mem['locus_id'].astype(str)
    al_mem['allele_id'] = al_mem['allele_id'].astype(str)
    medoid_by_locus = dict(zip(loc['locus_id'], loc['medoid_fid']))
    members_by_locus = {str(k): list(v['fid'].astype(str)) for k, v in mem.groupby('locus_id', sort=True)}
    allele_members_by_id = group_members(al_mem, 'allele_id')
    tasks = [(locus_id, members_by_locus[locus_id], str(medoid_by_locus[locus_id])) for locus_id in sorted(members_by_locus, key=lambda x: int(x) if str(x).isdigit() else str(x))]
    global _ARTIFACT_FEAT_IDX, _ARTIFACT_SH_WEIGHTS, _ARTIFACT_ARGS, _ARTIFACT_ALLELE_MEMBERS
    _ARTIFACT_FEAT_IDX = feat_idx
    _ARTIFACT_SH_WEIGHTS = sh_weights
    _ARTIFACT_ARGS = args
    _ARTIFACT_ALLELE_MEMBERS = allele_members_by_id
    rows: List[Dict[str, object]] = []
    n_threads = int(getattr(args, 'threads', 1))
    if n_threads > 1:
        backend = resolve_parallel_backend(args, default='process')
        chunksize = max(1, int(getattr(args, 'chunksize', 1)))
        if backend == 'thread':
            with ThreadPoolExecutor(max_workers=n_threads) as ex:
                for row in ex.map(_export_one_locus_artifact, tasks):
                    rows.append(row)
        else:
            ctx = mp.get_context(resolve_mp_start_method(args, default='fork'))
            with ctx.Pool(processes=n_threads, maxtasksperchild=int(getattr(args, 'maxtasksperchild', 0) or 0) or None) as pool:
                for row in pool.imap(_export_one_locus_artifact, tasks, chunksize=chunksize):
                    rows.append(row)
    else:
        for task in tasks:
            rows.append(_export_one_locus_artifact(task))
    summary_df = pd.DataFrame(rows)
    summary_tsv = os.path.join(args.out_dir, 'artifact_export_summary.tsv.gz')
    ensure_parent(summary_tsv)
    summary_df.to_csv(summary_tsv, sep='	', compression='gzip', index=False)
    if bool(getattr(args, 'write_parquet', False)):
        write_optional_parquet(summary_df, os.path.join(args.out_dir, 'artifact_export_summary.parquet'))
    print(json_pretty({'locus_dirs': int(len(rows)), 'ok_n': int((summary_df['status'] == 'ok').sum() if len(summary_df) else 0), 'out_dir': os.path.abspath(args.out_dir)}), file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description='PanCGI graph-native non-redundant CpG-island production pipeline')
    sub = parser.add_subparsers(dest='cmd', required=True)

    graph_unfold.add_cli_parser(sub)

    p = sub.add_parser('make-cpgi-fasta', help='Extract per-sample CpGI FASTA directly from HAL')
    p.add_argument('--hal', required=True)
    p.add_argument('--catalog', required=True)
    p.add_argument('--contigs', required=True)
    p.add_argument('--log-dir', required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--docker-bin', default='docker')
    p.add_argument('--docker-image', required=True)
    p.add_argument('--hal2fasta', default='/opt/hal/bin/hal2fasta')
    p.add_argument('--missing-report', required=False, default='')
    p.set_defaults(func=cmd_make_cpgi_fasta)

    p = sub.add_parser('build-features-prod', help='Build callable feature JSONL from BED + pathBED + CpGI FASTA + HAL/SV coordinate sidecars')
    p.add_argument('--threads', type=int, default=1, help='Concurrent genome workers; output preserves catalogue order')
    p.add_argument('--catalog', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--excluded', required=True)
    p.add_argument('--min-graph-cov', type=float, default=0.95)
    p.add_argument('--flank-bp', type=int, default=1000, help='Source-path flank window (bp) to cache on each side of every CpGI feature for strict genotype')
    p.add_argument('--flank-max-steps', type=int, default=32, help='Maximum number of source-path flank steps to cache on each side')
    p.add_argument('--hal-psl-dir', required=True)
    p.add_argument('--hal-min-coverage', type=float, default=0.5, help='Minimum CpGI query coverage for accepting a HAL primary coordinate')
    p.add_argument('--hal-min-identity', type=float, default=0.0, help='Minimum HAL identity for accepting a primary coordinate')
    p.add_argument('--hal-ambig-identity-delta', type=float, default=0.001, help='Multi-hit identity tie window that marks a CpGI ambiguous')
    p.add_argument('--hal-ambig-coverage-delta', type=float, default=0.01, help='Multi-hit query-coverage tie window that marks a CpGI ambiguous')
    p.add_argument('--hal-ambig-aligned-bp-delta', type=int, default=10, help='Multi-hit aligned-bp tie window that marks a CpGI ambiguous')
    p.add_argument('--sv-contig-coordinate-base', choices=[0], type=int, default=0, help='Assembly-contig coordinates are 0-based half-open')
    p.set_defaults(func=cmd_build_features_prod)

    p = sub.add_parser('cluster-loci-prod', help='First-pass graph-native locus clustering')
    p.add_argument('--features', required=True)
    p.add_argument('--out-locus', required=True)
    p.add_argument('--out-members', required=True)
    p.add_argument('--anchor-k', type=int, default=3)
    p.add_argument('--max-mid-anchors', type=int, default=8)
    p.add_argument('--max-anchor-bucket', type=int, default=128)
    p.add_argument('--exact-ro', type=float, default=0.95)
    p.add_argument('--exact-size', type=float, default=0.95)
    p.add_argument('--exact-ctx', type=float, default=0.80)
    p.add_argument('--ro-min', type=float, default=0.50)
    p.add_argument('--ro-ctx', type=float, default=0.30)
    p.add_argument('--szro-ro', type=float, default=0.20)
    p.add_argument('--szro-size', type=float, default=0.67)
    p.add_argument('--szro-ctx', type=float, default=0.50)
    p.add_argument('--weak-seq-gate', type=float, default=0.25)
    p.add_argument('--ambiguity-margin', type=float, default=0.25)
    p.add_argument('--kmer', type=int, default=9)
    p.set_defaults(func=cmd_cluster_loci_prod)

    p = sub.add_parser('polish-loci', help='Choose locus medoids and reassign all features against medoid reps')
    p.add_argument('--features', required=True)
    p.add_argument('--in-locus', required=True)
    p.add_argument('--in-members', required=True)
    p.add_argument('--out-locus', required=True)
    p.add_argument('--out-members', required=True)
    p.add_argument('--anchor-k', type=int, default=3)
    p.add_argument('--max-mid-anchors', type=int, default=8)
    p.add_argument('--max-anchor-bucket', type=int, default=128)
    p.add_argument('--max-medoid-members', type=int, default=500)
    p.add_argument('--exact-ro', type=float, default=0.95)
    p.add_argument('--exact-size', type=float, default=0.95)
    p.add_argument('--exact-ctx', type=float, default=0.80)
    p.add_argument('--ro-min', type=float, default=0.50)
    p.add_argument('--ro-ctx', type=float, default=0.30)
    p.add_argument('--szro-ro', type=float, default=0.20)
    p.add_argument('--szro-size', type=float, default=0.67)
    p.add_argument('--szro-ctx', type=float, default=0.50)
    p.add_argument('--weak-seq-gate', type=float, default=0.25)
    p.add_argument('--ambiguity-margin', type=float, default=0.25)
    p.add_argument('--kmer', type=int, default=9)
    p.add_argument('--threads', type=int, default=1, help='Parallelize medoid selection across provisional loci')
    p.add_argument('--chunksize', type=int, default=8)
    p.add_argument('--maxtasksperchild', type=int, default=0)
    p.add_argument('--feature-store-db', default='', help='Optional sqlite cache path for exact polishing backend')
    p.add_argument('--feature-store-rebuild', action='store_true', help='Force rebuild of the exact polishing SQLite cache')
    p.add_argument('--feature-store-batch-size', type=int, default=4096)
    p.add_argument('--feature-store-read-chunk', type=int, default=900)
    p.set_defaults(func=cmd_polish_loci)


    p = sub.add_parser('anchor-loci-primary', help='Assign polished loci using the primary reference, formal HAL projections and provided SV sites')
    p.add_argument('--anchor-cache-bytes', type=int, default=64 * 1024 * 1024)
    p.add_argument('--anchor-max-record-bytes', type=int, default=16 * 1024 * 1024)
    p.add_argument('--anchor-max-pending', type=int, default=16)
    p.add_argument('--anchor-group-bytes', type=int, default=32 * 1024 * 1024)
    p.add_argument('--anchor-weight-cache-bytes', type=int, default=8 * 1024 * 1024)
    p.add_argument('--anchor-resource-report', default='')
    p.add_argument('--features', required=True)
    p.add_argument('--in-locus', required=True)
    p.add_argument('--in-members', required=True)
    p.add_argument('--out-locus', required=True)
    p.add_argument('--out-members', required=True)
    p.add_argument('--anchor-k', type=int, default=3)
    p.add_argument('--max-mid-anchors', type=int, default=8)
    p.add_argument('--max-medoid-members', type=int, default=500)
    p.add_argument('--nonref-site-window-bp', type=int, default=100, help='Cluster non-ref insertion-only loci on PRIMARY when projected insertion intervals are within this many bp')
    p.add_argument('--hal-nonref-min-coverage', type=float, default=0.5, help='Minimum unique HAL CpGI coverage required for HAL-only and DEL non-ref site clustering')
    p.add_argument('--threads', type=int, default=1, help='Parallelize anchor-locus assignment across features')
    p.add_argument('--parallel-backend', choices=['thread', 'process'], default='process')
    p.add_argument('--mp-start-method', choices=['fork', 'spawn', 'forkserver'], default='fork')
    p.add_argument('--chunksize', type=int, default=8)
    p.add_argument('--maxtasksperchild', type=int, default=0)
    p.set_defaults(func=cmd_anchor_loci_primary)

    p = sub.add_parser('cluster-alleles-prod', help='Allele clustering inside polished loci using whole-CpGI source sequence')
    p.add_argument('--features', required=True)
    p.add_argument('--locus-catalog', required=True)
    p.add_argument('--locus-members', required=True)
    p.add_argument('--out-allele', required=True)
    p.add_argument('--out-allele-members', required=True)
    p.add_argument('--anchor-k', type=int, default=3)
    p.add_argument('--max-mid-anchors', type=int, default=8)
    p.add_argument('--identity', type=float, default=0.80)
    p.add_argument('--min-len-ratio', type=float, default=0.80)
    p.add_argument('--seq-backend', choices=['parasail'], default='parasail')
    p.add_argument('--very-long-threshold', type=int, default=100000, help='If either sequence length reaches this threshold, switch to the very-long backend')
    p.add_argument('--very-long-backend', choices=['external'], default='external', help='Declared WFA backend for very long sequence pairs')
    p.add_argument('--very-long-external-template', default='', help='Shell command template for very long pairs, e.g. "python wfa_longalign_wrapper.py --seq1 {seq1} --seq2 {seq2}". Must emit one JSON object with numeric identity between zero and one.')
    p.add_argument('--parasail-mode', choices=['sg', 'sg_qb', 'sg_qe', 'sg_qx', 'sg_db', 'sg_de', 'sg_dx', 'sg_qb_de', 'sg_qe_db', 'sg_qb_db', 'sg_qe_de'], default='sg')
    p.add_argument('--parasail-match', type=int, default=2)
    p.add_argument('--parasail-mismatch', type=int, default=-3)
    p.add_argument('--parasail-gap-open', type=int, default=5)
    p.add_argument('--parasail-gap-extend', type=int, default=2)
    p.add_argument('--allele-rep-strategy', choices=['asm_seq_medoid', 'seq_medoid', 'ref_first', 'asm_longest', 'longest', 'first'], default='asm_seq_medoid')
    p.add_argument('--allele-cluster-mode', choices=['greedy', 'greedy_clique', 'allpairs_clique'], default='allpairs_clique', help='greedy = representative-based threshold clustering; greedy_clique = greedy seed then all-pairs clique refine inside each seed cluster; allpairs_clique = direct all-pairs clique-constrained clustering within each locus')
    p.add_argument('--threads', type=int, default=1, help='Parallelize allele clustering across loci')
    p.add_argument('--parallel-backend', choices=['thread', 'process'], default='process', help='Use process/fork or threads for concurrent allele tasks')
    p.add_argument('--mp-start-method', choices=['fork', 'spawn', 'forkserver'], default='fork')
    p.add_argument('--chunksize', type=int, default=8)
    p.add_argument('--maxtasksperchild', type=int, default=0, help='0 means unlimited')
    p.add_argument('--max-allele-medoid-members', type=int, default=0, help='0 means no cap; if >0, seq_medoid uses at most this many members when picking the representative.')
    p.add_argument('--dump-locus-dir', default='')
    p.add_argument('--dump-min-members', type=int, default=0)
    p.add_argument('--dump-pairwise-matrix', action='store_true')
    p.add_argument('--dump-msa-backend', choices=['none', 'center_star', 'mafft'], default='none')
    p.add_argument('--dump-msa-name', default='center_star_msa.fa')
    p.add_argument('--dump-msa-min-members', type=int, default=2)
    p.add_argument('--dump-msa-max-members', type=int, default=0, help='0 means no cap')
    p.add_argument('--msa-mafft-bin', default='mafft')
    p.add_argument('--msa-match', type=float, default=2.0)
    p.add_argument('--msa-mismatch', type=float, default=-3.0)
    p.add_argument('--msa-gap-open', type=float, default=-5.0)
    p.add_argument('--msa-gap-extend', type=float, default=-2.0)
    p.add_argument('--flush-every-loci', type=int, default=100, help='Flush output gzip streams every N loci; 0 disables periodic flushing')
    p.add_argument('--locus-ids-file', default='', help='Run only locus IDs listed in this file. Header with locus_id is accepted; otherwise first column is used.')
    p.add_argument('--locus-ids', default='', help='Comma/newline separated locus IDs to run')
    p.add_argument('--locus-start-index', type=int, default=None, help='Optional 0-based start index after explicit locus filtering')
    p.add_argument('--locus-end-index', type=int, default=None, help='Optional 0-based half-open end index after explicit locus filtering')
    p.add_argument('--locus-shard-count', type=int, default=0, help='If >0, run one deterministic shard of the selected locus list')
    p.add_argument('--locus-shard-index', type=int, default=0, help='0-based shard index used with --locus-shard-count')
    p.add_argument('--locus-shard-mode', choices=['contiguous', 'modulo'], default='contiguous', help='contiguous enables simple range shards; modulo can improve load balance')
    p.add_argument('--write-locus-list', default='', help='Write the actual selected locus task list before clustering')
    p.add_argument('--locus-list-only', action='store_true', help='Only write/list selected loci; do not cluster')
    p.add_argument('--feature-load-mode', choices=['auto', 'jsonl', 'sqlite'], default='auto', help='auto uses sqlite for locus-subset/sharded runs and jsonl for full runs')
    p.add_argument('--feature-store-db', default='', help='Allele feature-store sqlite path; default is <features>.step07_allele.sqlite when sqlite mode is used')
    p.add_argument('--feature-store-rebuild', action='store_true', help='Rebuild the allele feature store even if a matching store exists')
    p.add_argument('--feature-store-batch-size', type=int, default=4096)
    p.add_argument('--feature-store-read-chunk', type=int, default=900)
    p.set_defaults(func=cmd_cluster_alleles_prod)

    p = sub.add_parser('build-allele-feature-store', help='Build/reuse the SQLite feature store used by sharded allele clustering')
    p.add_argument('--features', required=True)
    p.add_argument('--out-db', default='', help='Output sqlite path; default is <features>.step07_allele.sqlite')
    p.add_argument('--anchor-k', type=int, default=3)
    p.add_argument('--max-mid-anchors', type=int, default=8)
    p.add_argument('--batch-size', type=int, default=4096)
    p.add_argument('--rebuild', action='store_true')
    p.set_defaults(func=cmd_build_allele_feature_store)

    p = sub.add_parser('export-allele-loci', help='Export sorted locus IDs for Allele clustering and optionally split into shard locus lists')
    p.add_argument('--locus-catalog', required=True)
    p.add_argument('--locus-members', required=True)
    p.add_argument('--out', default='', help='Output all-loci TSV')
    p.add_argument('--shard-count', type=int, default=0)
    p.add_argument('--shard-dir', default='')
    p.add_argument('--shard-prefix', default='allele_loci')
    p.add_argument('--shard-mode', choices=['contiguous', 'modulo'], default='contiguous')
    p.set_defaults(func=cmd_export_allele_loci)

    p = sub.add_parser('export-locus-chunks', help='Split post-anchor locus catalog/members into deterministic chunks for parallel downstream locus jobs')
    p.add_argument('--features', default='', help='Full features.jsonl.gz; required only with --write-feature-subsets')
    p.add_argument('--locus-catalog', required=True)
    p.add_argument('--locus-members', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--prefix', default='chunk')
    p.add_argument('--chunk-count', type=int, default=0)
    p.add_argument('--chunk-size', type=int, default=0)
    p.add_argument('--chunk-mode', choices=['contiguous', 'modulo'], default='contiguous')
    p.add_argument('--write-feature-subsets', action='store_true', help='Also write per-chunk features.jsonl.gz by streaming full features once. Exact sharding should still prefer the global sqlite feature store for shingle weights.')
    p.add_argument('--max-open-files', type=int, default=64)
    p.set_defaults(func=cmd_export_locus_chunks)

    p = sub.add_parser('merge-allele-shards', help='Merge Allele shard outputs back into globally sorted allele_catalog/allele_members tables')
    p.add_argument('--locus-list', required=True, help='Full sorted locus list produced by export-allele-loci')
    p.add_argument('--allele-shards', nargs='+', required=True, help='Shard allele_catalog TSV/GZ files or shell/glob patterns')
    p.add_argument('--allele-member-shards', nargs='+', required=True, help='Shard allele_members TSV/GZ files or shell/glob patterns')
    p.add_argument('--out-allele', required=True)
    p.add_argument('--out-allele-members', required=True)
    p.add_argument('--allow-missing-loci', action='store_true')
    p.set_defaults(func=cmd_merge_allele_shards)

    p = sub.add_parser('strict-genotype-prod', help='Compute strict locus-specific 0/1/NA genotype matrices using pathBED sequence presence')
    p.add_argument('--features', required=True)
    p.add_argument('--locus-catalog', required=True)
    p.add_argument('--locus-members', required=True)
    p.add_argument('--allele-catalog', required=True)
    p.add_argument('--allele-members', required=True)
    p.add_argument('--catalog', required=True)
    p.add_argument('--out-prefix', required=True)
    p.add_argument('--threads', type=int, default=1)
    p.set_defaults(func=cmd_strict_genotype_prod)

    p = sub.add_parser('emit-merged', help='Emit final big merged allele table + merged FASTA')
    p.add_argument('--features', required=True)
    p.add_argument('--locus-catalog', required=True)
    p.add_argument('--locus-members', required=True)
    p.add_argument('--allele-catalog', required=True)
    p.add_argument('--allele-members', required=True)
    p.add_argument('--catalog', required=False, default='', help='Optional catalog to control haplotype column order in the merged genotype table')
    p.add_argument('--locus-genotype-tsv', required=True)
    p.add_argument('--allele-genotype-tsv', required=True)
    p.add_argument('--out-tsv', required=True)
    p.add_argument('--out-fasta', required=False, default='')
    p.add_argument('--out-parquet', required=False, default='')
    p.add_argument('--write-parquet', action='store_true')
    p.set_defaults(func=cmd_emit_merged)

    p = sub.add_parser('export-locus-artifacts-prod', help='Post-process per-locus normalized FASTA, cluster membership, pairwise similarity and optional MSA artifacts from completed results')
    p.add_argument('--features', required=True)
    p.add_argument('--locus-catalog', required=True)
    p.add_argument('--locus-members', required=True)
    p.add_argument('--allele-members', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--anchor-k', type=int, default=3)
    p.add_argument('--max-mid-anchors', type=int, default=8)
    p.add_argument('--min-members', type=int, default=0)
    p.add_argument('--max-members', type=int, default=0)
    p.add_argument('--pairwise-matrix', action='store_true')
    p.add_argument('--msa-backend', dest='dump_msa_backend', choices=['none', 'center_star', 'mafft'], default='center_star')
    p.add_argument('--dump-msa-name', default='center_star_msa.fa')
    p.add_argument('--dump-msa-min-members', type=int, default=2)
    p.add_argument('--dump-msa-max-members', type=int, default=0)
    p.add_argument('--msa-mafft-bin', default='mafft')
    p.add_argument('--msa-match', type=float, default=2.0)
    p.add_argument('--msa-mismatch', type=float, default=-3.0)
    p.add_argument('--msa-gap-open', type=float, default=-5.0)
    p.add_argument('--msa-gap-extend', type=float, default=-2.0)
    p.add_argument('--seq-backend', choices=['parasail'], default='parasail')
    p.add_argument('--very-long-threshold', type=int, default=100000)
    p.add_argument('--very-long-backend', choices=['external'], default='external')
    p.add_argument('--very-long-external-template', default='')
    p.add_argument('--parasail-mode', choices=['sg', 'sg_qb', 'sg_qe', 'sg_qx', 'sg_db', 'sg_de', 'sg_dx', 'sg_qb_de', 'sg_qe_db', 'sg_qb_db', 'sg_qe_de'], default='sg')
    p.add_argument('--parasail-match', type=int, default=2)
    p.add_argument('--parasail-mismatch', type=int, default=-3)
    p.add_argument('--parasail-gap-open', type=int, default=5)
    p.add_argument('--parasail-gap-extend', type=int, default=2)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--parallel-backend', choices=['thread', 'process'], default='process')
    p.add_argument('--mp-start-method', choices=['fork', 'spawn', 'forkserver'], default='fork')
    p.add_argument('--chunksize', type=int, default=1)
    p.add_argument('--maxtasksperchild', type=int, default=0)
    p.add_argument('--write-parquet', action='store_true')
    p.set_defaults(func=cmd_export_locus_artifacts_prod)

    parser.allow_abbrev = False
    for command in sub.choices.values():
        command.allow_abbrev = False
    args = parser.parse_args()
    if os.environ.get('PANCGI_PARAMETERS_DIR'):
        directory = Path(os.environ['PANCGI_PARAMETERS_DIR'])
        directory.mkdir(parents=True, exist_ok=True)
        parameters = {k:v for k,v in vars(args).items() if not callable(v)}
        (directory/(args.cmd+'.json')).write_text(json.dumps(parameters, indent=2)+'\n')
    args.func(args)


if __name__ == '__main__':
    main()
