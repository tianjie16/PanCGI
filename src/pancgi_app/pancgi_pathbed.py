import bisect
import csv
import os
import re
import threading
from collections import OrderedDict, defaultdict, deque
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import pancgi_core as core

WALK_TOKEN_RE = re.compile(r'([<>])([^<>,]+)')
PATHBED_TOKEN_RE = re.compile(r'^([<>])(.+)$')


def iter_pathbed(path: str) -> Iterator[Tuple[str, int, int, str]]:
    with core.open_text(path, 'rt') as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith('#'):
                continue
            f = raw.rstrip('\n').split('\t')
            if len(f) < 4:
                f = raw.rstrip('\n').split()
            if len(f) < 4:
                raise ValueError(f'Malformed pathBED row in {path}: {raw[:200]!r}')
            yield f[0], int(f[1]), int(f[2]), f[3]


_PATHBED_FILE_INDEX_CACHE: Dict[str, Dict[Tuple[str, str], str]] = {}


def _load_pathbed_file_index(path_bed_dir: str) -> Dict[Tuple[str, str], str]:
    root = os.path.abspath(path_bed_dir)
    cached = _PATHBED_FILE_INDEX_CACHE.get(root)
    if cached is not None:
        return cached
    inventory = os.path.join(root, 'pathbed_inventory.tsv')
    if not os.path.isfile(inventory):
        raise FileNotFoundError(f'PathBED inventory is required: {inventory}')
    index: Dict[Tuple[str, str], str] = {}
    with open(inventory, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'genome_id', 'contig_id', 'output_file'}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f'PathBED inventory is missing columns: {missing}')
        for line_number, row in enumerate(reader, start=2):
            key = (str(row['genome_id']).strip(), str(row['contig_id']).strip())
            output_file = str(row['output_file']).strip()
            if not all(key) or not output_file or key in index:
                raise ValueError(f'Invalid or duplicate pathBED inventory row at line {line_number}: {key}')
            resolved = os.path.join(root, output_file)
            if not os.path.isfile(resolved):
                raise FileNotFoundError(resolved)
            index[key] = resolved
    _PATHBED_FILE_INDEX_CACHE[root] = index
    return index


def resolve_pathbed_file(path_bed_dir: str, sample: str, hap: str, contig: str) -> Optional[str]:
    if str(hap) != '0':
        raise ValueError(f'Internal graph haplotype key must be 0 after manifest normalization: {hap!r}')
    return _load_pathbed_file_index(path_bed_dir).get((str(sample), str(contig)))


def iter_pathbed_files_for_genome(path_bed_dir: str, genome_id: str, internal_hap: str) -> Iterator[str]:
    if str(internal_hap) != '0':
        raise ValueError(f'Internal graph haplotype key must be 0 after manifest normalization: {internal_hap!r}')
    index = _load_pathbed_file_index(path_bed_dir)
    for (_genome_id, contig_id), path in sorted(index.items()):
        if _genome_id == str(genome_id):
            yield path


_PATHBED_ROWS_CACHE: "OrderedDict[str, Dict[str, object]]" = OrderedDict()
_PATHBED_ROWS_CACHE_LOCK = threading.Lock()
_PATHBED_ROWS_CACHE_MAX_FILES = 16


def set_pathbed_cache_max_files(n: int) -> None:
    global _PATHBED_ROWS_CACHE_MAX_FILES
    _PATHBED_ROWS_CACHE_MAX_FILES = max(1, int(n))
    with _PATHBED_ROWS_CACHE_LOCK:
        while len(_PATHBED_ROWS_CACHE) > _PATHBED_ROWS_CACHE_MAX_FILES:
            _PATHBED_ROWS_CACHE.popitem(last=False)


def _build_pathbed_cache_entry(pathbed_path: str) -> Dict[str, object]:
    rows = list(iter_pathbed(pathbed_path))
    starts = [int(r[1]) for r in rows]
    ends = [int(r[2]) for r in rows]
    return {
        'rows': rows,
        'starts': starts,
        'ends': ends,
    }


def load_pathbed_rows_indexed_cached(pathbed_path: str, max_cached_files: int = _PATHBED_ROWS_CACHE_MAX_FILES) -> Dict[str, object]:
    pathbed_path = str(pathbed_path)
    with _PATHBED_ROWS_CACHE_LOCK:
        entry = _PATHBED_ROWS_CACHE.get(pathbed_path)
        if entry is not None:
            _PATHBED_ROWS_CACHE.move_to_end(pathbed_path)
            if isinstance(entry, dict):
                return entry
            upgraded = _build_pathbed_cache_entry(pathbed_path)
            _PATHBED_ROWS_CACHE[pathbed_path] = upgraded
            _PATHBED_ROWS_CACHE.move_to_end(pathbed_path)
            return upgraded
    entry = _build_pathbed_cache_entry(pathbed_path)
    with _PATHBED_ROWS_CACHE_LOCK:
        _PATHBED_ROWS_CACHE[pathbed_path] = entry
        _PATHBED_ROWS_CACHE.move_to_end(pathbed_path)
        while len(_PATHBED_ROWS_CACHE) > max_cached_files:
            _PATHBED_ROWS_CACHE.popitem(last=False)
    return entry


def load_pathbed_rows_cached(pathbed_path: str, max_cached_files: int = _PATHBED_ROWS_CACHE_MAX_FILES) -> List[Tuple[str, int, int, str]]:
    return load_pathbed_rows_indexed_cached(pathbed_path, max_cached_files=max_cached_files)['rows']


def _token_to_step(token: str, row_s: int, row_e: int, seg_s: int, seg_e: int) -> Optional[List[int]]:
    if seg_s >= seg_e:
        return None
    m = PATHBED_TOKEN_RE.match(token)
    if not m:
        raise ValueError(f'Cannot parse path-bed token: {token}')
    ori = m.group(1)
    node_id = int(m.group(2))
    node_len = row_e - row_s
    if node_len <= 0:
        return None
    local_s = seg_s - row_s
    local_e = seg_e - row_s
    if ori == '>':
        node_s, node_e, rev = local_s, local_e, 0
    else:
        node_s, node_e, rev = node_len - local_e, node_len - local_s, 1
    return [node_id, rev, node_s, node_e]


def _segments_to_flank_steps(
    segments: List[Tuple[int, int, int, int, str]],
    side: str,
    flank_bp: int,
    flank_max_steps: int,
) -> Tuple[List[List[int]], int]:
    if not segments:
        return [], 0
    picked: List[Tuple[int, int, int, int, str]] = []
    bp = 0
    it = reversed(segments) if side == 'left' else iter(segments)
    for seg_s, seg_e, row_s, row_e, tok in it:
        if seg_s >= seg_e:
            continue
        seg_len = max(0, seg_e - seg_s)
        if flank_bp > 0 and bp + seg_len > flank_bp:
            keep = max(0, flank_bp - bp)
            if keep <= 0:
                break
            if side == 'left':
                seg_s = seg_e - keep
            else:
                seg_e = seg_s + keep
            seg_len = keep
        picked.append((seg_s, seg_e, row_s, row_e, tok))
        bp += seg_len
        if (flank_bp > 0 and bp >= flank_bp) or (flank_max_steps > 0 and len(picked) >= flank_max_steps):
            break
    if side == 'left':
        picked.reverse()
    out: List[List[int]] = []
    for seg_s, seg_e, row_s, row_e, tok in picked:
        step = _token_to_step(tok, row_s, row_e, seg_s, seg_e)
        if step is not None:
            out.append(step)
    return out, bp


def _map_single_feature_on_rows_indexed(
    contig: str,
    start0: int,
    end0: int,
    fid: str,
    entry: Dict[str, object],
    flank_bp: int = 1000,
    flank_max_steps: int = 32,
) -> Optional[Dict[str, object]]:
    rows: List[Tuple[str, int, int, str]] = list(entry.get('rows', []))
    starts: List[int] = list(entry.get('starts', []))
    ends: List[int] = list(entry.get('ends', []))
    if not rows:
        return None

    left_idx = bisect.bisect_right(ends, int(start0))
    right_excl = bisect.bisect_left(starts, int(end0))
    if right_excl < left_idx:
        right_excl = left_idx

    body_steps: List[List[int]] = []
    graph_bp = 0
    graph_path: List[str] = []
    for i in range(left_idx, min(right_excl, len(rows))):
        row_contig, row_s, row_e, token = rows[i]
        if row_contig != contig:
            continue
        ov_s = max(int(start0), int(row_s))
        ov_e = min(int(end0), int(row_e))
        if ov_s >= ov_e:
            continue
        step = _token_to_step(token, int(row_s), int(row_e), ov_s, ov_e)
        if step is not None:
            body_steps.append(step)
            graph_bp += int(ov_e) - int(ov_s)
            graph_path.append(token)

    if not body_steps:
        return {
            'steps': [],
            'nodeints': [],
            'graph_bp': 0,
            'graph_cov': 0.0,
            'graph_path': '',
            'graph_nsteps': 0,
            'left_flank_steps': [],
            'right_flank_steps': [],
            'left_flank_bp': 0,
            'right_flank_bp': 0,
            'left_flank_nsteps': 0,
            'right_flank_nsteps': 0,
        }

    left_segments: List[Tuple[int, int, int, int, str]] = []
    remain_bp = int(flank_bp)
    if left_idx < len(rows):
        row_contig, row_s, row_e, token = rows[left_idx]
        if row_contig == contig and int(start0) > int(row_s):
            seg_e = int(start0)
            seg_s = int(row_s)
            if remain_bp > 0 and seg_e - seg_s > remain_bp:
                seg_s = seg_e - remain_bp
            if seg_s < seg_e:
                left_segments.append((seg_s, seg_e, int(row_s), int(row_e), token))
                if remain_bp > 0:
                    remain_bp -= seg_e - seg_s
    i = left_idx - 1
    while i >= 0 and ((flank_max_steps <= 0 or len(left_segments) < flank_max_steps) and (flank_bp <= 0 or remain_bp > 0)):
        row_contig, row_s, row_e, token = rows[i]
        if row_contig != contig:
            i -= 1
            continue
        seg_s = int(row_s)
        seg_e = int(row_e)
        if flank_bp > 0 and seg_e - seg_s > remain_bp:
            seg_s = seg_e - remain_bp
        if seg_s < seg_e:
            left_segments.append((seg_s, seg_e, int(row_s), int(row_e), token))
            if flank_bp > 0:
                remain_bp -= seg_e - seg_s
        i -= 1
    left_segments.reverse()
    left_flank_steps, left_flank_bp = _segments_to_flank_steps(left_segments, side='left', flank_bp=flank_bp, flank_max_steps=flank_max_steps)

    right_segments: List[Tuple[int, int, int, int, str]] = []
    remain_bp = int(flank_bp)
    last_body_idx = min(right_excl, len(rows)) - 1
    if 0 <= last_body_idx < len(rows):
        row_contig, row_s, row_e, token = rows[last_body_idx]
        if row_contig == contig and int(end0) < int(row_e):
            seg_s = int(end0)
            seg_e = int(row_e)
            if remain_bp > 0 and seg_e - seg_s > remain_bp:
                seg_e = seg_s + remain_bp
            if seg_s < seg_e:
                right_segments.append((seg_s, seg_e, int(row_s), int(row_e), token))
                if remain_bp > 0:
                    remain_bp -= seg_e - seg_s
    i = right_excl
    while i < len(rows) and ((flank_max_steps <= 0 or len(right_segments) < flank_max_steps) and (flank_bp <= 0 or remain_bp > 0)):
        row_contig, row_s, row_e, token = rows[i]
        if row_contig != contig:
            i += 1
            continue
        seg_s = int(row_s)
        seg_e = int(row_e)
        if flank_bp > 0 and seg_e - seg_s > remain_bp:
            seg_e = seg_s + remain_bp
        if seg_s < seg_e:
            right_segments.append((seg_s, seg_e, int(row_s), int(row_e), token))
            if flank_bp > 0:
                remain_bp -= seg_e - seg_s
        i += 1
    right_flank_steps, right_flank_bp = _segments_to_flank_steps(right_segments, side='right', flank_bp=flank_bp, flank_max_steps=flank_max_steps)

    feat_len = int(end0) - int(start0)
    graph_cov = (graph_bp / feat_len) if feat_len > 0 else 0.0
    return {
        'steps': body_steps,
        'nodeints': core.merge_nodeints_from_steps(body_steps) if body_steps else [],
        'graph_bp': int(graph_bp),
        'graph_cov': graph_cov,
        'graph_path': ''.join(graph_path),
        'graph_nsteps': len(body_steps),
        'left_flank_steps': left_flank_steps,
        'right_flank_steps': right_flank_steps,
        'left_flank_bp': int(left_flank_bp),
        'right_flank_bp': int(right_flank_bp),
        'left_flank_nsteps': len(left_flank_steps),
        'right_flank_nsteps': len(right_flank_steps),
    }


def _map_features_on_rows(
    contig: str,
    feature_rows: List[Tuple[int, int, str]],
    rows: Sequence[Tuple[str, int, int, str]],
    flank_bp: int = 1000,
    flank_max_steps: int = 32,
) -> Dict[str, Dict[str, object]]:

    feature_rows = sorted(feature_rows)
    active: List[Tuple[int, int, str]] = []
    next_idx = 0
    n_feat = len(feature_rows)
    out_steps: Dict[str, List[List[int]]] = defaultdict(list)
    out_bp: Dict[str, int] = defaultdict(int)
    out_path: Dict[str, List[str]] = defaultdict(list)

    left_flank_steps: Dict[str, List[List[int]]] = {}
    left_flank_bp: Dict[str, int] = defaultdict(int)
    right_flank_steps: Dict[str, List[List[int]]] = defaultdict(list)
    right_flank_bp: Dict[str, int] = defaultdict(int)

    prev_segments: deque = deque()
    prev_bp = 0
    pending_right: Dict[str, Dict[str, object]] = {}

    for row_contig, row_s, row_e, token in rows:
        if row_contig != contig:
            continue

        while next_idx < n_feat and feature_rows[next_idx][0] < row_e:
            fs, fe, fid = feature_rows[next_idx]
            if fe > row_s:
                segs = list(prev_segments)
                if fs > row_s:
                    segs.append((row_s, fs, row_s, row_e, token))
                lf_steps, lf_bp = _segments_to_flank_steps(segs, side='left', flank_bp=flank_bp, flank_max_steps=flank_max_steps)
                left_flank_steps[fid] = lf_steps
                left_flank_bp[fid] = lf_bp
                active.append((fs, fe, fid))
            next_idx += 1

        ended_in_row: set = set()
        if active:
            m = PATHBED_TOKEN_RE.match(token)
            if not m:
                raise ValueError(f'Cannot parse path-bed token: {token}')
            ori = m.group(1)
            node_id = int(m.group(2))
            node_len = row_e - row_s
            if node_len > 0:
                for fs, fe, fid in active:
                    ov_s = max(fs, row_s)
                    ov_e = min(fe, row_e)
                    if ov_s < ov_e:
                        local_s = ov_s - row_s
                        local_e = ov_e - row_s
                        if ori == '>':
                            node_s, node_e, rev = local_s, local_e, 0
                        else:
                            node_s, node_e, rev = node_len - local_e, node_len - local_s, 1
                        out_steps[fid].append([node_id, rev, node_s, node_e])
                        out_bp[fid] += (ov_e - ov_s)
                        out_path[fid].append(token)
                    if fe <= row_e:
                        ended_in_row.add(fid)
                        if fid not in pending_right:
                            pending_right[fid] = {'segments': [], 'bp': 0}
                        if fe < row_e:
                            pending_right[fid]['segments'].append((fe, row_e, row_s, row_e, token))
                            pending_right[fid]['bp'] += max(0, row_e - fe)
                if ended_in_row:
                    active = [x for x in active if x[2] not in ended_in_row]

        if pending_right:
            done = []
            for fid, state in pending_right.items():
                if fid not in ended_in_row:
                    state['segments'].append((row_s, row_e, row_s, row_e, token))
                    state['bp'] += max(0, row_e - row_s)
                if (flank_bp > 0 and state['bp'] >= flank_bp) or (flank_max_steps > 0 and len(state['segments']) >= flank_max_steps):
                    rf_steps, rf_bp = _segments_to_flank_steps(state['segments'], side='right', flank_bp=flank_bp, flank_max_steps=flank_max_steps)
                    right_flank_steps[fid] = rf_steps
                    right_flank_bp[fid] = rf_bp
                    done.append(fid)
            for fid in done:
                del pending_right[fid]

        prev_segments.append((row_s, row_e, row_s, row_e, token))
        prev_bp += max(0, row_e - row_s)
        while (flank_max_steps > 0 and len(prev_segments) > flank_max_steps) or (flank_bp > 0 and prev_bp > flank_bp and len(prev_segments) > 1):
            seg_s, seg_e, _row_s0, _row_e0, _tok = prev_segments.popleft()
            prev_bp -= max(0, seg_e - seg_s)

    for fid, state in pending_right.items():
        rf_steps, rf_bp = _segments_to_flank_steps(state['segments'], side='right', flank_bp=flank_bp, flank_max_steps=flank_max_steps)
        right_flank_steps[fid] = rf_steps
        right_flank_bp[fid] = rf_bp

    mapped: Dict[str, Dict[str, object]] = {}
    for fs, fe, fid in feature_rows:
        steps = out_steps.get(fid, [])
        graph_bp = out_bp.get(fid, 0)
        feat_len = fe - fs
        graph_cov = graph_bp / feat_len if feat_len > 0 else 0.0
        lfs = left_flank_steps.get(fid, [])
        rfs = right_flank_steps.get(fid, [])
        mapped[fid] = {
            'steps': steps,
            'nodeints': core.merge_nodeints_from_steps(steps) if steps else [],
            'graph_bp': graph_bp,
            'graph_cov': graph_cov,
            'graph_path': ''.join(out_path.get(fid, [])),
            'graph_nsteps': len(steps),
            'left_flank_steps': lfs,
            'right_flank_steps': rfs,
            'left_flank_bp': int(left_flank_bp.get(fid, sum(max(0, int(s[3]) - int(s[2])) for s in lfs))),
            'right_flank_bp': int(right_flank_bp.get(fid, sum(max(0, int(s[3]) - int(s[2])) for s in rfs))),
            'left_flank_nsteps': len(lfs),
            'right_flank_nsteps': len(rfs),
        }
    return mapped


def map_features_on_contig(
    contig: str,
    feature_rows: List[Tuple[int, int, str]],
    pathbed_path: str,
    flank_bp: int = 1000,
    flank_max_steps: int = 32,
) -> Dict[str, Dict[str, object]]:


    rows = iter_pathbed(pathbed_path)
    return _map_features_on_rows(contig, feature_rows, rows, flank_bp=flank_bp, flank_max_steps=flank_max_steps)


def map_single_feature_on_contig_cached(
    contig: str,
    start0: int,
    end0: int,
    fid: str,
    pathbed_path: str,
    flank_bp: int = 1000,
    flank_max_steps: int = 32,
    max_cached_files: int = _PATHBED_ROWS_CACHE_MAX_FILES,
) -> Optional[Dict[str, object]]:
    entry = load_pathbed_rows_indexed_cached(pathbed_path, max_cached_files=max_cached_files)
    return _map_single_feature_on_rows_indexed(
        contig,
        int(start0),
        int(end0),
        str(fid),
        entry,
        flank_bp=flank_bp,
        flank_max_steps=flank_max_steps,
    )
