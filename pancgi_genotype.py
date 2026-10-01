import json
from typing import Dict, Tuple


def insertion_envelope(feature: Dict[str, object]) -> Tuple[int, int, int]:
    start0 = int(feature.get('asm_start0') or 0)
    end0 = int(feature.get('asm_end0') or start0)
    feature_start0 = start0
    feature_end0 = end0
    contig = str(feature.get('contig') or '')
    raw = feature.get('sv_ins_detail_json')
    if raw in (None, '', '[]'):
        return start0, end0, 0
    try:
        records = json.loads(str(raw))
    except Exception as exc:
        raise ValueError(f'Invalid sv_ins_detail_json for {feature.get("fid", "")!r}') from exc
    used = 0
    for record in records:
        if str(record.get('svtype') or '').upper() != 'INS':
            continue
        if str(record.get('contig') or contig) != contig:
            continue
        ins_start0 = record.get('eff_start0', record.get('asm_start0'))
        ins_end0 = record.get('eff_end0', record.get('asm_end0'))
        if ins_start0 in (None, '') or ins_end0 in (None, ''):
            continue
        ins_start0 = int(ins_start0)
        ins_end0 = int(ins_end0)
        if ins_end0 <= feature_start0 or ins_start0 >= feature_end0:
            continue
        start0 = min(start0, ins_start0)
        end0 = max(end0, ins_end0)
        used += 1
    return start0, end0, used
