from __future__ import annotations

import gzip
import re
from collections import defaultdict
from typing import Dict, List

FID_RE = re.compile(
    r"^(?P<sample>[^#]+)#(?P<hap>[^#]+)#(?P<contig>.+):"
    r"(?P<start>\d+)-(?P<end>\d+)$"
)


def open_text(path: str, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_fid(fid: str) -> Dict[str, object]:
    match = FID_RE.match(str(fid))
    if not match:
        raise ValueError(f"Cannot parse feature ID: {fid}")
    return {
        "sample": match.group("sample"),
        "hap": match.group("hap"),
        "contig": match.group("contig"),
        "start0": int(match.group("start")),
        "end0": int(match.group("end")),
    }


def merge_nodeints_from_steps(steps: List[List[int]]) -> List[List[int]]:
    per_node: Dict[int, List[List[int]]] = defaultdict(list)
    for node_id, _reverse, start, end in steps:
        per_node[int(node_id)].append([int(start), int(end)])
    merged: List[List[int]] = []
    for node_id, intervals in per_node.items():
        intervals.sort()
        current_start, current_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append([node_id, current_start, current_end])
                current_start, current_end = start, end
        merged.append([node_id, current_start, current_end])
    merged.sort(key=lambda item: (item[0], item[1], item[2]))
    return merged


def kind_priority(kind: str) -> int:
    if kind == "reference_primary":
        return 3
    if kind == "reference_comparison":
        return 2
    return 1
