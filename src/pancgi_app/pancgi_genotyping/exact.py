from __future__ import annotations

from collections import defaultdict


def step_bp(step):
    return max(0, int(step[3]) - int(step[2]))


def prefix_bp(steps):
    values = [0]
    for step in steps:
        values.append(values[-1] + step_bp(step))
    return values


def _matching_occurrence(lookup, step, chrom, ordinal, direction):
    node_id = int(step[0])
    source_rev = int(step[1])
    target_rev = source_rev if int(direction) == 1 else 1 - source_rev
    return lookup.get((node_id, str(chrom), int(ordinal), int(target_rev)))


def _project_step(step, occurrence):
    node_start = int(step[2])
    node_end = int(step[3])
    target_start = int(occurrence[1])
    target_end = int(occurrence[2])
    target_rev = int(occurrence[3])
    if target_rev == 0:
        return target_start + node_start, target_start + node_end
    return target_end - node_end, target_end - node_start


def _distance_from_envelope(prefix, index, side):
    if side == "right":
        return int(prefix[index])
    if side == "left":
        return int(prefix[-1] - prefix[index + 1])
    return 0


def _eligible_indices(steps, side, max_distance_bp):
    if side == "body":
        return 0, len(steps)
    prefix = prefix_bp(steps)
    indices = []
    for index in range(len(steps)):
        if _distance_from_envelope(prefix, index, side) <= int(max_distance_bp):
            indices.append(index)
    if not indices:
        return 0, 0
    return min(indices), max(indices) + 1


def exhaustive_maximal_exact_runs(
    steps,
    side,
    target_index,
    max_distance_bp=1000000,
    max_occurrences_per_step=0,
    minimum_run_bp=1,
):
    if side not in {"left", "right", "body"}:
        raise ValueError(side)
    normalized = [tuple(int(value) for value in step) for step in steps]
    start_index, end_index = _eligible_indices(normalized, side, max_distance_bp)
    if start_index == end_index:
        return [], {
            "source_step_n": len(normalized),
            "eligible_step_n": 0,
            "present_step_n": 0,
            "saturated_step_n": 0,
            "raw_placement_n": 0,
            "source_run_n": 0,
            "run_n": 0,
            "exact_complete": True,
        }

    by_node = target_index["by_node"]
    lookup = target_index["by_key"]
    occurrence_count_by_node = target_index.get("occurrence_count_by_node", {})
    prefix = prefix_bp(normalized)
    raw = {}
    occurrence_counts = {
        source_index: int(occurrence_count_by_node.get(
            int(normalized[source_index][0]),
            len(by_node.get(int(normalized[source_index][0]), [])),
        ))
        for source_index in range(start_index, end_index)
    }
    present_step_n = sum(count > 0 for count in occurrence_counts.values())
    saturated_indices = {
        source_index
        for source_index, count in occurrence_counts.items()
        if int(max_occurrences_per_step) > 0 and count > int(max_occurrences_per_step)
    }
    saturated_step_n = len(saturated_indices)
    for source_index in range(start_index, end_index):
        step = normalized[source_index]
        occurrences = by_node.get(int(step[0]), [])
        if not occurrences:
            continue
        if source_index in saturated_indices:
            continue
        for occurrence in occurrences:
            chrom = str(occurrence[0])
            direction = 1 if int(step[1]) == int(occurrence[3]) else -1
            ordinal = int(occurrence[4])
            if source_index > start_index and source_index - 1 not in saturated_indices:
                previous = _matching_occurrence(
                    lookup,
                    normalized[source_index - 1],
                    chrom,
                    ordinal - direction,
                    direction,
                )
                if previous is not None:
                    continue
            hi = source_index
            while hi + 1 < end_index:
                following = _matching_occurrence(
                    lookup,
                    normalized[hi + 1],
                    chrom,
                    ordinal + direction * (hi + 1 - source_index),
                    direction,
                )
                if following is None:
                    break
                hi += 1
            run_support_bp = int(prefix[hi + 1] - prefix[source_index])
            if run_support_bp < int(minimum_run_bp):
                continue
            matched = [
                _matching_occurrence(
                    lookup,
                    normalized[index],
                    chrom,
                    ordinal + direction * (index - source_index),
                    direction,
                )
                for index in range(source_index, hi + 1)
            ]
            if any(item is None for item in matched):
                raise RuntimeError("maximal exact run lost an internally verified step")
            projected = [
                _project_step(source_step, target_step)
                for source_step, target_step in zip(normalized[source_index : hi + 1], matched)
            ]
            target_start0 = min(item[0] for item in projected)
            target_end0 = max(item[1] for item in projected)
            ordinals = [int(item[4]) for item in matched]
            source_start = float(prefix[source_index])
            source_end = float(prefix[hi + 1])
            if side == "left":
                source_start -= float(prefix[-1])
                source_end -= float(prefix[-1])
            key = (
                int(source_index),
                int(hi),
                chrom,
                int(direction),
                int(target_start0),
                int(target_end0),
                min(ordinals),
                max(ordinals),
            )
            raw[key] = {
                "source_lo": int(source_index),
                "source_hi": int(hi),
                "source_start_bp": source_start,
                "source_end_bp": source_end,
                "source_center_bp": (source_start + source_end) / 2.0,
                "source_support_bp": run_support_bp,
                "source_step_n": int(hi - source_index + 1),
                "source_distance_bp": _distance_from_envelope(prefix, source_index if side == "right" else hi, side),
                "chrom": chrom,
                "direction": int(direction),
                "target_start1": int(target_start0) + 1,
                "target_end1": int(target_end0),
                "target_center1": (float(target_start0) + float(target_end0) + 1.0) / 2.0,
                "target_ord_start": min(ordinals),
                "target_ord_end": max(ordinals),
            }

    placements_by_source = defaultdict(list)
    for row in raw.values():
        placements_by_source[(int(row["source_lo"]), int(row["source_hi"]))].append(row)
    runs = []
    for run_id, ((_lo, _hi), placements) in enumerate(sorted(placements_by_source.items()), start=1):
        placement_n = len(placements)
        for placement in placements:
            row = dict(placement)
            row["source_run_id"] = f"{side[0].upper()}{run_id}"
            row["placement_n"] = int(placement_n)
            runs.append(row)
    runs.sort(key=lambda row: (
        int(row["source_distance_bp"]),
        int(row["placement_n"]),
        -int(row["source_support_bp"]),
        str(row["chrom"]),
        int(row["target_start1"]),
    ))
    qc = {
        "source_step_n": len(normalized),
        "eligible_step_n": int(end_index - start_index),
        "present_step_n": int(present_step_n),
        "saturated_step_n": int(saturated_step_n),
        "raw_placement_n": len(raw),
        "source_run_n": len(placements_by_source),
        "run_n": len(runs),
        "exact_complete": saturated_step_n == 0,
    }
    return runs, qc


def model_exact_evidence(
    model,
    target_index,
    max_distance_bp=1000000,
    max_occurrences_per_step=0,
    minimum_run_bp=1,
):
    body_runs, body_qc = exhaustive_maximal_exact_runs(
        model.get("steps", []),
        "body",
        target_index,
        max_distance_bp=max_distance_bp,
        max_occurrences_per_step=max_occurrences_per_step,
        minimum_run_bp=minimum_run_bp,
    )
    left_runs, left_qc = exhaustive_maximal_exact_runs(
        model.get("left_flank_steps", []),
        "left",
        target_index,
        max_distance_bp=max_distance_bp,
        max_occurrences_per_step=max_occurrences_per_step,
        minimum_run_bp=minimum_run_bp,
    )
    right_runs, right_qc = exhaustive_maximal_exact_runs(
        model.get("right_flank_steps", []),
        "right",
        target_index,
        max_distance_bp=max_distance_bp,
        max_occurrences_per_step=max_occurrences_per_step,
        minimum_run_bp=minimum_run_bp,
    )
    return {
        "model_fid": str(model.get("model_fid") or ""),
        "model_role": str(model.get("model_role") or ""),
        "feat_len": int(model.get("feat_len") or 0),
        "envelope_len": int(model.get("envelope_len") or model.get("feat_len") or 0),
        "insertion_envelope_n": int(model.get("insertion_envelope_n") or 0),
        "body_runs": body_runs,
        "left_runs": left_runs,
        "right_runs": right_runs,
        "body_qc": body_qc,
        "left_qc": left_qc,
        "right_qc": right_qc,
        "exact_complete": all(qc["exact_complete"] for qc in (body_qc, left_qc, right_qc)),
    }
