import argparse
import csv
import json
import shutil
import sys
import tempfile
from pathlib import Path

from pancgi_mapping import ProgressLog
from pancgi_parallel import completed_tasks


def build_task(task):
    import cpgi_nr_prod as prod
    ordinal, row, options, directory = task
    args = argparse.Namespace(**options)
    root = Path(directory)
    feature_path = root / (f'{ordinal:06d}.features' + ('.gz' if args.out.endswith('.gz') else ''))
    excluded_path = root / (f'{ordinal:06d}.excluded' + ('.gz' if args.excluded.endswith('.gz') else ''))
    try:
        with prod.open_text(str(feature_path), 'wt') as out, prod.open_text(str(excluded_path), 'wt') as excluded:
            writer = csv.writer(excluded, delimiter='\t', lineterminator='\n')
            report = prod.build_feature_genome(row, args, out, writer)
        report.update(label=row['label'], feature_path=str(feature_path), excluded_path=str(excluded_path))
        return report
    except Exception as exc:
        raise RuntimeError(f"Feature construction failed for {row['label']}: {exc}") from exc
    finally:
        prod._SV_INS_CACHE.clear()


def run(args):
    import cpgi_nr_prod as prod
    workers = getattr(args, 'threads', 1)
    if type(workers) is not int or workers < 1:
        raise ValueError('Worker count must be a positive integer')
    rows = list(prod.iter_internal_catalog(args.catalog))
    labels = [row['label'] for row in rows]
    if not rows or len(labels) != len(set(labels)):
        raise ValueError('Feature catalogue must contain unique, nonempty genome rows')
    out, excluded = Path(args.out), Path(args.excluded)
    if out.resolve() == excluded.resolve():
        raise ValueError('Feature and exclusion output paths must differ')
    temporary_out = Path(prod.temp_output_path(str(out)))
    temporary_excluded = Path(prod.temp_output_path(str(excluded)))
    progress_path = Path(str(out) + '.progress.jsonl')
    for path in (out, excluded, temporary_out, temporary_excluded, progress_path):
        if path.exists():
            raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
    progress = ProgressLog(progress_path)
    counts = dict(features_written=0, excluded=0, partial_graph_coverage_excluded=0, input_records=0)
    options = {key: value for key, value in vars(args).items() if not callable(value)}
    effective_workers = min(workers, len(rows))
    progress.emit('features', 'start', genomes=len(rows), workers=effective_workers)
    try:
        with tempfile.TemporaryDirectory(prefix='.feature-parts-', dir=out.parent) as directory:
            tasks = [(i, row, options, directory) for i, row in enumerate(rows)]
            if effective_workers == 1:
                results = ((i, build_task(task)) for i, task in enumerate(tasks))
            else:
                results = completed_tasks(build_task, tasks, effective_workers)
            ordered = {}
            try:
                for ordinal, result in results:
                    ordered[ordinal] = result
                    for key in counts:
                        counts[key] += result[key]
                    progress.emit('features', 'genome_complete', ordinal=ordinal, label=result['label'],
                        genomes_completed=len(ordered), genomes_total=len(rows), **counts)
            finally:
                results.close()
            if set(ordered) != set(range(len(rows))):
                raise RuntimeError('Incomplete feature task set')
            with prod.open_text(str(temporary_excluded), 'wt') as handle:
                handle.write('label\tkind\tfid\treason\n')
            with temporary_out.open('wb') as merged, temporary_excluded.open('ab') as rejected:
                for ordinal in range(len(rows)):
                    record = ordered[ordinal]
                    for key, destination in [('feature_path', merged), ('excluded_path', rejected)]:
                        with open(record[key], 'rb') as source:
                            shutil.copyfileobj(source, destination, 1024 * 1024)
                    progress.emit('features', 'genome_merged', ordinal=ordinal, label=record['label'])
            prod.finalize_output_path(str(temporary_out), str(out))
            prod.finalize_output_path(str(temporary_excluded), str(excluded))
        progress.emit('features', 'complete', genomes=len(rows), workers=effective_workers, **counts)
    except BaseException as exc:
        progress.emit('features', 'error', message=str(exc))
        raise
    finally:
        progress.close()
    print(json.dumps(dict(features_written=counts['features_written'], excluded=counts['excluded'],
        partial_graph_coverage_excluded=counts['partial_graph_coverage_excluded'],
        missing_source_sequence=0, coordinate_backend='hal_sv'), indent=2, sort_keys=True), file=sys.stderr)
