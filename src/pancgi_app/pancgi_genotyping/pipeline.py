import csv
import gzip
import hashlib
import json
import multiprocessing
import sqlite3
import tempfile
import zlib
from array import array
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

import pancgi_core as core
import pancgi_pathbed as pathbed
from pancgi_genotype import insertion_envelope
from .caller import CONFIG, evaluate, propagate
from .exact import model_exact_evidence


def rows(path):
    with core.open_text(str(path)) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def read_walk(path, contig):
    walk = (array('q'), array('q'), array('q'), array('b'))
    previous = 0
    for name, start, end, token in pathbed.iter_pathbed(str(path)):
        if name != contig or start != previous or end <= start or token[0] not in '<>':
            raise ValueError('Source path is not a complete contiguous walk')
        walk[0].append(start)
        walk[1].append(end)
        walk[2].append(int(token[1:]))
        walk[3].append(int(token[0] == '<'))
        previous = end
    if not previous:
        raise ValueError('Source path is empty')
    return walk


def interval_steps(walk, start, end):
    starts, ends, nodes, reverse = walk
    first = bisect_right(ends, start)
    last = bisect_left(starts, end)
    result = []
    for i in range(first, last):
        lo, hi = max(start, starts[i]), min(end, ends[i])
        if hi <= lo:
            continue
        offset0, offset1 = (lo - starts[i], hi - starts[i]) if not reverse[i] else (ends[i] - hi, ends[i] - lo)
        result.append([int(nodes[i]), int(reverse[i]), int(offset0), int(offset1)])
    return result


def build_models(args, catalog, loci, alleles, destination):
    import cpgi_nr_prod as prod
    by_locus, wanted = defaultdict(list), set()
    for allele in alleles:
        by_locus[allele['locus_id']].append(allele['allele_rep_fid'])
        wanted.add(allele['allele_rep_fid'])
    for locus in loci:
        if locus['anchor_type'] == 'primary_ref':
            wanted.add(locus['anchor_ref_fid'])
    if str(args.features).endswith('.sqlite'):
        from pancgi_features import selected_raw
        features = selected_raw(args.features, wanted)
    else:
        features = {r['fid']: r for r in prod.iter_feature_jsonl(args.features) if r['fid'] in wanted}
    if set(features) != wanted:
        raise ValueError('Genotype source features are missing')
    entries = {r['label']: r for r in catalog}
    primary = [r for r in catalog if r['kind'] == 'reference_primary']
    if len(primary) != 1:
        raise ValueError('Exactly one primary reference is required')
    groups = defaultdict(list)
    model_counts = {}
    for locus in loci:
        lid = locus['locus_id']
        specs = []
        if locus['anchor_type'] == 'primary_ref':
            specs.append(('primary_reference_anchor', locus['anchor_ref_fid']))
        specs.extend(('allele_representative', fid) for fid in by_locus[lid])
        seen, rank = set(), 0
        for role, fid in specs:
            if fid in seen:
                continue
            seen.add(fid)
            f = features[fid]
            pf = core.parse_fid(fid)
            entry = entries[f['label']]
            source = pathbed.resolve_pathbed_file(entry['path_bed_dir'], entry['graph_sample'], entry['graph_hap'], pf['contig'])
            if source is None:
                raise FileNotFoundError(fid)
            start, end, n = insertion_envelope(f)
            model = dict(model_role=role, model_fid=fid, steps=f['steps'], feat_len=int(f['source_seq_len']),
                         envelope_start0=start, envelope_end0=end, envelope_len=end-start, insertion_envelope_n=n)
            groups[(source, pf['contig'])].append((lid, rank, model, (pf['start0'], pf['end0'])))
            rank += 1
        if locus['anchor_type'] != 'primary_ref':
            chrom = locus['nonref_site_chrom']
            start, end = int(locus['nonref_site_start1']) - 1, int(locus['nonref_site_end1'])
            entry = primary[0]
            source = pathbed.resolve_pathbed_file(entry['path_bed_dir'], entry['graph_sample'], entry['graph_hap'], chrom)
            if source is None:
                raise FileNotFoundError(f'Primary reference context: {lid}')
            model = dict(model_role='main_reference_context', model_fid='main_reference_context:' + lid,
                         feat_len=end-start, envelope_start0=start, envelope_end0=end,
                         envelope_len=end-start, insertion_envelope_n=0)
            groups[(source, chrom)].append((lid, rank, model, (start, end)))
            rank += 1
        if not rank:
            raise ValueError(f'No source models for {lid}')
        model_counts[lid] = rank
    connection = sqlite3.connect(destination)
    connection.execute('CREATE TABLE models (locus TEXT, rank INT, payload BLOB, PRIMARY KEY(locus,rank))')
    try:
        for (source, contig), tasks in sorted(groups.items()):
            walk = read_walk(source, contig)
            for lid, rank, model, body in tasks:
                start, end = model['envelope_start0'], model['envelope_end0']
                if not 0 <= start < end <= walk[1][-1]:
                    raise ValueError('Source envelope is outside the path')
                body_steps = interval_steps(walk, *body)
                if 'steps' in model and model['steps'] != body_steps:
                    raise ValueError('Genotype body changed relative to the feature')
                model['steps'] = body_steps
                for side, lo, hi in [('left', max(0, start-CONFIG['source_context_bp']), start),
                                     ('right', end, min(walk[1][-1], end+CONFIG['source_context_bp']))]:
                    steps = interval_steps(walk, lo, hi)
                    bp = sum(s[3]-s[2] for s in steps)
                    if bp != hi-lo:
                        raise ValueError('Incomplete source context')
                    model[side + '_flank_steps'] = steps
                    model[side + '_flank_len'] = bp
                model['expanded_flank_bp'] = CONFIG['source_context_bp']
                model['expanded_flank_max_steps'] = None
                connection.execute('INSERT INTO models VALUES (?,?,?)', (lid, rank, zlib.compress(json.dumps(model, separators=(',', ':')).encode(), 1)))
            connection.commit()
        actual = dict(connection.execute('SELECT locus,COUNT(*) FROM models GROUP BY locus'))
        if actual != model_counts:
            raise RuntimeError('Incomplete source model inventory')
    finally:
        connection.close()


def build_target(entry, destination):
    files = list(pathbed.iter_pathbed_files_for_genome(entry['path_bed_dir'], entry['graph_sample'], entry['graph_hap']))
    if not files:
        raise FileNotFoundError('No target paths')
    connection = sqlite3.connect(destination)
    connection.execute('PRAGMA cache_size=-32768')
    connection.execute('CREATE TABLE occurrences (node INT, chrom TEXT, start0 INT, end0 INT, rev INT, ordinal INT)')
    batch, seen = [], set()
    for file in files:
        previous, name = 0, None
        for ordinal, (chrom, start, end, token) in enumerate(pathbed.iter_pathbed(file)):
            if start != previous or end <= start or (name is not None and name != chrom) or token[0] not in '<>':
                raise ValueError('Incomplete or malformed target walk')
            name, previous = chrom, end
            batch.append((int(token[1:]), chrom, start, end, int(token[0] == '<'), ordinal))
            if len(batch) >= 100000:
                connection.executemany('INSERT INTO occurrences VALUES (?,?,?,?,?,?)', batch)
                connection.commit()
                batch.clear()
        if name is None or name in seen:
            raise ValueError('Empty or duplicated target contig')
        seen.add(name)
    connection.executemany('INSERT INTO occurrences VALUES (?,?,?,?,?,?)', batch)
    connection.execute('CREATE INDEX node_lookup ON occurrences(node)')
    connection.commit()
    return connection


def target_evidence(connection, models):
    wanted = sorted({int(step[0]) for m in models for key in ('steps', 'left_flank_steps', 'right_flank_steps') for step in m[key]})
    by_node, by_key = defaultdict(list), {}
    for offset in range(0, len(wanted), 500):
        batch = wanted[offset:offset+500]
        for node, chrom, start, end, reverse, ordinal in connection.execute(
                'SELECT * FROM occurrences WHERE node IN (' + ','.join('?' for _ in batch) + ') ORDER BY node,chrom,ordinal,rev,start0,end0', batch):
            key = node, chrom, ordinal, reverse
            if key in by_key:
                raise ValueError('Duplicate graph occurrence identity')
            value = (chrom, start, end, reverse, ordinal)
            by_node[node].append(value)
            by_key[key] = value
    target = dict(by_node=by_node, by_key=by_key, occurrence_count_by_node={n: len(v) for n, v in by_node.items()})
    result = []
    for rank, m in enumerate(models, 1):
        evidence = model_exact_evidence(m, target, max_distance_bp=CONFIG['source_context_bp'], max_occurrences_per_step=0, minimum_run_bp=1)
        if not evidence['exact_complete']:
            raise RuntimeError('Incomplete exact evidence')
        evidence.update(model_rank=rank, left_flank_len=m['left_flank_len'], right_flank_len=m['right_flank_len'])
        result.append(evidence)
    return dict(candidate_evidence=result)


def run_target(task):
    entry, model_file, locus_ids, allele_map, observed_loci, observed_alleles, output = task
    output = Path(output)
    with tempfile.TemporaryDirectory(dir=output, prefix=entry['label'] + '.') as scratch:
        target = None
        source = sqlite3.connect(f'file:{model_file}?mode=ro', uri=True)
        loci, alleles = {}, {}
        try:
            with gzip.open(output / (entry['label'] + '.evidence.jsonl.gz'), 'wt') as log:
                for lid in locus_ids:
                    observed = set(allele_map[lid]) & observed_alleles
                    if lid in observed_loci:
                        decision = dict(call='0', reason='observed_member')
                    else:
                        if target is None:
                            target = build_target(entry, Path(scratch) / 'target.sqlite')
                        models = [json.loads(zlib.decompress(r[0])) for r in source.execute('SELECT payload FROM models WHERE locus=? ORDER BY rank', (lid,))]
                        if not models:
                            raise ValueError('Missing locus model')
                        decision = evaluate(target_evidence(target, models), scratch)
                    loci[lid], values = propagate(lid in observed_loci, observed, allele_map[lid], decision['call'])
                    alleles.update(values)
                    log.write(json.dumps(dict(locus_id=lid, hal_genome=entry['hal_genome'], genotype=loci[lid], evidence=decision), separators=(',', ':')) + '\n')
        finally:
            if target is not None:
                target.close()
            source.close()
    return entry['label'], loci, alleles


def run(args):
    import cpgi_nr_prod as prod
    if args.threads < 1:
        raise ValueError('Worker count must be positive')
    catalog = list(prod.iter_internal_catalog(args.catalog))
    loci, alleles = rows(args.locus_catalog), rows(args.allele_catalog)
    locus_ids, allele_ids = [r['locus_id'] for r in loci], [r['allele_id'] for r in alleles]
    if len(set(locus_ids)) != len(locus_ids) or len(set(allele_ids)) != len(allele_ids):
        raise ValueError('Duplicated locus or allele IDs')
    allele_map = {lid: [] for lid in locus_ids}
    for row in alleles:
        allele_map[row['locus_id']].append(row['allele_id'])
    observed_loci, observed_alleles = defaultdict(set), defaultdict(set)
    allele_to_locus = {r['allele_id']: r['locus_id'] for r in alleles}
    for row in rows(args.locus_members):
        if row['locus_id'] not in allele_map:
            raise ValueError('Unknown locus in member table')
        observed_loci[row['label']].add(row['locus_id'])
    for row in rows(args.allele_members):
        if row['allele_id'] not in allele_to_locus or row['locus_id'] != allele_to_locus[row['allele_id']]:
            raise ValueError('Unknown or incorrectly assigned allele in member table')
        observed_alleles[row['label']].add(row['allele_id'])
    labels = [r['label'] for r in catalog]
    if not (set(observed_loci) | set(observed_alleles)) <= set(labels):
        raise ValueError('Membership table contains undeclared genomes')
    if len(set(labels)) != len(labels):
        raise ValueError('Duplicated genome identity')
    for label in labels:
        if observed_loci[label] != {allele_to_locus[a] for a in observed_alleles[label]}:
            raise ValueError('Locus and allele observed memberships disagree')
    root = Path(args.out_prefix + '.genotyping')
    if root.exists():
        raise FileExistsError(root)
    root.mkdir(parents=True)
    models = root / 'models.sqlite'
    build_models(args, catalog, loci, alleles, models)
    tasks = [(r, str(models.resolve()), locus_ids, allele_map, observed_loci[r['label']], observed_alleles[r['label']], str(root)) for r in catalog]
    results = {}
    if args.threads > 1:
        with multiprocessing.get_context('spawn').Pool(args.threads, maxtasksperchild=1) as pool:
            for label, loc, al in pool.imap_unordered(run_target, tasks):
                results[label] = (loc, al)
    else:
        for task in tasks:
            label, loc, al = run_target(task)
            results[label] = (loc, al)
    if set(results) != set(labels):
        raise RuntimeError('Missing target results; refusing to fill with NA')
    for kind, ids, column in [('locus', locus_ids, 0), ('allele', allele_ids, 1)]:
        filename = args.out_prefix + f'.{kind}_genotype_strict.tsv.gz'
        with gzip.open(filename + '.partial', 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow([kind + '_id'] + labels)
            for identifier in ids:
                writer.writerow([identifier] + [results[label][column][identifier] for label in labels])
        Path(filename + '.partial').replace(filename)
    (root / 'completed.json').write_text(json.dumps(dict(status='pass', config=CONFIG, target_n=len(labels), locus_n=len(locus_ids), allele_n=len(allele_ids)), indent=2) + '\n')
