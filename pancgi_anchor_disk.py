import codecs
import csv
import json
import marshal
import os
import sqlite3
import threading
import time
import zlib
from collections import OrderedDict
from collections.abc import Mapping
from pathlib import Path

import cpgi_nr_prod as prod
from pancgi_anchor_bounded import deep_bytes, ordered_assignments


def key_text(value):
    return json.dumps(value, separators=(',', ':'), ensure_ascii=True)


def compare_keys(left, right):
    a, b = json.loads(left), json.loads(right)
    return (a > b) - (a < b)


def connect(path, readonly=False):
    uri = Path(path).resolve().as_uri() + ('?mode=ro' if readonly else '?mode=rwc')
    db = sqlite3.connect(uri, uri=True, check_same_thread=False)
    db.create_collation('PYKEY', compare_keys)
    db.execute('PRAGMA cache_size=-4096')
    db.execute('PRAGMA mmap_size=0')
    db.execute('PRAGMA temp_store=FILE')
    if readonly:
        db.execute('PRAGMA query_only=ON')
    return db


class Lookup:
    def __init__(self, path):
        self.path = str(path)
        self.pid = os.getpid()
        self.lock = threading.RLock()
        self.db = connect(path, True)

    def local(self):
        if self.pid != os.getpid():
            self.db.close()
            self.pid = os.getpid()
            self.lock = threading.RLock()
            self.db = connect(self.path, True)

    def one(self, sql, parameters=()):
        self.local()
        with self.lock:
            return self.db.execute(sql, parameters).fetchone()

    def rows(self, sql, parameters=()):
        self.local()
        cursor = self.db.execute(sql, parameters)
        try:
            while True:
                with self.lock:
                    row = cursor.fetchone()
                if row is None:
                    break
                yield row
        finally:
            cursor.close()

    def close(self):
        self.db.close()


class ByteCache:
    def __init__(self, limit):
        self.limit = limit
        self.pid = os.getpid()
        self.data = OrderedDict()
        self.size = self.peak = 0
        self.lock = threading.RLock()

    def get_or_load(self, key, loader):
        if self.pid != os.getpid():
            self.__init__(self.limit)
        with self.lock:
            if key in self.data:
                self.data.move_to_end(key)
                return self.data[key][0]
            value = loader()
            size = deep_bytes((key, value)) + 192
            while self.data and self.size + size > self.limit:
                _, (_, removed) = self.data.popitem(last=False)
                self.size -= removed
            if size <= self.limit:
                self.data[key] = (value, size)
                self.size += size
                self.peak = max(self.peak, self.size)
            return value


class KeyMap(Mapping):
    def __init__(self, lookup, table, cache_bytes=0):
        if table not in ('weights', 'seeds'):
            raise ValueError(table)
        self.lookup, self.table = lookup, table
        self.count = lookup.one('SELECT count(*) FROM ' + table)[0]
        self.cache = ByteCache(cache_bytes)

    def __len__(self):
        return self.count

    def __iter__(self):
        for row in self.lookup.rows('SELECT k FROM ' + self.table + ' ORDER BY k'):
            yield row[0]

    def __getitem__(self, key):
        def read():
            row = self.lookup.one('SELECT v FROM ' + self.table + ' WHERE k=?', (key,))
            if row is None:
                raise KeyError(key)
            return row[0]
        return self.cache.get_or_load(key, read)


class ReferenceMap(Mapping):
    def __init__(self, lookup, features):
        self.lookup, self.features = lookup, features
        self.count = lookup.one('SELECT count(*) FROM refs')[0]

    def __len__(self):
        return self.count

    def __iter__(self):
        for row in self.lookup.rows('SELECT fid FROM refs ORDER BY ordinal'):
            yield row[0]

    def __getitem__(self, fid):
        if self.lookup.one('SELECT 1 FROM refs WHERE fid=?', (fid,)) is None:
            raise KeyError(fid)
        return self.features[fid]


class MetadataMap:

    def __init__(self, lookup, cache_bytes):
        self.lookup = lookup
        self.cache = ByteCache(cache_bytes)

    def __getitem__(self, fid):
        def read():
            row = self.lookup.one('SELECT metadata FROM tasks WHERE fid=?', (fid,))
            if row is None:
                raise KeyError(fid)
            return marshal.loads(row[0])
        return self.cache.get_or_load(fid, read)


class PhaseProgress:
    def __init__(self, path, features):
        self.handle = Path(path).open('x')
        self.features = features
        self.started = time.monotonic()

    def emit(self, phase, event, processed=None, total=None):
        record = dict(phase=phase, event=event, elapsed_seconds=time.monotonic()-self.started,
                      parent_cpu_seconds=time.process_time(), processed=processed, total=total,
                      parent_feature_access=self.features.statistics())
        self.handle.write(json.dumps(record, sort_keys=True) + '\n')
        self.handle.flush()

    def close(self):
        self.handle.close()


def emit(progress, phase, event, processed=None, total=None):
    if progress is not None:
        progress.emit(phase, event, processed, total)


def bounded_list(items, limit, label, metrics=None):
    output, size = [], 64
    for item in items:
        
        size += deep_bytes(item) + 16
        if size > limit:
            raise MemoryError(label + ' exceeds scientific group byte limit; no truncation allowed')
        output.append(item)
    if metrics is not None:
        metrics[label + '_peak_accounted_bytes'] = max(metrics.get(label + '_peak_accounted_bytes', 0), size)
        metrics[label + '_peak_records'] = max(metrics.get(label + '_peak_records', 0), len(output))
    return output


class IntervalMap:
    def __init__(self, lookup, limit):
        self.lookup, self.limit = lookup, limit
        self.cache = ByteCache(limit)

    def get(self, chrom, default=None):
        def read():
            return bounded_list((marshal.loads(row[0]) for row in self.lookup.rows(
                'SELECT payload FROM intervals WHERE chrom=? ORDER BY sk COLLATE PYKEY', (chrom,))),
                self.limit, 'reference chromosome')
        return self.cache.get_or_load(chrom, read)


def compressed_text(blob):
    decoder, utf8 = zlib.decompressobj(), codecs.getincrementaldecoder('utf-8')()
    while True:
        block = blob.read(65536)
        if not block:
            break
        while block:
            raw = decoder.decompress(block, 65536)
            block = decoder.unconsumed_tail
            if decoder.unused_data:
                raise ValueError('Trailing compressed frequency data')
            if raw:
                yield utf8.decode(raw)
    if not decoder.eof:
        raise ValueError('Incomplete compressed frequency data')
    yield utf8.decode(b'', final=True)


def object_pairs(chunks, token_limit):
    iterator, decoder = iter(chunks), json.JSONDecoder()
    buffer, position, ended = '', 0, False

    def fill():
        nonlocal buffer, position, ended
        buffer, position = buffer[position:], 0
        try:
            block = next(iterator)
        except StopIteration:
            ended = True
            return False
        if len(buffer) + len(block) > token_limit:
            raise MemoryError('Frequency JSON token exceeds byte-resource limit')
        buffer += block
        return True

    def peek():
        nonlocal position
        while True:
            while position < len(buffer) and buffer[position].isspace():
                position += 1
            if position < len(buffer):
                return buffer[position]
            if ended or not fill():
                return ''

    def expect(character):
        nonlocal position
        if peek() != character:
            raise ValueError('Malformed frequency JSON; expected ' + character)
        position += 1

    def value():
        nonlocal position
        if not peek():
            raise ValueError('Truncated frequency JSON')
        while True:
            try:
                item, stop = decoder.raw_decode(buffer, position)
            except json.JSONDecodeError:
                if ended or not fill():
                    raise ValueError('Malformed frequency JSON value') from None
                continue
            
            if stop == len(buffer) and not ended:
                fill()
                continue
            position = stop
            return item

    expect('{')
    if peek() != '}':
        while True:
            key = value()
            if not isinstance(key, str):
                raise ValueError('Frequency key must be a string')
            expect(':')
            count = value()
            if type(count) is not int:
                raise ValueError('Frequency must be an integer')
            separator = peek()
            if separator not in (',', '}'):
                raise ValueError('Malformed frequency JSON separator')
            yield key, count
            if separator == '}':
                break
            expect(',')
    expect('}')
    if peek():
        raise ValueError('Trailing frequency JSON')


def build_lookup(path, args, features, metrics, progress=None):
    db = connect(path)
    try:
        db.executescript('''
            CREATE TABLE weights(k TEXT PRIMARY KEY, v REAL NOT NULL) WITHOUT ROWID;
            CREATE TABLE seeds(k TEXT PRIMARY KEY, v TEXT NOT NULL) WITHOUT ROWID;
            CREATE TABLE tasks(fid TEXT PRIMARY KEY, sk TEXT NOT NULL, metadata BLOB NOT NULL) WITHOUT ROWID;
            CREATE TABLE refs(fid TEXT PRIMARY KEY, ordinal INTEGER NOT NULL) WITHOUT ROWID;
            CREATE TABLE intervals(chrom TEXT, sk TEXT, payload BLOB);
        ''')
        emit(progress, 'weights', 'start')
        rows = features.conn.execute('SELECT rowid FROM shared_frequencies').fetchmany(2)
        if len(rows) != 1:
            raise ValueError('Expected one complete-background frequency payload')
        with features.conn.blobopen('shared_frequencies', 'payload', rows[0][0], readonly=True) as blob:
            count = 0
            for key, frequency in object_pairs(compressed_text(blob), args.anchor_max_record_bytes):
                if not 1 <= frequency <= len(features):
                    raise ValueError('Invalid full-background frequency')
                weight = prod.shingle_weights(len(features), {key: frequency})[key]
                db.execute('INSERT INTO weights VALUES (?,?)', (key, weight))
                count += 1
                if count % 10000 == 0:
                    db.commit()
        metrics['full_background_weight_count'] = count
        emit(progress, 'weights', 'complete', count, count)
        emit(progress, 'seeds', 'start')
        with prod.open_text(args.in_members, 'rt') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not {'fid', 'locus_id'} <= set(reader.fieldnames or []):
                raise ValueError('Seed member columns missing')
            for i, row in enumerate(reader, 1):
                if deep_bytes(row) > args.anchor_max_record_bytes:
                    raise MemoryError('Seed row exceeds byte limit')
                db.execute('INSERT OR REPLACE INTO seeds VALUES (?,?)', (str(row['fid']), str(row['locus_id'])))
                if i % 10000 == 0:
                    db.commit()
        emit(progress, 'seeds', 'complete')
        emit(progress, 'metadata', 'start', 0, len(features))
        prepared_before = features.loads
        count = 0
        for ordinal, (fid, rec) in enumerate(features.iter_metadata(), 1):
            db.execute('INSERT INTO tasks VALUES (?,?,?)',
                       (fid, key_text(prod.anchor_task_sort_key(fid)), marshal.dumps(rec)))
            if str(rec.get('kind')) == 'reference_primary':
                db.execute('INSERT INTO refs VALUES (?,?)', (fid, ordinal))
                interval = prod.ref_interval1_from_feat(rec)
                if interval is not None:
                    chrom, start, end = interval
                    entry = (int(start), int(end), str(fid))
                    db.execute('INSERT INTO intervals VALUES (?,?,?)', (str(chrom), key_text(entry), marshal.dumps(entry)))
            if ordinal % 10000 == 0:
                db.commit()
                emit(progress, 'metadata', 'progress', ordinal, len(features))
            count = ordinal
        metrics['metadata_records'] = count
        metrics['lookup_prepared_loads'] = features.loads - prepared_before
        if count != len(features) or metrics['lookup_prepared_loads']:
            raise RuntimeError('Lookup metadata must cover the full background without feature preparation')
        emit(progress, 'metadata', 'complete', count, len(features))
        emit(progress, 'lookup_sort', 'start')
        db.execute('CREATE INDEX task_order ON tasks(sk COLLATE PYKEY)')
        db.execute('CREATE INDEX interval_order ON intervals(chrom, sk COLLATE PYKEY)')
        db.commit()
        metrics['seed_map_count'] = db.execute('SELECT count(*) FROM seeds').fetchone()[0]
        metrics['primary_reference_count'] = db.execute('SELECT count(*) FROM refs').fetchone()[0]
        emit(progress, 'lookup_sort', 'complete')
    finally:
        db.close()


def create_spool(path):
    db = connect(path)
    db.executescript('''
        CREATE TABLE assignments(ordinal INTEGER PRIMARY KEY, anchor_key TEXT, payload BLOB, cluster_token TEXT);
        CREATE TABLE sites(ordinal INTEGER PRIMARY KEY, sk TEXT NOT NULL);
        CREATE TABLE clusters(token TEXT PRIMARY KEY, sk TEXT NOT NULL) WITHOUT ROWID;
        CREATE TABLE loci(anchor_key TEXT PRIMARY KEY, sk TEXT NOT NULL, locus_id INTEGER, payload BLOB) WITHOUT ROWID;
    ''')
    return db


def add_assignment(db, ordinal, row, args):
    if deep_bytes(row) > args.anchor_max_record_bytes:
        raise MemoryError('Assignment row exceeds byte limit')
    db.execute('INSERT INTO assignments VALUES (?,?,?,NULL)', (ordinal, str(row['anchor_key']), marshal.dumps(row)))
    seed = prod.nonref_site_seed_from_row(row, args) if str(row.get('anchor_type')) == 'nonref' else None
    if seed is not None:
        seed_id = row.get('seed_locus_id') or ''
        sk = int(seed_id) if str(seed_id).isdigit() else ordinal
        key = [str(seed['chrom']), int(seed['start1']), int(seed['end1']), sk, ordinal]
        db.execute('INSERT INTO sites VALUES (?,?)', (ordinal, key_text(key)))


def cluster_sites(db, args, metrics):
    db.execute('CREATE INDEX site_order ON sites(sk COLLATE PYKEY)')
    db.commit()
    window = int(getattr(args, 'nonref_site_window_bp', 100) or 0)
    component, accounted, previous_chrom, high = [], 64, None, None
    component_id = 0

    def flush():
        nonlocal component, accounted, component_id
        if not component:
            return
        component_id += 1
        component.sort(key=lambda x: x[0])
        rows = [row for _, row in component]
        prod.cluster_nonref_site_rows(rows, args)
        clusters = {}
        for ordinal, row in component:
            token = str(component_id) + ':' + str(row['nonref_site_cluster_ord'])
            if token not in clusters:
                clusters[token] = [str(row['nonref_site_chrom']), int(row['nonref_site_start1']),
                                   int(row['nonref_site_end1']), ordinal]
            db.execute('UPDATE assignments SET payload=?,cluster_token=? WHERE ordinal=?',
                       (marshal.dumps(row), token, ordinal))
        for token, sk in clusters.items():
            db.execute('INSERT INTO clusters VALUES (?,?)', (token, key_text(sk)))
        metrics['site_component_peak_accounted_bytes'] = max(metrics.get('site_component_peak_accounted_bytes', 0), accounted)
        metrics['site_component_peak_records'] = max(metrics.get('site_component_peak_records', 0), len(component))
        db.commit()
        component, accounted = [], 64

    for ordinal, sk, payload in db.execute('SELECT s.ordinal,s.sk,a.payload FROM sites s JOIN assignments a USING(ordinal) ORDER BY s.sk COLLATE PYKEY'):
        chrom, start, end, _, _ = json.loads(sk)
        
        if component and (chrom != previous_chrom or start > high + window):
            flush()
            high = None
        row = marshal.loads(payload)
        accounted += deep_bytes((ordinal, row)) + 32
        if accounted > args.anchor_group_bytes:
            raise MemoryError('Nonref connected work group exceeds byte limit; no cluster truncation allowed')
        component.append((ordinal, row))
        previous_chrom, high = chrom, max(end, high) if high is not None else end
    flush()
    db.execute('CREATE INDEX cluster_members ON assignments(cluster_token,ordinal)')
    db.execute('CREATE INDEX cluster_order ON clusters(sk COLLATE PYKEY)')
    db.commit()
    for cluster_ord, (token, sk) in enumerate(db.execute('SELECT token,sk FROM clusters ORDER BY sk COLLATE PYKEY'), 1):
        chrom = json.loads(sk)[0]
        key = f'NRSITE::{chrom}::{cluster_ord}'
        for ordinal, payload in db.execute('SELECT ordinal,payload FROM assignments WHERE cluster_token=? ORDER BY ordinal', (token,)):
            row = marshal.loads(payload)
            row.update(anchor_key=key, nonref_site_cluster_ord=cluster_ord)
            db.execute('UPDATE assignments SET anchor_key=?,payload=? WHERE ordinal=?', (key, marshal.dumps(row), ordinal))
        if cluster_ord % 1000 == 0:
            db.commit()
    db.execute('CREATE INDEX assignment_locus ON assignments(anchor_key,ordinal)')
    db.commit()
    metrics['site_work_components'] = component_id
    metrics['site_clusters'] = db.execute('SELECT count(*) FROM clusters').fetchone()[0]


def locus_sort_key(key, rows, features):
    if key.startswith('REF::'):
        fid = key.split('::', 1)[1]
        rec = features[fid]
        return [0, str(rec.get('primary_chr') or ''), int(rec.get('primary_start1') or 0),
                int(rec.get('primary_end1') or 0), fid, key]
    if key.startswith('NRSITE::'):
        parts = key.split('::')
        chrom = parts[1] if len(parts) > 1 else ''
        try:
            ordinal = int(parts[2])
        except Exception:
            ordinal = 10**12
        starts = [int(r.get('nonref_site_start1') or 0) for r in rows if r.get('nonref_site_start1') not in ('', None)]
        ends = [int(r.get('nonref_site_end1') or 0) for r in rows if r.get('nonref_site_end1') not in ('', None)]
        return [1, chrom, min(starts) if starts else 0, max(ends) if ends else 0, ordinal, key]
    seed = key.split('::', 1)[1]
    try:
        return [2, 0, int(seed), key]
    except Exception:
        return [2, 1, seed, key]


def group_rows(db, key, args, metrics):
    return bounded_list((marshal.loads(row[0]) for row in db.execute(
        'SELECT payload FROM assignments WHERE anchor_key=? ORDER BY ordinal', (key,))),
        args.anchor_group_bytes, 'locus_group', metrics)


def aggregate_and_write(db, args, features, weights, metrics, metadata, progress=None):
    from pancgi_anchor_tables import describe_locus, member_values, locus_values, MEMBER_HEADER, LOCUS_HEADER
    db.execute('INSERT INTO loci(anchor_key,sk) SELECT DISTINCT anchor_key,\'\' FROM assignments')
    emit(progress, 'locus_order', 'start')
    for (key,) in db.execute('SELECT anchor_key FROM loci ORDER BY anchor_key'):
        rows = group_rows(db, key, args, metrics)
        db.execute('UPDATE loci SET sk=? WHERE anchor_key=?', (key_text(locus_sort_key(key, rows, metadata)), key))
        del rows
    db.execute('CREATE INDEX locus_order ON loci(sk COLLATE PYKEY)')
    db.commit()
    emit(progress, 'locus_order', 'complete')
    emit(progress, 'locus_science', 'start')
    for idx, (key,) in enumerate(db.execute('SELECT anchor_key FROM loci ORDER BY sk COLLATE PYKEY'), 1):
        rows = group_rows(db, key, args, metrics)
        fids = [str(r['fid']) for r in rows]
        insert_count = sum(int(r.get('insert_anchor_flag', 0) or 0) == 1 for r in rows)
        locus_id, info = describe_locus(idx, key, fids, rows, features, weights, args, {key: insert_count})
        if deep_bytes(info) > args.anchor_group_bytes:
            raise MemoryError('Locus metadata exceeds byte limit')
        db.execute('UPDATE loci SET locus_id=?,payload=? WHERE anchor_key=?', (int(locus_id), marshal.dumps(info), key))
        del rows, fids, info
        if idx % 1000 == 0:
            db.commit()
            emit(progress, 'locus_science', 'progress', idx)
    db.execute('CREATE UNIQUE INDEX locus_ids ON loci(locus_id)')
    db.commit()
    emit(progress, 'locus_science', 'complete')
    
    cache = ByteCache(args.anchor_group_bytes)
    handle, writer = prod.write_tsv_header(args.out_members, MEMBER_HEADER)
    prepared_before = features.loads if hasattr(features, 'loads') else 0
    emit(progress, 'member_output', 'start')
    member_count = 0
    try:
        for key, payload in db.execute('SELECT anchor_key,payload FROM assignments ORDER BY ordinal'):
            def read():
                lid, packed = db.execute('SELECT locus_id,payload FROM loci WHERE anchor_key=?', (key,)).fetchone()
                return str(lid), marshal.loads(packed)
            locus_id, info = cache.get_or_load(key, read)
            row = marshal.loads(payload)
            writer.writerow(member_values(row, locus_id, info, metadata[row['fid']]))
            member_count += 1
            if member_count % 10000 == 0:
                emit(progress, 'member_output', 'progress', member_count)
    finally:
        handle.close()
    metrics['member_output_records'] = member_count
    metrics['member_output_prepared_loads'] = (features.loads if hasattr(features, 'loads') else 0) - prepared_before
    if metrics['member_output_prepared_loads']:
        raise RuntimeError('Member output must not prepare complete features')
    emit(progress, 'member_output', 'complete', member_count, member_count)
    handle, writer = prod.write_tsv_header(args.out_locus, LOCUS_HEADER)
    try:
        for locus_id, payload in db.execute('SELECT locus_id,payload FROM loci ORDER BY locus_id'):
            writer.writerow(locus_values(str(locus_id), marshal.loads(payload)))
    finally:
        handle.close()
    metrics['locus_count'] = db.execute('SELECT count(*) FROM loci').fetchone()[0]
    metrics['locus_output_cache_peak_bytes'] = cache.peak


def execute(args, features):
    for name in ('anchor_group_bytes', 'anchor_weight_cache_bytes'):
        if getattr(args, name) <= 0:
            raise ValueError(name + ' must be positive')
    work = Path(args.out_locus).parent
    lookup_path, spool_path = work / 'lookup.sqlite', work / 'assignments.sqlite'
    metrics = {'strategy': 'disk_backed_metadata', 'phase_seconds': {}}
    lookup, spool = None, None
    progress = PhaseProgress(getattr(args, 'anchor_progress_path', work / 'progress.jsonl'), features)
    try:
        started = time.monotonic()
        build_lookup(lookup_path, args, features, metrics, progress)
        metrics['phase_seconds']['lookup'] = time.monotonic() - started
        lookup = Lookup(lookup_path)
        weights = KeyMap(lookup, 'weights', args.anchor_weight_cache_bytes)
        references = ReferenceMap(lookup, features)
        intervals = IntervalMap(lookup, args.anchor_group_bytes)
        metadata = MetadataMap(lookup, args.anchor_cache_bytes)
        prod._ANCHOR_FEAT_IDX = features
        prod._ANCHOR_SEED_LOCUS_BY_FID = KeyMap(lookup, 'seeds')
        prod._ANCHOR_PRIMARY_REFS = references


        prod._ANCHOR_REF_NODE_INDEX = {}
        prod._ANCHOR_REF_INTERVAL_INDEX = intervals
        prod._ANCHOR_SH_WEIGHTS = weights
        prod._ANCHOR_ARGS = args
        spool = create_spool(spool_path)
        tasks = (row[0] for row in lookup.rows('SELECT fid FROM tasks ORDER BY sk COLLATE PYKEY'))
        started = time.monotonic()
        count = 0
        emit(progress, 'assignment', 'start', 0, len(features))
        for count, row in enumerate(ordered_assignments(prod._anchor_assign_one_fid, tasks, args), 1):
            add_assignment(spool, count, row, args)
            if count % 10000 == 0:
                spool.commit()
                emit(progress, 'assignment', 'progress', count, len(features))
        spool.commit()
        if count != len(features):
            raise ValueError('Assignment count differs from full background')
        metrics['assignment_count'] = count
        metrics['phase_seconds']['assignment'] = time.monotonic() - started
        emit(progress, 'assignment', 'complete', count, len(features))
        started = time.monotonic()
        emit(progress, 'site_clustering', 'start')
        cluster_sites(spool, args, metrics)
        metrics['phase_seconds']['site_clustering'] = time.monotonic() - started
        emit(progress, 'site_clustering', 'complete')
        started = time.monotonic()
        aggregate_and_write(spool, args, features, weights, metrics, metadata, progress)
        metrics['phase_seconds']['locus_and_output'] = time.monotonic() - started
        metrics['parent_weight_cache_peak_bytes'] = weights.cache.peak
        metrics['parent_interval_cache_peak_bytes'] = intervals.cache.peak
        metrics['parent_metadata_cache_peak_bytes'] = metadata.cache.peak
        return metrics
    finally:
        progress.close()
        if spool is not None:
            spool.close()
        if lookup is not None:
            lookup.close()
        metrics['scratch_file_bytes'] = {p.name: p.stat().st_size for p in work.glob('*.sqlite')}
        
        with (work / 'disk_metrics.json').open('x') as handle:
            json.dump(metrics, handle, indent=2, sort_keys=True)
            handle.write('\n')
