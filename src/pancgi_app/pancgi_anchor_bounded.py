import copy
import json
import multiprocessing
import os
import shutil
import sqlite3
import sys
import tempfile
import threading
import time
import zlib
from collections import OrderedDict, deque
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import pancgi_features as shared
from pancgi_contract import digest, file_state

METADATA_FIELDS = ('fid', 'kind', 'primary_chr', 'primary_start1',
                   'primary_end1', 'nonref_only')


def deep_bytes(value):
    seen = set()
    stack = [value]
    size = 0
    while stack:
        item = stack.pop()
        if id(item) in seen:
            continue
        seen.add(id(item))
        size += sys.getsizeof(item)
        if isinstance(item, dict):
            stack.extend(item.keys())
            stack.extend(item.values())
        elif isinstance(item, (list, tuple, set, frozenset)):
            stack.extend(item)
    return size


class FeatureMap(Mapping):

    def __init__(self, path, anchor_k, max_mid_anchors, cache_bytes, max_record_bytes):
        if cache_bytes <= 0 or max_record_bytes <= 0:
            raise ValueError('Feature resource limits must be positive')
        self.path = str(Path(path).resolve())
        self.state = file_state(self.path)
        self.anchor_k = anchor_k
        self.max_mid_anchors = max_mid_anchors
        self.cache_limit = cache_bytes
        self.max_record_bytes = max_record_bytes
        self.pid = os.getpid()
        self.lock = threading.RLock()
        self.conn = self._connect()
        try:
            self.info = shared.metadata(self.conn, anchor_k, max_mid_anchors)
        except BaseException:
            self.conn.close()
            raise
        self.data_version = self.conn.execute('PRAGMA data_version').fetchone()[0]
        self.cache = OrderedDict()
        self.cache_bytes = self.peak_cache_bytes = self.peak_cache_records = 0
        self.loads = self.hits = self.evictions = self.oversize_uncached = 0
        self.raw_loads = self.decoded_raw_bytes = self.metadata_reads = 0
        actual = self.conn.execute('SELECT count(*) FROM polish_feature_store').fetchone()[0]
        if actual != self.info['record_count']:
            self.close()
            raise ValueError('Shared feature count mismatch')

    def _connect(self):
        conn = shared.connect(self.path, cache_kib=4096, check_same_thread=False)
        conn.execute('PRAGMA mmap_size=0')
        conn.execute('PRAGMA temp_store=FILE')
        return conn

    def _local(self):
        
        if self.pid != os.getpid():
            self.conn.check_files()
            self.conn.close()
            self.pid = os.getpid()
            self.lock = threading.RLock()
            self.conn = self._connect()
            shared.metadata(self.conn, self.anchor_k, self.max_mid_anchors)
            self.data_version = self.conn.execute('PRAGMA data_version').fetchone()[0]
            self.cache = OrderedDict()
            self.cache_bytes = self.peak_cache_bytes = self.peak_cache_records = 0
            self.loads = self.hits = self.evictions = self.oversize_uncached = 0
            self.raw_loads = self.decoded_raw_bytes = self.metadata_reads = 0

    def __len__(self):
        return int(self.info['record_count'])

    def __iter__(self):
        self._local()
        
        cursor = self.conn.execute('SELECT fid FROM polish_feature_store ORDER BY ordinal')
        try:
            for row in cursor:
                yield row[0]
        finally:
            cursor.close()

    def __getitem__(self, fid):
        import cpgi_nr_prod as prod
        self._local()
        with self.lock:
            if fid in self.cache:
                self.hits += 1
                self.cache.move_to_end(fid)
                return self.cache[fid][0]
            row = self.conn.execute('SELECT raw FROM polish_feature_store WHERE fid=? AND length(raw)<=?',
                                    (fid, 2 * self.max_record_bytes + 1024)).fetchone()
            if row is None:
                if self.conn.execute('SELECT 1 FROM polish_feature_store WHERE fid=?', (fid,)).fetchone():
                    raise MemoryError('Feature exceeds compressed record byte limit')
                raise KeyError(fid)
            record = self._decode_raw(fid, row[0])
            feature = prod.prep_feature(record, self.anchor_k, self.max_mid_anchors)
            size = deep_bytes(feature) + sys.getsizeof(fid) + 128
            self.loads += 1
            while self.cache and self.cache_bytes + size > self.cache_limit:
                _, (_, removed) = self.cache.popitem(last=False)
                self.cache_bytes -= removed
                self.evictions += 1
            if size <= self.cache_limit:
                self.cache[fid] = (feature, size)
                self.cache_bytes += size
                self.peak_cache_bytes = max(self.peak_cache_bytes, self.cache_bytes)
                self.peak_cache_records = max(self.peak_cache_records, len(self.cache))
            else:
                self.oversize_uncached += 1
            return feature

    def _decode_raw(self, fid, packed):
        if packed is None or len(packed) > 2 * self.max_record_bytes + 1024:
            raise MemoryError('Feature exceeds compressed record byte limit')
        decoder = zlib.decompressobj()
        raw = decoder.decompress(packed, self.max_record_bytes + 1)
        if len(raw) > self.max_record_bytes or not decoder.eof:
            raise MemoryError('Feature exceeds raw record byte limit; no truncation allowed')
        if decoder.unused_data:
            raise ValueError('Trailing data in compressed feature')
        record = json.loads(raw)
        if str(record['fid']) != fid:
            raise ValueError('Feature index/raw identity mismatch')
        self.raw_loads += 1
        self.decoded_raw_bytes += len(raw)
        return record

    def _metadata_from_payload(self, fid, packed):
        record = self._decode_raw(fid, packed)
        self.metadata_reads += 1
        return {key: record[key] for key in METADATA_FIELDS if key in record}

    def iter_metadata(self):
        self._local()
        cursor = self.conn.execute(
            'SELECT fid, CASE WHEN length(raw)<=? THEN raw ELSE NULL END '
            'FROM polish_feature_store ORDER BY ordinal',
            (2 * self.max_record_bytes + 1024,))
        count = 0
        try:
            while True:
                with self.lock:
                    row = cursor.fetchone()
                    if row is None:
                        break
                    fid, record = row[0], self._metadata_from_payload(*row)
                count += 1
                yield fid, record
            if count != len(self):
                raise ValueError('Metadata scan count differs from full background')
            self.check()
        finally:
            cursor.close()

    def check(self):
        self._local()
        if (file_state(self.path) != self.state or
                self.conn.execute('PRAGMA data_version').fetchone()[0] != self.data_version):
            raise RuntimeError('Shared feature index changed during anchor')
        shared.metadata(self.conn, self.anchor_k, self.max_mid_anchors)

    def statistics(self):
        return {k: getattr(self, k) for k in ('cache_limit', 'cache_bytes',
            'peak_cache_bytes', 'peak_cache_records', 'loads', 'hits', 'evictions',
            'oversize_uncached', 'raw_loads', 'decoded_raw_bytes', 'metadata_reads', 'pid')}

    def close(self):
        self.conn.close()


def ordered_assignments(function, tasks, args):
    workers = int(args.threads)
    limit = int(args.anchor_max_pending)
    if workers < 1 or limit < 1:
        raise ValueError('Workers and pending-task limit must be positive')
    if workers == 1:
        yield from map(function, tasks)
        return
    if args.parallel_backend == 'process':
        if args.mp_start_method != 'fork' or args.maxtasksperchild != 0:
            raise ValueError('Anchor process backend requires fork and maxtasksperchild=0')
        executor = ProcessPoolExecutor(max_workers=workers,
            mp_context=multiprocessing.get_context('fork'))
    else:
        executor = ThreadPoolExecutor(max_workers=workers)
    pending = deque()
    source = iter(tasks)
    try:
        for _ in range(limit):
            try:
                task = next(source)
            except StopIteration:
                break
            pending.append(executor.submit(function, task))
        while pending:
            yield pending.popleft().result()
            try:
                task = next(source)
            except StopIteration:
                continue
            pending.append(executor.submit(function, task))
    finally:
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)


def run(args, implementation):
    if not str(args.features).endswith('.sqlite'):
        raise ValueError('Anchor requires a validated shared SQLite feature index')
    if args.threads < 1 or args.anchor_max_pending < 1:
        raise ValueError('Workers and pending-task limit must be positive')
    targets = [Path(args.out_members).resolve(), Path(args.out_locus).resolve()]
    report_path = Path(args.anchor_resource_report or str(targets[1]) + '.anchor.json').resolve()
    paths = [Path(args.features).resolve(), Path(args.in_locus).resolve(), Path(args.in_members).resolve()]
    progress_path = Path(str(report_path) + '.progress.jsonl')
    if len(set(targets + [report_path, progress_path])) != 4 or set(targets + [report_path, progress_path]) & set(paths):
        raise ValueError('Anchor output paths must be distinct from one another and inputs')
    for path in targets + [report_path, progress_path]:
        if path.exists():
            raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
    states = {str(path): file_state(path) for path in paths}
    started = time.time()
    work = Path(tempfile.mkdtemp(prefix='anchor_', dir=report_path.parent))
    previous_sqlite_tmp = os.environ.get('SQLITE_TMPDIR')
    os.environ['SQLITE_TMPDIR'] = str(work)
    staged_args = copy.copy(args)
    staged_args.anchor_progress_path = str(progress_path)
    staged_args.out_members = str(work / targets[0].name)
    staged_args.out_locus = str(work / targets[1].name)
    if targets[0].name == targets[1].name:
        staged_args.out_members = str(work / ('members_' + targets[0].name))
    reader = None
    published = []
    report = dict(status='running', schema='pancgi_anchor_execution_v1',
        scope='anchor',
        source_states=states, arguments=vars(args).copy(), work_dir=str(work),
        progress_path=str(progress_path))
    report['arguments'].pop('func', None)
    try:
        reader = FeatureMap(args.features, args.anchor_k, args.max_mid_anchors,
                            args.anchor_cache_bytes, args.anchor_max_record_bytes)
        if Path(reader.info['source']).resolve() in set(targets + [report_path]):
            raise ValueError('Anchor output aliases the original source')
        report['disk_pipeline'] = implementation(staged_args, reader)
        reader.check()
        for path, state in states.items():
            if file_state(path) != state:
                raise RuntimeError('Anchor input changed: ' + path)
        report['record_count'] = len(reader)
        report['outputs_sha256'] = {}
        
        for staged, target in zip([staged_args.out_members, staged_args.out_locus], targets):
            os.link(staged, target)
            published.append((Path(staged), target))
            report['outputs_sha256'][str(target)] = digest(target)
        report['status'] = 'completed'
    except BaseException as exc:
        report.update(status='failed', error=type(exc).__name__ + ': ' + str(exc))
        for staged, target in published:
            if target.exists() and target.samefile(staged):
                target.unlink()
        raise
    finally:
        if reader is not None:
            report['parent_feature_cache'] = reader.statistics()
            reader.close()
        import cpgi_nr_prod as prod
        for name in ('_ANCHOR_FEAT_IDX', '_ANCHOR_SEED_LOCUS_BY_FID', '_ANCHOR_PRIMARY_REFS',
                     '_ANCHOR_REF_NODE_INDEX', '_ANCHOR_REF_INTERVAL_INDEX',
                     '_ANCHOR_SH_WEIGHTS', '_ANCHOR_ARGS'):
            setattr(prod, name, None)
        report['elapsed_seconds'] = time.time() - started
        with report_path.open('x') as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write('\n')
        if report['status'] != 'failed':
            shutil.rmtree(work)
        if previous_sqlite_tmp is None:
            os.environ.pop('SQLITE_TMPDIR', None)
        else:
            os.environ['SQLITE_TMPDIR'] = previous_sqlite_tmp
