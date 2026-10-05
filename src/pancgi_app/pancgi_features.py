import argparse
import gzip
import hashlib
import io
import json
import sqlite3
import zlib
from collections import Counter
from contextlib import closing
from functools import lru_cache
from pathlib import Path

from pancgi_contract import digest, file_state, open_text


def _index_state(path):
    result = [tuple(file_state(path))]
    for suffix in ('-wal', '-journal'):
        try:
            result.append(tuple(file_state(str(path) + suffix)))
        except FileNotFoundError:
            result.append(None)
    return tuple(result)


class _ReadConnection(sqlite3.Connection):
    def check_files(self):
        if _index_state(self.index_path) != self.index_state:
            raise RuntimeError('Shared feature index changed during reading')
        if self.source_path is not None and file_state(self.source_path) != self.source_state:
            raise ValueError('Source features changed during reading')

    def check(self):
        self.check_files()
        if self.execute('PRAGMA data_version').fetchone()[0] != self.index_version:
            raise RuntimeError('Shared feature index changed during reading')


class _HashingReader(io.RawIOBase):
    def __init__(self, source, hasher):
        super().__init__()
        self.source, self.hasher = source, hasher

    def readable(self):
        return True

    def readinto(self, buffer):
        size = self.source.readinto(buffer)
        if size:
            self.hasher.update(memoryview(buffer)[:size])
        return size


def connect(path, *, cache_kib=32768, check_same_thread=True):
    path = Path(path).resolve()
    before = _index_state(path)
    db = sqlite3.connect(path.as_uri() + '?mode=ro', uri=True,
                         factory=_ReadConnection, check_same_thread=check_same_thread)
    db.index_path = str(path)
    db.source_path = db.source_state = None
    try:
        db.execute('PRAGMA query_only=ON')
        db.execute(f'PRAGMA cache_size=-{int(cache_kib)}')
        db.index_version = db.execute('PRAGMA data_version').fetchone()[0]
        db.index_state = _index_state(path)
        if db.index_state[0] != before[0]:
            raise RuntimeError('Shared feature index changed while opening')
        db.check()
    except BaseException:
        db.close()
        raise
    return db


def _compare_source_records(db, info):
    cursor = db.execute('SELECT ordinal, fid, raw FROM polish_feature_store ORDER BY ordinal')
    count = 0
    try:
        with open_text(info['source']) as handle:
            for line in handle:
                if not line.strip():
                    continue
                count += 1
                row = cursor.fetchone()
                if row is None or row[0] != count:
                    raise ValueError('Source features changed: shared-store order or count mismatch')
                expected = line.encode('utf-8')
                decoder = zlib.decompressobj()
                actual = decoder.decompress(row[2], len(expected) + 1)
                if (actual != expected or not decoder.eof or decoder.unused_data or
                        str(json.loads(line)['fid']) != row[1]):
                    raise ValueError('Source features changed: shared-store content mismatch')
        if count != info['record_count'] or cursor.fetchone() is not None:
            raise ValueError('Source features changed: shared-store count mismatch')
    finally:
        cursor.close()


@lru_cache(maxsize=8)
def _verify_source(index_path, index_state, payload, source_state):
    info = json.loads(payload)
    if tuple(file_state(info['source'])) != source_state:
        raise ValueError('Source features changed before identity verification')
    if _index_state(index_path) != index_state:
        raise RuntimeError('Shared feature index changed before identity verification')
    if list(source_state) != info['source_state']:
        if info['schema_version'] == 2:
            if digest(info['source']) != info['source_sha256']:
                raise ValueError('Source features changed: SHA256 mismatch')
        else:
            with closing(connect(index_path)) as verification:
                _compare_source_records(verification, info)
                verification.check()
    if tuple(file_state(info['source'])) != source_state:
        raise ValueError('Source features changed during identity verification')
    if _index_state(index_path) != index_state:
        raise RuntimeError('Shared feature index changed during identity verification')


def metadata(db, anchor_k=None, max_mid_anchors=None):
    if isinstance(db, _ReadConnection):
        db.check()
    index_path = db.execute('PRAGMA database_list').fetchone()[2]
    before = _index_state(index_path)
    version = db.execute('PRAGMA data_version').fetchone()[0]
    row = db.execute('SELECT payload FROM shared_metadata').fetchone()
    if row is None:
        raise ValueError('Incomplete shared feature store')
    result = json.loads(row[0])
    if result['schema_version'] not in (1, 2) or result['status'] != 'complete':
        raise ValueError('Unsupported shared feature store')
    if (not isinstance(result['source_state'], list) or len(result['source_state']) != 5 or
            any(type(value) is not int for value in result['source_state'])):
        raise ValueError('Invalid shared feature source state')
    if result['schema_version'] == 2:
        fingerprint = result.get('source_sha256')
        if (not isinstance(fingerprint, str) or len(fingerprint) != 64 or
                any(c not in '0123456789abcdef' for c in fingerprint)):
            raise ValueError('Missing or invalid shared feature source SHA256')
    for key, value in [('anchor_k', anchor_k), ('max_mid_anchors', max_mid_anchors)]:
        if value is not None and result[key] != value:
            raise ValueError(f'Shared feature parameter mismatch: {key}')
    observed = file_state(result['source'])
    _verify_source(index_path, tuple(before), row[0], tuple(observed))
    if (_index_state(index_path) != before or
            db.execute('PRAGMA data_version').fetchone()[0] != version):
        raise RuntimeError('Shared feature index changed during identity verification')
    if file_state(result['source']) != observed:
        raise ValueError('Source features changed during identity verification')
    if isinstance(db, _ReadConnection):
        db.source_path, db.source_state = result['source'], observed
        db.check()
    return result


def frequencies(path, anchor_k, max_mid_anchors):
    with closing(connect(path)) as db:
        info = metadata(db, anchor_k, max_mid_anchors)
        row = db.execute('SELECT payload FROM shared_frequencies').fetchone()
        if row is None:
            raise ValueError('Missing global feature frequencies')
        df = json.loads(zlib.decompress(row[0]))
        db.check()
        return info['record_count'], df


def iter_raw(path):
    with closing(connect(path)) as db:
        info = metadata(db)
        count = 0
        try:
            for row in db.execute('SELECT raw FROM polish_feature_store ORDER BY ordinal'):
                count += 1
                yield json.loads(zlib.decompress(row[0]))
            if count != info['record_count']:
                raise ValueError('Shared feature count mismatch')
        finally:
            db.check()


def selected_raw(path, fids):
    wanted = list(dict.fromkeys(fids))
    result = {}
    order = {}
    with closing(connect(path)) as db:
        metadata(db)
        for first in range(0, len(wanted), 900):
            part = wanted[first:first+900]
            sql = 'SELECT fid,raw,ordinal FROM polish_feature_store WHERE fid IN (' + ','.join('?' for _ in part) + ')'
            for fid, raw, ordinal in db.execute(sql, part):
                result[fid] = json.loads(zlib.decompress(raw))
                order[fid] = ordinal
        db.check()
    if set(result) != set(wanted):
        raise ValueError('Missing selected features in shared store')
    return {fid: result[fid] for fid in sorted(result, key=order.__getitem__)}


def build(source, output, anchor_k, max_mid_anchors):
    from cpgi_nr_prod import prep_feature, _polish_store_pack_record
    source, output = Path(source).resolve(), Path(output).resolve()
    partial = output.with_name(output.name + '.partial')
    if output.exists() or partial.exists():
        raise FileExistsError(output)
    state = file_state(source)
    output.parent.mkdir(parents=True, exist_ok=True)
    n, df = 0, Counter()
    with closing(sqlite3.connect(partial)) as db:
        db.execute('PRAGMA cache_size=-32768')
        db.execute('CREATE TABLE polish_feature_store (ordinal INTEGER PRIMARY KEY, fid TEXT NOT NULL UNIQUE, blob BLOB NOT NULL, raw BLOB NOT NULL)')
        db.execute('CREATE TABLE shared_metadata (payload TEXT NOT NULL)')
        db.execute('CREATE TABLE shared_frequencies (payload BLOB NOT NULL)')
        hasher = hashlib.sha256()
        with source.open('rb') as raw_source, io.BufferedReader(_HashingReader(raw_source, hasher)) as hashed:
            stream = gzip.GzipFile(fileobj=hashed, mode='rb') if str(source).endswith('.gz') else hashed
            with io.TextIOWrapper(stream, encoding='utf-8', newline='') as handle:
                for raw in handle:
                    if not raw.strip():
                        continue
                    rec = json.loads(raw)
                    prepared = prep_feature(rec, anchor_k, max_mid_anchors)
                    df.update(set(prepared['_shingles']))
                    n += 1
                    db.execute('INSERT INTO polish_feature_store VALUES (?,?,?,?)',
                        (n, str(rec['fid']), _polish_store_pack_record(prepared), zlib.compress(raw.encode(), 1)))
                    if n % 10000 == 0:
                        db.commit()
        if file_state(source) != state:
            raise RuntimeError('Source features changed during shared-store construction')
        report = dict(schema_version=2, status='complete', source=str(source), source_state=state,
            source_sha256=hasher.hexdigest(),
            anchor_k=anchor_k, max_mid_anchors=max_mid_anchors, record_count=n,
            universe='all accepted features in original input order')
        db.execute('INSERT INTO shared_metadata VALUES (?)', (json.dumps(report),))
        db.execute('INSERT INTO shared_frequencies VALUES (?)', (zlib.compress(json.dumps(dict(df)).encode(), 1),))
        db.commit()
    partial.replace(output)
    return report


def main():
    p = argparse.ArgumentParser(description='Build one complete feature index and global shingle-frequency artifact')
    p.add_argument('--features', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--anchor-k', type=int, required=True)
    p.add_argument('--max-mid-anchors', type=int, required=True)
    a = p.parse_args()
    print(json.dumps(build(a.features, a.out, a.anchor_k, a.max_mid_anchors), indent=2))


if __name__ == '__main__':
    main()
