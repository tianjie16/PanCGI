import argparse
import json
import sqlite3
import zlib
from collections import Counter
from contextlib import closing
from pathlib import Path

from pancgi_contract import file_state, open_text


def connect(path):
    db = sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True)
    db.execute('PRAGMA cache_size=-32768')
    return db


def metadata(db, anchor_k=None, max_mid_anchors=None):
    row = db.execute('SELECT payload FROM shared_metadata').fetchone()
    if row is None:
        raise ValueError('Incomplete shared feature store')
    result = json.loads(row[0])
    if result['schema_version'] != 1 or result['status'] != 'complete':
        raise ValueError('Unsupported shared feature store')
    if file_state(result['source']) != result['source_state']:
        raise ValueError('Source features changed after shared-store construction')
    for key, value in [('anchor_k', anchor_k), ('max_mid_anchors', max_mid_anchors)]:
        if value is not None and result[key] != value:
            raise ValueError(f'Shared feature parameter mismatch: {key}')
    return result


def frequencies(path, anchor_k, max_mid_anchors):
    with closing(connect(path)) as db:
        info = metadata(db, anchor_k, max_mid_anchors)
        row = db.execute('SELECT payload FROM shared_frequencies').fetchone()
        if row is None:
            raise ValueError('Missing global feature frequencies')
        df = json.loads(zlib.decompress(row[0]))
        return info['record_count'], df


def iter_raw(path):
    with closing(connect(path)) as db:
        info = metadata(db)
        count = 0
        for row in db.execute('SELECT raw FROM polish_feature_store ORDER BY ordinal'):
            count += 1
            yield json.loads(zlib.decompress(row[0]))
        if count != info['record_count']:
            raise ValueError('Shared feature count mismatch')


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
        with open_text(source) as handle:
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
        report = dict(schema_version=1, status='complete', source=str(source), source_state=state,
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
