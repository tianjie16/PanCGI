import json
import math
import sqlite3
import tempfile
from pathlib import Path

from .directions import choose, iter_pair_candidates, placements, placements_sorted


class CandidateStore:
    def __init__(self, path, records, buffer_n=10000):
        self.path = Path(path)
        self.conn = sqlite3.connect(self.path)
        self.conn.execute('PRAGMA journal_mode=OFF')
        self.conn.execute('PRAGMA synchronous=OFF')
        self.conn.execute('PRAGMA temp_store=FILE')
        self.conn.execute('PRAGMA cache_size=-32768')
        self.conn.execute('CREATE TABLE candidates (contig TEXT, direction INT, lo REAL, hi REAL, model TEXT, window INT, score REAL, seq INT, payload TEXT, PRIMARY KEY(contig,direction,lo,hi,model,window)) WITHOUT ROWID')
        self.raw_n = 0
        self.groups_cache = {}
        pending = {}
        query = ('INSERT INTO candidates VALUES (?,?,?,?,?,?,?,?,?) '
                 'ON CONFLICT(contig,direction,lo,hi,model,window) DO UPDATE SET '
                 'score=excluded.score,seq=excluded.seq,payload=excluded.payload '
                 'WHERE excluded.score>candidates.score')

        def flush():
            self.conn.executemany(query, [(*key, c['score'], seq, json.dumps(c, separators=(',', ':')))
                                          for key, (c, seq) in pending.items()])
            self.conn.commit()
            pending.clear()

        for seq, c in enumerate(records):
            if not all(math.isfinite(c[k]) for k in ('score', 'lo', 'hi')):
                raise RuntimeError('nonfinite candidate evidence')
            self.raw_n += 1
            key = tuple(c[k] for k in ('contig', 'direction', 'lo', 'hi', 'model', 'window'))
            previous = pending.get(key)
            if previous is None or c['score'] > previous[0]['score']:
                pending[key] = (c, seq)
            if len(pending) >= buffer_n:
                flush()
        if pending:
            flush()
        self.unique_n = self.conn.execute('SELECT COUNT(*) FROM candidates').fetchone()[0]

    def ranked(self, bounds=None):
        query = 'SELECT payload FROM candidates'
        args = []
        if bounds is not None:
            contig, direction, lo, hi, tolerance = bounds
            query += ' WHERE contig=? AND direction=? AND lo BETWEEN ? AND ? AND hi BETWEEN ? AND ?'
            args = [contig, direction, lo - tolerance, lo + tolerance, hi - tolerance, hi + tolerance]
        query += ' ORDER BY score DESC,contig,direction,lo,hi,model,seq'
        for (payload,) in self.conn.execute(query, args):
            yield json.loads(payload)

    def groups(self, tolerance, bounds=None):
        key = (tolerance, bounds)
        if key not in self.groups_cache:
            self.groups_cache[key] = placements_sorted(self.ranked(bounds), tolerance)
        return self.groups_cache[key]

    def close(self):
        self.conn.close()


class CandidateCache:
    def __init__(self, scratch=None, buffer_n=10000):
        self.temporary = tempfile.TemporaryDirectory(prefix='pancgi_candidates.', dir=scratch)
        self.stores = {}
        self.buffer_n = buffer_n

    def get(self, record, window, params):
        key = (window,) + params
        if key not in self.stores:
            stream = (p for m in record['candidate_evidence'] for p in iter_pair_candidates(m, window, *params))
            self.stores[key] = CandidateStore(Path(self.temporary.name) / f'{len(self.stores)}.sqlite', stream, self.buffer_n)
        return self.stores[key]

    def diagnostics(self):
        return {'stores': len(self.stores), 'candidates_enumerated': sum(s.raw_n for s in self.stores.values()),
                'distinct_mapping_template_window_keys': sum(s.unique_n for s in self.stores.values()),
                'disk_bytes': sum(s.path.stat().st_size for s in self.stores.values())}

    def __enter__(self):
        return self

    def __exit__(self, *unused):
        for store in self.stores.values():
            store.close()
        self.temporary.cleanup()


