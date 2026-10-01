import json
import sqlite3
from collections.abc import Mapping
from functools import lru_cache


class RowIndex(Mapping):
    def __init__(self, connection, table, count):
        self.connection, self.table, self.count = connection, table, count
        self._get = lru_cache(maxsize=512)(self._lookup)

    def _lookup(self, key):
        row = self.connection.execute(f'SELECT payload FROM {self.table} WHERE id=?', (key,)).fetchone()
        if row is None:
            raise KeyError(key)
        return json.loads(row[0])

    def __getitem__(self, key):
        return self._get(key)

    def __len__(self):
        return self.count

    def __iter__(self):
        for row in self.connection.execute(f'SELECT id FROM {self.table} ORDER BY ordinal'):
            yield row[0]

    def items(self):
        for key, value in self.connection.execute(f'SELECT id,payload FROM {self.table} ORDER BY ordinal'):
            yield key, json.loads(value)

    def values(self):
        for row in self.connection.execute(f'SELECT payload FROM {self.table} ORDER BY ordinal'):
            yield json.loads(row[0])


class RowSequence:
    def __init__(self, index, key):
        self.index, self.key = index, key

    def __iter__(self):
        return self.index.values()

    def __len__(self):
        return len(self.index)


class ValidationStore:
    def __init__(self, path):
        self.connection = sqlite3.connect(path)
        self.connection.execute('PRAGMA cache_size=-32768')
        self.tables = 0

    def add(self, rows, key, check):
        table = f't{self.tables}'
        self.tables += 1
        self.connection.execute(f'CREATE TABLE {table} (ordinal INTEGER PRIMARY KEY, id TEXT UNIQUE, payload TEXT)')
        count = 0
        for count, row in enumerate(rows, 1):
            check(row)
            identifier = row[key] if key is not None else str(count)
            try:
                self.connection.execute(f'INSERT INTO {table} VALUES (?,?,?)', (count, identifier, json.dumps(row)))
            except sqlite3.IntegrityError as exc:
                raise ValueError(f'Duplicate {key}') from exc
            if count % 10000 == 0:
                self.connection.commit()
        self.connection.commit()
        return RowSequence(RowIndex(self.connection, table, count), key)

    def close(self):
        self.connection.close()
