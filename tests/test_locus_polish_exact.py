from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def write_jsonl_gz(path: Path, rows) -> None:
    with gzip.open(path, 'wt') as fh:
        for row in rows:
            fh.write(json.dumps(row) + '\n')


def write_tsv(path: Path, header, rows) -> None:
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def make_feat(fid, label, kind, steps, *, chr1='', start1='', end1=''):
    seq_len = len(steps) * 10
    return {
        'fid': fid,
        'label': label,
        'kind': kind,
        'steps': [[int(x), 0] for x in steps],
        'nodeints': [[int(x), 0, 10] for x in steps],
        'source_seq': 'A' * seq_len,
        'source_seq_len': seq_len,
        'primary_chr': chr1,
        'primary_start1': start1,
        'primary_end1': end1,
        'nonref_only': 0 if chr1 else 1,
    }


def read_text(path: Path) -> str:
    if str(path).endswith('.gz'):
        with gzip.open(path, 'rt') as fh:
            return fh.read()
    return path.read_text()


def test_locus_polishing_is_deterministic() -> None:
    here = Path(__file__).resolve().parents[1]
    script = here / 'cpgi_nr_prod.py'
    with tempfile.TemporaryDirectory(prefix='pancgi_polishing_') as td:
        td = Path(td)
        features = td / 'features.jsonl.gz'
        in_locus = td / 'locus_seed.tsv.gz'
        in_members = td / 'locus_seed_members.tsv.gz'
        out_first_locus = td / 'first.locus.tsv.gz'
        out_first_members = td / 'first.members.tsv.gz'
        out_second_locus = td / 'second.locus.tsv.gz'
        out_second_members = td / 'second.members.tsv.gz'
        feature_store = td / 'polish.sqlite'

        rows = [
            make_feat('PRIMARY#0#chr1:100-129', 'PRIMARY_hap0', 'reference_primary', [1, 2, 3], chr1='chr1', start1=100, end1=129),
            make_feat('ASM1#1#ctg1:0-29', 'ASM1_hap1', 'assembly', [1, 2, 3]),
            make_feat('ASM2#1#ctg1:0-29', 'ASM2_hap1', 'assembly', [1, 2, 4]),
            make_feat('COMPARISON#0#chr2:200-229', 'COMPARISON_hap0', 'reference_comparison', [10, 11, 12]),
            make_feat('ASM3#1#ctg2:0-29', 'ASM3_hap1', 'assembly', [10, 11, 12]),
            make_feat('ASM4#1#ctgX:0-29', 'ASM4_hap1', 'assembly', [1, 2, 3]),
        ]
        write_jsonl_gz(features, rows)
        write_tsv(in_locus, ['locus_id', 'lead_fid'], [
            ['1', 'PRIMARY#0#chr1:100-129'],
            ['2', 'COMPARISON#0#chr2:200-229'],
        ])
        write_tsv(in_members, ['locus_id', 'fid'], [
            ['1', 'PRIMARY#0#chr1:100-129'],
            ['1', 'ASM1#1#ctg1:0-29'],
            ['1', 'ASM2#1#ctg1:0-29'],
            ['2', 'COMPARISON#0#chr2:200-229'],
            ['2', 'ASM3#1#ctg2:0-29'],
        ])

        common = [
            sys.executable, str(script), 'polish-loci',
            '--features', str(features),
            '--in-locus', str(in_locus),
            '--in-members', str(in_members),
            '--anchor-k', '3',
            '--max-mid-anchors', '8',
            '--threads', '1',
            '--chunksize', '1',
            '--feature-store-db', str(feature_store),
        ]
        subprocess.run(common + ['--out-locus', str(out_first_locus), '--out-members', str(out_first_members)], check=True)
        subprocess.run(common + ['--out-locus', str(out_second_locus), '--out-members', str(out_second_members)], check=True)

        assert feature_store.exists()
        assert read_text(out_first_locus) == read_text(out_second_locus)
        assert read_text(out_first_members) == read_text(out_second_members)
        assert 'ga_cov' not in read_text(out_first_members).splitlines()[0]


if __name__ == '__main__':
    test_locus_polishing_is_deterministic()
    print('OK: locus polishing is deterministic on synthetic input')
