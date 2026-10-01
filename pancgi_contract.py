from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from decimal import Decimal
from pathlib import Path

import pancgi_mapping as mapping

GENOMES = ['hal_genome', 'role', 'cpgi_bed', 'sv_tsv']
PATHS = ['gfa_record_type', 'gfa_sample', 'gfa_haplotype', 'gfa_sequence', 'gfa_path_name', 'hal_genome', 'hal_sequence']
HAL_ONLY = ['hal_genome', 'hal_sequence']
HAL_ONLY_REPORT = HAL_ONLY + ['length_bp', 'reason', 'cgi_rows', 'sv_assembly_rows', 'sv_reference_rows']
SV = ['ref_contig', 'ref_pos1', 'ref_end1', 'sv_id', 'sv_type', 'asm_contig', 'asm_start0', 'asm_end0', 'sv_length']


def open_text(path, mode='rt'):
    return gzip.open(path, mode, encoding='utf-8', newline='') if str(path).endswith('.gz') else open(path, mode, encoding='utf-8', newline='')


def read_exact(path, columns, allow_empty=False):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, None)
        if header != columns:
            raise ValueError(f'{path}: expected exactly {columns}, observed {header}')
        result = []
        for line, fields in enumerate(reader, 2):
            if len(fields) != len(columns):
                raise ValueError(f'{path}:{line}: expected {len(columns)} fields')
            if any(value != value.strip() or '\x00' in value for value in fields):
                raise ValueError(f'{path}:{line}: whitespace or NUL in a field')
            result.append(dict(zip(columns, fields)))
    if not result and not allow_empty:
        raise ValueError(f'{path}: no records')
    return result


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def file_state(path):
    s = Path(path).stat()
    return [s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns]


def verify_lock(path):
    lock = json.loads(Path(path).read_text())
    if lock['status'] != 'validated':
        raise ValueError('Unvalidated input lock')
    for name, state in lock['file_states'].items():
        if file_state(name) != state:
            raise RuntimeError(f'Input changed since validation: {name}')
    root = Path(__file__).resolve().parent
    actual = {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob('*')) if p.is_file() and p.suffix in ('.py', '.sh')}
    if actual != lock['code']:
        raise RuntimeError('Pipeline source changed during this run')
    for name, sha in lock['internal_files'].items():
        if digest(Path(path).parent / name) != sha:
            raise RuntimeError('Validated identity mapping changed')
    return lock


def storage_id(prefix, *names):
    return prefix + hashlib.sha256(json.dumps(names, ensure_ascii=True, separators=(',', ':')).encode()).hexdigest()


def integer(value, context):
    if not value or not value.isascii() or not value.isdecimal():
        raise ValueError(f'{context}: expected a non-negative decimal integer')
    return int(value)


def parse_sv_length(value, sv_type, context):
    digits = value[1:] if value.startswith(('-', '+')) else value
    if not digits or not digits.isascii() or not digits.isdecimal():
        raise ValueError(f'{context}: sv_length must be the supplied signed integer SVLEN')
    length = int(value)
    if sv_type not in ('INS', 'DEL') or (sv_type == 'INS' and length <= 0) or (sv_type == 'DEL' and length >= 0):
        raise ValueError(f'{context}: sv_length must be positive for INS and negative for DEL')
    return length


def read_bed(path, lengths):
    seen = set()
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            f = raw.rstrip('\r\n').split('\t')
            if len(f) != 10:
                raise ValueError(f'{path}:{number}: expected headerless tab-delimited BED10')
            if f[0] not in lengths:
                raise ValueError(f'{path}:{number}: unknown HAL sequence {f[0]!r}')
            start, end = integer(f[1], path), integer(f[2], path)
            length, cpg, gc = (integer(f[i], path) for i in (4, 5, 6))
            if not 0 <= start < end <= lengths[f[0]] or length != end - start:
                raise ValueError(f'{path}:{number}: inconsistent length or out-of-bounds CGI')
            if not 0 <= 2 * cpg <= gc <= length:
                raise ValueError(f'{path}:{number}: inconsistent CpG/GC counts')
            metrics = [float(f[i]) for i in (7, 8, 9)]
            if not all(math.isfinite(v) and v >= 0 for v in metrics) or max(metrics[:2]) > 100:
                raise ValueError(f'{path}:{number}: invalid CGI metrics')
            expected = [200 * cpg / length, 100 * gc / length]
            for column, actual, target in zip((7, 8), metrics[:2], expected):
                unit = Decimal(f[column]).as_tuple().exponent
                tolerance = 0.5 * 10.0 ** min(0, unit) + 1e-9
                if abs(actual - target) > tolerance:
                    raise ValueError(f'{path}:{number}: percentage column {column + 1} contradicts counts')
            key = f[0], start, end
            if key in seen:
                raise ValueError(f'{path}:{number}: duplicate CGI interval')
            seen.add(key)
            yield f


def validate_sv(path, lengths, reference_lengths):
    records = read_exact(path, SV, allow_empty=True)
    seen = set()
    for row in records:
        if not row['sv_id'] or row['sv_id'] in seen:
            raise ValueError(f'{path}: empty or duplicated SV ID {row["sv_id"]!r}')
        seen.add(row['sv_id'])
        if row['sv_type'] not in ('INS', 'DEL'):
            raise ValueError(f'{path}: unsupported SV type {row["sv_type"]!r}; explicit INS/DEL events required')
        parse_sv_length(row['sv_length'], row['sv_type'], f'{path}:{row["sv_id"]}')
        if row['asm_contig'] not in lengths or row['ref_contig'] not in reference_lengths:
            raise ValueError(f'{path}: undeclared assembly or primary reference sequence')
        p, e, s, t = (integer(row[k], path) for k in ('ref_pos1', 'ref_end1', 'asm_start0', 'asm_end0'))
        if not 1 <= p <= e <= reference_lengths[row['ref_contig']]:
            raise ValueError(f'{path}: invalid 1-based reference anchor/end')
        if not 0 <= s <= t <= lengths[row['asm_contig']]:
            raise ValueError(f'{path}: invalid 0-based assembly interval')
        if row['sv_type'] == 'INS' and not (p == e and s < t):
            raise ValueError(f'{path}: INS requires ref_pos1 == ref_end1 and a positive assembly span')
        if row['sv_type'] == 'DEL' and not (p < e and s == t):
            raise ValueError(f'{path}: DEL requires ref_pos1 < ref_end1 and a point assembly junction')
    return records


def validate_tables(genomes_file, paths_file, hal_inventory, gfa_inventory, hal_only_exclusions=None):
    genomes = read_exact(genomes_file, GENOMES)
    paths = read_exact(paths_file, PATHS)
    root = Path(genomes_file).resolve().parent
    hal_rows = read_exact(hal_inventory, mapping.HAL_INVENTORY_COLUMNS)
    gfa_rows = read_exact(gfa_inventory, mapping.GFA_INVENTORY_COLUMNS)
    hal = {(r['hal_genome'], r['hal_sequence']): int(r['length_bp']) for r in hal_rows}
    graph = {r['gfa_path_id']: r for r in gfa_rows}
    graph_names = {}
    for r in gfa_rows:
        key = (r['record_type'], r['w_sample_id'], r['w_haplotype_index'], r['sequence_id']) if r['record_type'] == 'W' else (r['record_type'], r['raw_path_name'])
        if key in graph_names:
            raise ValueError('Duplicate literal GFA path identity')
        graph_names[key] = r['gfa_path_id']
    if len(hal) != len(hal_rows) or len(graph) != len(gfa_rows):
        raise ValueError('Duplicate sequence/path inventory entries')
    names = [r['hal_genome'] for r in genomes]
    if len(names) != len(set(names)) or any(not n for n in names):
        raise ValueError('HAL genome names must be unique and nonempty')
    if set(names) & {'locus_id', 'allele_id'}:
        raise ValueError('HAL genome names locus_id and allele_id are reserved output column names')
    roles = [r['role'] for r in genomes]
    if not set(roles) <= mapping.VALID_ROLES or roles.count('primary_reference') != 1 or roles.count('comparison_reference') > 1 or roles.count('sample') == 0:
        raise ValueError('Require one primary reference, at most one comparison reference, and at least one sample')
    for row in genomes:
        if not any(g == row['hal_genome'] for g, _ in hal):
            raise ValueError(f'No HAL sequences for {row["hal_genome"]!r}')
        for key in ('cpgi_bed', 'sv_tsv'):
            if not row[key]:
                raise ValueError(f'Missing {key} for {row["hal_genome"]!r}')
            path = Path(row[key])
            row[key] = str((root / path).resolve())
            if not Path(row[key]).is_file():
                raise FileNotFoundError(row[key])
    seen, assigned = set(), set()
    for row in paths:
        if row['gfa_record_type'] == 'W':
            if row['gfa_path_name'] or not all(row[k] for k in ('gfa_sample', 'gfa_haplotype', 'gfa_sequence')):
                raise ValueError('W mappings require sample, haplotype and sequence; path name must be empty')
            gfa_key = ('W', row['gfa_sample'], row['gfa_haplotype'], row['gfa_sequence'])
        elif row['gfa_record_type'] == 'P':
            if not row['gfa_path_name'] or any(row[k] for k in ('gfa_sample', 'gfa_haplotype', 'gfa_sequence')):
                raise ValueError('P mappings require the complete path name; W fields must be empty')
            gfa_key = ('P', row['gfa_path_name'])
        else:
            raise ValueError('gfa_record_type must be W or P')
        if gfa_key not in graph_names:
            raise ValueError(f'Unknown literal GFA path identity: {gfa_key}')
        row['gfa_path_id'] = graph_names[gfa_key]
        key = row['hal_genome'], row['hal_sequence']
        if key in seen or key not in hal or key[0] not in names:
            raise ValueError(f'Unknown, unselected or duplicate HAL sequence mapping: {key}')
        path = graph.get(row['gfa_path_id'])
        if path is None or row['gfa_path_id'] in assigned:
            raise ValueError('Unknown or multiply assigned GFA path')
        if path['overlap_status'] not in mapping.ACCEPTED_PATH_STATUSES:
            raise ValueError('Unsupported GFA path overlap representation')
        if int(path['start0']) != 0:
            raise ValueError(f'GFA path must cover the complete HAL sequence from zero: {key}')
        for field in ('end0', 'path_length_bp'):
            if path[field] and int(path[field]) != hal[key]:
                raise ValueError(f'GFA path must cover the complete HAL sequence from zero: {key}')
        seen.add(key)
        assigned.add(row['gfa_path_id'])
    expected = {key for key in hal if key[0] in names}
    excluded_rows = read_exact(hal_only_exclusions, HAL_ONLY, allow_empty=True) if hal_only_exclusions else []
    excluded = {(r['hal_genome'], r['hal_sequence']) for r in excluded_rows}
    if len(excluded) != len(excluded_rows) or not excluded <= expected or excluded & seen:
        raise ValueError('HAL-only exclusions must be unique selected HAL sequences without a mapping')
    if seen | excluded != expected:
        raise ValueError(f'Missing mappings without an explicit HAL-only declaration: {sorted(expected - seen - excluded)}')
    if {g for g, _ in seen} != set(names):
        raise ValueError('Every selected genome requires at least one mapped GFA path')
    primary = next(r['hal_genome'] for r in genomes if r['role'] == 'primary_reference')
    reference_lengths = {c: n for (g, c), n in hal.items() if g == primary}
    for row in genomes:
        lengths = {c: n for (g, c), n in hal.items() if g == row['hal_genome']}
        for bed in read_bed(row['cpgi_bed'], lengths):
            if (row['hal_genome'], bed[0]) in excluded:
                raise ValueError('HAL-only exclusion has CGI input: ' + str((row['hal_genome'], bed[0])))
        events = validate_sv(row['sv_tsv'], lengths, reference_lengths)
        if row['role'] != 'sample' and events:
            raise ValueError('Reference-role SV tables must be header-only; SV annotation is defined for sample haplotypes')
        for event in events:
            if (row['hal_genome'], event['asm_contig']) in excluded or (primary, event['ref_contig']) in excluded:
                raise ValueError('HAL-only exclusion has SV assembly or reference input: ' + event['sv_id'])
    return genomes, paths, hal


def prepare(args):
    exclusions = getattr(args, 'hal_only_exclusions', None)
    genomes, paths, hal = validate_tables(args.genomes, args.paths, args.hal_inventory, args.gfa_inventory, exclusions)
    sources = {str(Path(p).resolve()) for p in (args.gfa, args.hal, args.genomes, args.paths, args.hal_inventory, args.gfa_inventory)}
    if exclusions:
        sources.add(str(Path(exclusions).resolve()))
    sources.update(row[k] for row in genomes for k in ('cpgi_bed', 'sv_tsv'))
    states = {p: file_state(p) for p in sources}
    large_inputs = {str(Path(p).resolve()) for p in (args.gfa, args.hal)}
    hashes = {p: digest(p) for p in sorted(sources - large_inputs)}
    output = Path(args.out_dir).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    internal_genomes, internal_paths, identities = [], [], []
    for row in genomes:
        gid = storage_id('G', row['hal_genome'])
        internal_genomes.append(dict(genome_id=gid, role=row['role'], individual_id='', haplotype_id='',
            hal_genome=row['hal_genome'], cpgi_bed=row['cpgi_bed'], sv_file=row['sv_tsv'],
            sv_format='canonical_tsv', sv_schema='', sv_evidence_status='available'))
    for row in paths:
        gid = storage_id('G', row['hal_genome'])
        cid = storage_id('C', row['hal_sequence'])
        internal_paths.append(dict(genome_id=gid, contig_id=cid, gfa_path_id=row['gfa_path_id'],
            hal_sequence=row['hal_sequence'], cpgi_contig=row['hal_sequence'], sv_contig=row['hal_sequence'],
            mapping_status='confirmed', mapping_method='explicit_user_assignment', mapping_note=''))
        identities.append(dict(row, genome_id=gid, contig_id=cid))
    mapping.write_tsv(output / 'genomes.validated.tsv', internal_genomes, mapping.GENOME_COLUMNS)
    mapping.write_tsv(output / 'contigs.validated.tsv', internal_paths, mapping.CONTIG_COLUMNS)
    mapping.write_tsv(output / 'identities.tsv', identities, PATHS + ['gfa_path_id', 'genome_id', 'contig_id'])
    selected = {r['hal_genome'] for r in genomes}
    mapped = {(r['hal_genome'], r['hal_sequence']) for r in paths}
    excluded = [dict(hal_genome=g, hal_sequence=c, length_bp=hal[(g, c)],
        reason='user_declared_no_gfa_path_no_cgi_or_sv_input', cgi_rows=0, sv_assembly_rows=0, sv_reference_rows=0)
        for g, c in sorted(k for k in hal if k[0] in selected and k not in mapped)]
    mapping.write_tsv(output / 'hal_only_exclusions.tsv', excluded, HAL_ONLY_REPORT)
    expected_lengths = {r['gfa_path_id']: hal[(r['hal_genome'], r['hal_sequence'])] for r in paths}
    (output / 'path_lengths.json').write_text(json.dumps(expected_lengths, indent=2) + '\n')
    if any(file_state(p) != states[p] for p in sources):
        raise RuntimeError('Input changed during validation')
    code = Path(__file__).resolve().parent
    lock = dict(schema_version=2, status='validated', files=hashes, file_states=states,
        code={str(p.relative_to(code)): digest(p) for p in sorted(code.rglob('*')) if p.is_file() and p.suffix in ('.py', '.sh')},
        genomes=genomes, paths=paths, hal_only_exclusions=excluded,
        validation_scope='Explicit mapping and input contract; computed path lengths checked during unfolding',
        sequence_identity_checked=False,
        source_paths=dict(gfa=str(Path(args.gfa).resolve()), hal=str(Path(args.hal).resolve())),
        selected_haplotype_n=sum(r['role'] == 'sample' for r in genomes),
        internal_files={p.name: digest(p) for p in sorted(output.iterdir()) if p.suffix in ('.tsv', '.json')})
    logical = [(name, str(Path(getattr(args, name)).resolve())) for name in ('gfa','hal','genomes','paths')]
    if exclusions:
        logical.append(('hal_only_exclusions', str(Path(exclusions).resolve())))
    logical.extend((f'{r["hal_genome"]}:{field}', r[field]) for r in genomes for field in ('cpgi_bed','sv_tsv'))
    lock['logical_inputs'] = [dict(role=role, name=Path(file).name, sha256=hashes.get(file),
        integrity_method='file_metadata' if file in large_inputs else 'sha256',
        size_bytes=states[file][2], mtime_ns=states[file][3]) for role,file in logical]
    (output / 'inputs.lock.json').write_text(json.dumps(lock, indent=2, ensure_ascii=True) + '\n')
    return dict(status='pass', selected_haplotype_n=lock['selected_haplotype_n'], hal_only_excluded_n=len(excluded), output=str(output))


def main():
    p = argparse.ArgumentParser()
    for name in ('gfa', 'hal', 'genomes', 'paths', 'hal-inventory', 'gfa-inventory', 'out-dir'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--hal-only-exclusions', help='Reviewed two-column hal_genome/hal_sequence TSV declaring absent GFA paths; CGI/SV usage is rejected')
    args = p.parse_args()
    print(json.dumps(prepare(args), indent=2))


if __name__ == '__main__':
    main()
