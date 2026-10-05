from __future__ import annotations

import argparse
import json
from pathlib import Path

import pancgi_mapping as mapping
from pancgi_contract import SV, digest, open_text, parse_sv_length, read_exact

CATALOG_COLUMNS = ['label', 'kind', 'graph_sample', 'graph_hap', 'hal_genome', 'bed', 'cpgi_fa', 'path_bed_dir', 'sv_tsv']
SV_CANONICAL_COLUMNS = SV
SV_OUTPUT_COLUMNS = ['ID', 'VARID', 'contig', 'start', 'end', 'SVTYPE', '#CHROM', 'POS', 'VCF_END', 'SVLEN']


def prepare_inputs(args):
    genomes = mapping.read_tsv(args.genomes, mapping.GENOME_COLUMNS)
    contigs = mapping.read_tsv(args.contigs, mapping.CONTIG_COLUMNS)
    lock_path = Path(args.genomes).parent / 'inputs.lock.json'
    lock = json.loads(lock_path.read_text())
    if lock['status'] != 'validated':
        raise ValueError('Inputs were not validated')
    for file in (args.genomes, args.contigs):
        if digest(file) != lock['internal_files'][Path(file).name]:
            raise ValueError('Prepared mapping has changed')
    root = Path(args.out_dir).resolve()
    if root.exists():
        raise FileExistsError(root)
    root.mkdir(parents=True)
    pathbed = Path(args.pathbed_dir).resolve()
    if not (pathbed / 'pathbed_inventory.tsv').is_file():
        raise FileNotFoundError(pathbed / 'pathbed_inventory.tsv')
    index = {(r['genome_id'], r['hal_sequence']): r['contig_id'] for r in contigs}
    primary = next(r['genome_id'] for r in genomes if r['role'] == 'primary_reference')
    catalog, qc = [], []
    kinds = {'primary_reference': 'reference_primary', 'comparison_reference': 'reference_comparison', 'sample': 'assembly'}
    def member_rows():
        for row in genomes:
            gid = row['genome_id']
            for key in ('cpgi_bed', 'sv_file'):
                if digest(row[key]) != lock['files'][row[key]]:
                    raise ValueError(f'Source changed after validation: {row[key]}')
            bed = root / (gid + '.bed')
            n = 0
            with open_text(row['cpgi_bed']) as source, bed.open('w') as out:
                for raw in source:
                    fields = raw.rstrip('\r\n').split('\t')
                    if len(fields) != 10:
                        raise ValueError('Prepared CGI source is not BED10')
                    literal = fields[0]
                    cid = index[(gid, literal)]
                    fields[0] = f'{gid}#0#{cid}'
                    fid = f'{fields[0]}:{fields[1]}-{fields[2]}'
                    out.write('\t'.join(fields) + '\n')
                    yield dict(fid=fid, hal_genome=row['hal_genome'], hal_sequence=literal,
                        start0=fields[1], end0=fields[2], input_cgi_name=fields[3])
                    n += 1
            sv = root / (gid + '.sv.tsv')
            sv_records = []
            for r in read_exact(row['sv_file'], SV, allow_empty=True):
                p, e, s, t = (int(r[k]) for k in ('ref_pos1', 'ref_end1', 'asm_start0', 'asm_end0'))
                sv_records.append({'ID': r['sv_id'], 'VARID': r['sv_id'], 'contig': index[(gid, r['asm_contig'])],
                    'start': s, 'end': t, 'SVTYPE': r['sv_type'], '#CHROM': index[(primary, r['ref_contig'])],
                    'POS': p, 'VCF_END': e, 'SVLEN': parse_sv_length(r['sv_length'], r['sv_type'], r['sv_id'])})
            mapping.write_tsv(sv, sv_records, SV_OUTPUT_COLUMNS)
            catalog.append(dict(label=gid, kind=kinds[row['role']], graph_sample=gid, graph_hap='0',
                hal_genome=row['hal_genome'], bed=str(bed), cpgi_fa=str(root / (gid + '.fa.gz')),
                path_bed_dir=str(pathbed), sv_tsv=str(sv)))
            qc.append(dict(hal_genome=row['hal_genome'], role=row['role'], cpgi_n=n, sv_n=len(sv_records)))
    mapping.write_tsv(root / 'member_identity.tsv', member_rows(), ['fid', 'hal_genome', 'hal_sequence', 'start0', 'end0', 'input_cgi_name'])
    mapping.write_tsv(root / 'sample_catalog.tsv', catalog, CATALOG_COLUMNS)
    mapping.write_tsv(root / 'genome_input_qc.tsv', qc, ['hal_genome', 'role', 'cpgi_n', 'sv_n'])
    report = dict(status='pass', n_genomes=len(catalog), n_samples=sum(r['kind'] == 'assembly' for r in catalog),
                  inputs_lock_sha256=digest(lock_path), files={p.name: digest(p) for p in root.iterdir() if p.is_file()})
    (root / 'prepared_inputs.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    p = argparse.ArgumentParser()
    for option in ('genomes', 'contigs', 'out-dir', 'pathbed-dir'):
        p.add_argument('--' + option, required=True)
    print(json.dumps(prepare_inputs(p.parse_args()), indent=2))


if __name__ == '__main__':
    main()
