import argparse
import csv
import gzip
import json
import importlib.metadata
import platform
import os
import gc
from collections import Counter, defaultdict
from pathlib import Path

import pancgi_mapping as mapping
from pancgi_contract import digest, open_text, verify_lock
from pancgi_annotations import mechanism, RELATIONSHIPS
from pancgi_output_schema import export_schema, VERSION
from pancgi_validate_results import validate
from pancgi_evidence import export_record
from pancgi_public_parameters import public_parameters
from pancgi_identifiers import allele_id_map


def read(path):
    return list(iter_read(path))


def iter_read(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write(path, rows, columns):
    with open_text(path, 'wt') as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


class TableWriter:
    def __init__(self, path, columns):
        self.handle = open_text(path, 'wt')
        self.writer = csv.DictWriter(self.handle, columns, delimiter='\t', lineterminator='\n', extrasaction='raise')
        self.writer.writeheader()
        self.count = 0

    def append(self, row):
        self.writer.writerow(row)
        self.count += 1

    def close(self):
        self.handle.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def __len__(self):
        return self.count


def relative_strand(value):
    if value in ('', '+', '-'):
        return value
    if value in ('++', '--', '+-', '-+'):
        return '+' if value[0] == value[1] else '-'
    raise ValueError(f'Invalid PSL strand: {value!r}')


def main():
    parser = argparse.ArgumentParser()
    for key in ['internal-results', 'validated', 'inputs', 'out-dir']:
        parser.add_argument('--'+key, required=True)
    args = parser.parse_args()
    root, validated, inputs, output = map(Path, [args.internal_results, args.validated, args.inputs, args.out_dir])
    lock = verify_lock(validated / 'inputs.lock.json')
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    genomes = read(validated / 'genomes.validated.tsv')
    labels = {r['genome_id']: r['hal_genome'] for r in genomes}
    roles = {r['hal_genome']: r['role'] for r in genomes}
    samples = [r['hal_genome'] for r in genomes if r['role'] == 'sample']
    members = {}
    for i, r in enumerate(iter_read(inputs / 'member_identity.tsv'), 1):
        if r['fid'] in members:
            raise ValueError('Duplicate member identity')
        members[r['fid']] = dict(member_id=f'M{i:09d}', **{k: v for k, v in r.items() if k != 'fid'})
    identities = {r['contig_id']: r for r in read(validated / 'identities.tsv')}
    sequences = {key: r['hal_sequence'] for key, r in identities.items()}
    primary = next(name for name, role in roles.items() if role == 'primary_reference')
    allele_catalog = read(root / 'allele_catalog.tsv.gz')
    allele_ids = allele_id_map(allele_catalog, members)
    statistics = {}
    for kind in ['locus', 'allele']:
        filename = root / f'pancgi.{kind}_genotype_strict.tsv.gz'
        with open_text(filename) as handle:
            columns = next(csv.reader(handle, delimiter='\t'))
        if columns != [kind+'_id'] + list(labels):
            raise ValueError('Genotype columns do not match validated genome order')
        stats = {}
        def public_rows():
            for row in iter_read(filename):
                identifier = row[kind+'_id']
                if identifier in stats:
                    raise ValueError('Duplicate genotype identifier')
                values = {labels[g]: row[g] for g in labels}
                if not set(values.values()) <= {'0', '1', 'NA'}:
                    raise ValueError('Invalid genotype')
                calls = [values[g] for g in samples]
                ac = calls.count('1')
                n = ac+calls.count('0')
                af = ac/n if n else None
                stats[identifier] = dict(sample_haplotype_n=len(samples), ac=ac, n_called=n,
                    n_missing=len(samples)-n, call_rate=n/len(samples), af=af,
                    maf=min(af, 1-af) if af is not None else None, in_sample_set=int(ac > 0))
                public_id = allele_ids[identifier] if kind == 'allele' else identifier
                yield {kind+'_id': public_id, **values}
        write(output / f'pancgi.{kind}_genotypes.tsv.gz', public_rows(), [kind+'_id']+list(labels.values()))
        statistics[kind] = stats
    by_member = {}
    anchor_members = {r['fid']: r for r in iter_read(root / 'locus_anchored_members.tsv.gz')}
    for row in iter_read(root / 'allele_members.tsv.gz'):
        fid = row['fid']
        if fid in by_member:
            raise ValueError('Member assigned to multiple alleles')
        by_member[fid] = (row['locus_id'], allele_ids[row['allele_id']])
    feature_file = root.parent / 'features.jsonl.gz'
    loci = read(root / 'locus_anchored.tsv.gz')
    wanted_features = {r['allele_rep_fid'] for r in allele_catalog} | {r['anchor_ref_fid'] for r in loci if r['anchor_ref_fid']}
    features, accepted_ids = {}, set()
    member_columns = ['member_id', 'hal_genome', 'hal_sequence', 'start0', 'end0', 'input_cgi_name',
        'role', 'locus_id', 'allele_id', 'length_bp', 'cpg_n', 'gc_n', 'pct_cpg', 'pct_gc', 'observed_expected',
        'graph_coverage', 'graph_steps', 'mechanism_group', 'mechanism_class', 'mechanism_status', 'hal_mapping_status', 'primary_sequence',
        'primary_start0', 'primary_end0', 'primary_strand', 'hal_query_coverage', 'hal_identity', 'reference_ins_id',
        'longest_overlapping_ins_id', 'anchor_assignment', 'all_ins_matched_reference_members']
    member_rows = TableWriter(output / 'pancgi.members.tsv.gz', member_columns)
    overlaps = TableWriter(output / 'pancgi.member_sv.tsv.gz', ['member_id', 'locus_id', 'allele_id', 'hal_genome',
        'sv_id', 'sv_type', 'ref_contig', 'ref_pos1', 'ref_end1', 'asm_contig', 'asm_start0', 'asm_end0', 'sv_length',
        'relationship', 'overlap_bp', 'overlap_fraction_cgi', 'overlap_fraction_sv'])
    with member_rows, overlaps, open_text(feature_file) as handle:
        for raw in handle:
            f = json.loads(raw)
            fid = f['fid']
            if fid in accepted_ids:
                raise ValueError('Duplicate accepted feature')
            accepted_ids.add(fid)
            if fid in wanted_features:
                features[fid] = f
            lid, aid = by_member[fid]
            member = members[fid]
            role = roles[member['hal_genome']]
            events = json.loads(f['sv_overlap_detail_json'])
            row = dict(member, role=role, locus_id=lid, allele_id=aid,
                       length_bp=f['source_seq_len'], cpg_n=f['cpg_n'], gc_n=f['gc_n'],
                       pct_cpg=f['pct_cpg'], pct_gc=f['pct_gc'], observed_expected=f['oe'],
                       graph_coverage=f['graph_cov'], graph_steps=f['graph_nsteps'],
                       **mechanism(events, role), hal_mapping_status=f['hal_multimap_label'],
                       primary_sequence=sequences[f['hal_primary_chr']] if f['hal_primary_chr'] else '',
                       primary_start0=f['hal_primary_start0'], primary_end0=f['hal_primary_end0'],
                       primary_strand=relative_strand(f['hal_strand']),
                       hal_query_coverage=f['hal_query_cov'], hal_identity=f['hal_identity'])
            anchor = anchor_members[fid]
            row.update(reference_ins_id=anchor['reference_ins_id'], longest_overlapping_ins_id=f['sv_ins_longest_id'],
                anchor_assignment=anchor['anchor_assign_rule'],
                all_ins_matched_reference_members=json.dumps([members[n]['member_id'] for n in json.loads(anchor['all_ins_matched_ref_ids_json'])]))
            member_rows.append(row)
            for sv in events:
                overlaps.append(dict(member_id=member['member_id'], locus_id=lid, allele_id=aid,
                    hal_genome=member['hal_genome'], sv_id=sv['id'], sv_type=sv['svtype'],
                    ref_contig=sequences[sv['chrom']], ref_pos1=sv['pos1'], ref_end1=sv['vcf_end1'],
                    asm_contig=sequences[sv['contig']], asm_start0=sv['asm_start0'], asm_end0=sv['asm_end0'],
                    sv_length=sv['svlen'],
                    relationship=RELATIONSHIPS[sv['class']], overlap_bp=sv['overlap_bp'],
                    overlap_fraction_cgi=sv['overlap_pct_cpgi'], overlap_fraction_sv=sv['overlap_pct_sv']))
    if accepted_ids != set(by_member):
        raise ValueError('Accepted features and final allele members disagree')
    allele_mechanisms = {}
    composition = defaultdict(Counter)
    for allele in allele_catalog:
        f = features[allele['allele_rep_fid']]
        role = roles[members[f['fid']]['hal_genome']]
        annotation = mechanism(json.loads(f['sv_overlap_detail_json']), role)
        allele_mechanisms[allele['allele_id']] = annotation
        composition[allele['locus_id']][annotation['mechanism_class'] or 'not_applicable_reference_role'] += 1
    locus_rows = []
    for r in loci:
        fid = r['anchor_ref_fid']
        if fid:
            f = features[fid]
            sequence, start, end = sequences[f['contig']], f['asm_start0'], f['asm_end0']
        else:
            sequence = sequences[r['nonref_site_chrom']] if r['nonref_site_chrom'] else ''
            start = int(r['nonref_site_start1'])-1 if r['nonref_site_start1'] else ''
            end = int(r['nonref_site_end1']) if r['nonref_site_end1'] else ''
        locus_rows.append(dict(locus_id=r['locus_id'], locus_type={'primary_ref':'Reference', 'nonref':'Novel'}[r['anchor_type']],
            primary_sequence=sequence, primary_start0=start, primary_end0=end,
            coordinate_kind='reference_cgi_interval' if fid else ('novel_site_cluster' if sequence else 'unavailable'),
            coordinate_source='primary_reference_cgi' if fid else r['nonref_site_source'],
            allele_mechanism_counts=json.dumps(dict(sorted(composition[r['locus_id']].items()))),
            anchor_member_id=members[fid]['member_id'] if fid else '', **statistics['locus'][r['locus_id']]))
    stat_columns = ['sample_haplotype_n', 'ac', 'n_called', 'n_missing', 'call_rate', 'af', 'maf', 'in_sample_set']
    write(output / 'pancgi.loci.tsv.gz', locus_rows,
          ['locus_id', 'locus_type', 'primary_sequence', 'primary_start0', 'primary_end0', 'coordinate_kind',
           'coordinate_source', 'allele_mechanism_counts', 'anchor_member_id']+stat_columns)
    allele_rows = []
    with gzip.open(output / 'pancgi.alleles.fa.gz', 'wt') as fasta:
        for r in allele_catalog:
            f = features[r['allele_rep_fid']]
            member = members[r['allele_rep_fid']]
            public_id = allele_ids[r['allele_id']]
            allele_rows.append(dict(allele_id=public_id, locus_id=r['locus_id'], representative_member_id=member['member_id'],
                length_bp=f['source_seq_len'], cpg_n=f['cpg_n'], gc_n=f['gc_n'], pct_gc=f['pct_gc'],
                observed_expected=f['oe'], **allele_mechanisms[r['allele_id']], **statistics['allele'][r['allele_id']]))
            fasta.write('>'+public_id+'\n'+f['source_seq']+'\n')
    write(output / 'pancgi.alleles.tsv.gz', allele_rows,
          ['allele_id', 'locus_id', 'representative_member_id', 'length_bp', 'cpg_n', 'gc_n', 'pct_gc', 'observed_expected',
           'mechanism_group', 'mechanism_class', 'mechanism_status']+stat_columns)
    write(output / 'genomes.tsv', [dict(hal_genome=n, role=roles[n]) for n in labels.values()], ['hal_genome', 'role'])
    export_schema(output, list(labels.values()))
    excluded = read(root.parent/'features.excluded.tsv.gz')
    excluded_ids = [r['fid'] for r in excluded]
    if len(set(excluded_ids)) != len(excluded_ids) or set(excluded_ids) & accepted_ids or set(excluded_ids) | accepted_ids != set(members):
        raise ValueError('Input/accepted/excluded member counts do not close')
    reasons = Counter()
    evidence_n = 0
    public_loci = {r['locus_id']: r for r in locus_rows}
    with gzip.open(output/'pancgi.genotype_evidence.jsonl.gz', 'wt') as evidence_output:
        for gid in labels:
            with gzip.open(root/'pancgi.genotyping'/f'{gid}.evidence.jsonl.gz', 'rt') as evidence_input:
                for line in evidence_input:
                    record = json.loads(line)
                    if record['hal_genome'] != labels[gid] or record['locus_id'] not in statistics['locus']:
                        raise ValueError('Genotype evidence has an unknown key')
                    public_record = export_record(record, members, by_member, public_loci, identities, primary)
                    evidence_output.write(json.dumps(public_record, separators=(',', ':'), allow_nan=False)+'\n')
                    reasons[(record['genotype'], record['evidence']['reason'])] += 1
                    evidence_n += 1
    if evidence_n != len(locus_rows)*len(labels):
        raise ValueError('Genotype evidence cell count mismatch')
    qc = dict(input_member_n=len(members), accepted_member_n=len(accepted_ids), excluded_member_n=len(excluded),
              excluded_reasons=dict(Counter(r['reason'] for r in excluded)),
              genotype_evidence_n=evidence_n,
              genotype_reasons=[dict(genotype=gt, reason=reason, cells=count) for (gt,reason),count in sorted(reasons.items())],
              hal_only_exclusions=lock['hal_only_exclusions'])
    (output/'qc.json').write_text(json.dumps(qc, indent=2)+'\n')
    profile = json.loads((Path(__file__).parent/'RELEASE_MANIFEST.json').read_text())
    run = dict(schema_version=VERSION, package=profile, logical_inputs=lock['logical_inputs'],
               source_sha256=lock['code'], genotype_configuration=json.loads((root/'pancgi.genotyping/completed.json').read_text())['config'],
               python_version=platform.python_version(), dependency_versions={name:importlib.metadata.version(name)
                   for name in ('numpy','pandas','pyarrow','biopython','parasail','pywfa')})
    run['effective_stage_parameters'] = public_parameters(root.parent/'parameters')
    run['sequences'] = [dict(hal_genome=r['hal_genome'], hal_sequence=r['hal_sequence']) for r in identities.values()]
    run['hal_runtime'] = os.environ.get('PANCGI_HAL_RUNTIME', 'docker')
    (output/'run_manifest.json').write_text(json.dumps(run, indent=2)+'\n')
    del members, anchor_members, by_member, features, accepted_ids, excluded_ids, excluded, statistics
    gc.collect()
    validation = validate(output, require_complete=False)
    (output/'validation.json').write_text(json.dumps(validation, indent=2)+'\n')
    report = dict(status='complete', sample_haplotype_n=len(samples), genome_n=len(labels),
                  locus_n=len(locus_rows), allele_n=len(allele_rows), member_n=len(member_rows),
                  inputs_lock_sha256=digest(validated / 'inputs.lock.json'),
                  output_sha256={p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()})
    (output / 'completed.json.partial').write_text(json.dumps(report, indent=2)+'\n')
    (output / 'completed.json.partial').replace(output / 'completed.json')


if __name__ == '__main__':
    main()
