import argparse
import csv
import json
import math
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from pancgi_annotations import RELATIONSHIPS, mechanism
from pancgi_contract import digest, open_text, parse_sv_length
from pancgi_output_schema import FIELDS, VERSION
from pancgi_validation_store import RowSequence, ValidationStore
from pancgi_evidence import EVIDENCE_SCHEMA, validate_record
from pancgi_identifiers import representative_allele_id


def require(condition, message):
    if not condition:
        raise ValueError(message)


def table(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        header = reader.fieldnames
        require(header is not None and len(header) == len(set(header)), f'{path.name}: invalid header')
        rows = list(reader)
        require(all(None not in r and None not in r.values() for r in rows), f'{path.name}: invalid field count')
        return header, rows


def unique(rows, key):
    if isinstance(rows, RowSequence):
        require(rows.key == key, 'Incorrect validation index key')
        return rows.index
    result = {r[key]: r for r in rows}
    require(len(result) == len(rows), f'Duplicate {key}')
    return result


def check_cell(value, specification, context):
    if value == '':
        require(specification['nullable'], f'{context}: missing value')
        return
    kind = specification['type']
    if kind == 'integer':
        require(value.isascii() and value.isdecimal(), f'{context}: invalid nonnegative integer')
        parsed = int(value)
    elif kind == 'signed_integer':
        digits = value[1:] if value.startswith(('-', '+')) else value
        require(digits.isascii() and digits.isdecimal(), f'{context}: invalid signed integer')
        parsed = int(value)
    elif kind == 'number':
        parsed = float(value)
        require(math.isfinite(parsed) and parsed >= 0, f'{context}: invalid number')
        if specification['unit'] == 'fraction':
            require(parsed <= 1, f'{context}: fraction outside 0..1')
    elif kind in ('json_array', 'json_object'):
        parsed = json.loads(value)
        require(isinstance(parsed, list if kind == 'json_array' else dict), f'{context}: invalid JSON container')
    else:
        parsed = value
    if 'values' in specification:
        require(parsed in specification['values'], f'{context}: invalid value {value!r}')


def validate(directory, require_complete=True):
    with tempfile.TemporaryDirectory(prefix='pancgi-validation-') as scratch:
        store = ValidationStore(Path(scratch) / 'tables.sqlite')
        try:
            return _validate(directory, require_complete, store)
        finally:
            store.close()


def _validate(directory, require_complete, store):
    root = Path(directory)
    schema = json.loads((root/'schema.json').read_text())
    require(schema['schema_version'] == VERSION, 'Unsupported results schema')
    required = ['genomes.tsv', 'pancgi.loci.tsv.gz', 'pancgi.alleles.tsv.gz', 'pancgi.members.tsv.gz',
                'pancgi.member_sv.tsv.gz', 'pancgi.locus_genotypes.tsv.gz', 'pancgi.allele_genotypes.tsv.gz']
    data = {}
    for name in required:
        def check(r):
            require(None not in r and None not in r.values(), f'{name}: invalid field count')
            for field, value in r.items():
                spec = schema['tables'][name]['fields'][field] if name.endswith('_genotypes.tsv.gz') else FIELDS[field]
                check_cell(value, spec, f'{name}:{field}')
        key = {'genomes.tsv': 'hal_genome', 'pancgi.loci.tsv.gz': 'locus_id',
               'pancgi.alleles.tsv.gz': 'allele_id', 'pancgi.members.tsv.gz': 'member_id',
               'pancgi.member_sv.tsv.gz': None, 'pancgi.locus_genotypes.tsv.gz': 'locus_id',
               'pancgi.allele_genotypes.tsv.gz': 'allele_id'}[name]
        with open_text(root/name) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            header = reader.fieldnames
            require(header is not None and len(header) == len(set(header)), f'{name}: invalid header')
            require(header == schema['tables'][name]['columns'], f'{name}: schema/header disagreement')
            rows = store.add(reader, key, check)
        data[name] = (header, rows)
    genomes = unique(data['genomes.tsv'][1], 'hal_genome')
    samples = [g for g, r in genomes.items() if r['role'] == 'sample']
    require(samples and sum(r['role'] == 'primary_reference' for r in genomes.values()) == 1, 'Invalid genome roles')
    require(sum(r['role'] == 'comparison_reference' for r in genomes.values()) <= 1, 'Multiple comparison references')
    loci = unique(data['pancgi.loci.tsv.gz'][1], 'locus_id')
    alleles = unique(data['pancgi.alleles.tsv.gz'][1], 'allele_id')
    members = unique(data['pancgi.members.tsv.gz'][1], 'member_id')
    matrices = {}
    for kind, annotation in [('locus', loci), ('allele', alleles)]:
        header, rows = data[f'pancgi.{kind}_genotypes.tsv.gz']
        require(header == [kind+'_id']+list(genomes), f'{kind}: genome order mismatch')
        matrix = unique(rows, kind+'_id')
        require(set(matrix) == set(annotation), f'{kind}: missing or extra genotype rows')
        for identifier, row in matrix.items():
            calls = [row[g] for g in samples]
            require(set(row[g] for g in genomes) <= {'0','1','NA'}, 'Invalid genotype')
            ac, n = calls.count('1'), sum(v != 'NA' for v in calls)
            a = annotation[identifier]
            expected = dict(ac=ac, n_called=n, n_missing=len(samples)-n, sample_haplotype_n=len(samples), in_sample_set=int(ac > 0))
            require(all(int(a[k]) == v for k, v in expected.items()), 'Frequency count mismatch')
            require(math.isclose(float(a['call_rate']), n/len(samples), abs_tol=1e-9), 'Call rate mismatch')
            if n:
                require(math.isclose(float(a['af']), ac/n, abs_tol=1e-9), 'AF mismatch')
                require(math.isclose(float(a['maf']), min(ac/n, 1-ac/n), abs_tol=1e-9), 'MAF mismatch')
            else:
                require(a['af'] == a['maf'] == '', 'All-missing AF/MAF must be empty')
        matrices[kind] = matrix
    by_locus, observed = defaultdict(list), set()
    for aid, allele in alleles.items():
        lid, rep = allele['locus_id'], allele['representative_member_id']
        require(lid in loci and rep in members, 'Allele foreign key missing')
        require(members[rep]['allele_id'] == aid, 'Representative belongs to another allele')
        require(aid == representative_allele_id(members[rep]), 'Allele ID does not match representative identity')
        for field in ('mechanism_group', 'mechanism_class', 'mechanism_status', 'length_bp', 'cpg_n', 'gc_n', 'pct_gc', 'observed_expected'):
            require(allele[field] == members[rep][field], f'Allele representative mismatch: {field}')
        by_locus[lid].append(aid)
    events = defaultdict(list)
    inverse = {v:k for k,v in RELATIONSHIPS.items()}
    event_keys = set()
    for event in data['pancgi.member_sv.tsv.gz'][1]:
        mid = event['member_id']
        require(mid in members, 'SV member foreign key missing')
        member = members[mid]
        require(all(event[k] == member[k] for k in ('allele_id','locus_id','hal_genome')), 'SV event keys disagree')
        key = (mid, event['sv_id'])
        require(key not in event_keys, 'Duplicate member SV event')
        event_keys.add(key)
        start, end = int(event['asm_start0']), int(event['asm_end0'])
        cs, ce = int(member['start0']), int(member['end0'])
        require(event['asm_contig'] == member['hal_sequence'], 'SV is on a different assembly sequence')
        p, e = int(event['ref_pos1']), int(event['ref_end1'])
        parse_sv_length(event['sv_length'], event['sv_type'], f'{mid}:{event["sv_id"]}')
        if event['sv_type'] == 'DEL':
            require(p < e and start == end and cs < start < ce, 'Invalid DEL-junction geometry')
            require(event['overlap_fraction_sv'] == '' and int(event['overlap_bp']) == 0, 'DEL point overlap is undefined/zero')
            require(event['relationship'] == 'DEL-junction', 'DEL relationship mismatch')
        else:
            ov = min(end, ce) - max(start, cs)
            require(p == e and ov > 0 and end > start, 'Invalid INS geometry')
            require(int(event['overlap_bp']) == ov, 'INS overlap mismatch')
            require(math.isclose(float(event['overlap_fraction_sv']), ov/(end-start), abs_tol=1e-9), 'INS fraction mismatch')
            relation = 'INS-internal' if start <= cs and ce <= end else ('INS-spanning' if cs <= start and end <= ce else 'INS-junction')
            require(event['relationship'] == relation, 'INS relationship mismatch')
        require(math.isclose(float(event['overlap_fraction_cgi']), int(event['overlap_bp'])/(ce-cs), abs_tol=1e-9), 'CGI overlap fraction mismatch')
        events[mid].append(dict(svtype=event['sv_type'], **{'class':inverse[event['relationship']]}))
    for mid, member in members.items():
        aid, lid, g = (member[k] for k in ('allele_id','locus_id','hal_genome'))
        require(aid in alleles and alleles[aid]['locus_id'] == lid and g in genomes, 'Member foreign keys disagree')
        require(member['role'] == genomes[g]['role'], 'Member role mismatch')
        require(int(member['end0'])-int(member['start0']) == int(member['length_bp']) > 0, 'Member length mismatch')
        require(matrices['allele'][aid][g] == '1', 'Observed member is not genotype 1')
        observed.add((aid, g))
        annotation = mechanism(events[mid], member['role'])
        require(all(member[k] == v for k,v in annotation.items()), 'Member mechanism contradicts event table')
        for ref in json.loads(member['all_ins_matched_reference_members']):
            require(ref in members and members[ref]['role'] == 'primary_reference', 'Invalid reference member link')
    for aid, row in matrices['allele'].items():
        for g in genomes:
            require((row[g] == '1') == ((aid, g) in observed), 'Positive genotype lacks member evidence')
    for lid, locus in loci.items():
        require(by_locus[lid], 'Locus has no alleles')
        for g in genomes:
            calls = [matrices['allele'][aid][g] for aid in by_locus[lid]]
            value = matrices['locus'][lid][g]
            require((value == '1' and '1' in calls and 'NA' not in calls) or
                    (value == '0' and set(calls) == {'0'}) or
                    (value == 'NA' and set(calls) == {'NA'}), 'Locus/allele propagation mismatch')
        expected = Counter(alleles[aid]['mechanism_class'] or 'not_applicable_reference_role' for aid in by_locus[lid])
        require(dict(expected) == json.loads(locus['allele_mechanism_counts']), 'Locus mechanism composition mismatch')
        if locus['anchor_member_id']:
            require(locus['anchor_member_id'] in members, 'Missing locus anchor member')
    evidence_keys, reasons = set(), Counter()
    require(schema.get('jsonl') == {'pancgi.genotype_evidence.jsonl.gz': EVIDENCE_SCHEMA}, 'Evidence schema mismatch')
    run = json.loads((root/'run_manifest.json').read_text())
    require(run['schema_version'] == VERSION, 'Run manifest schema mismatch')
    require(isinstance(run['sequences'], list), 'Invalid sequence inventory')
    sequences = set()
    for item in run['sequences']:
        require(isinstance(item, dict) and set(item) == {'hal_genome', 'hal_sequence'}, 'Invalid sequence inventory fields')
        require(item['hal_genome'] in genomes and isinstance(item['hal_sequence'], str) and item['hal_sequence'], 'Unknown sequence identity')
        key = item['hal_genome'], item['hal_sequence']
        require(key not in sequences, 'Duplicate sequence identity')
        sequences.add(key)
    margin = run['genotype_configuration']['margin']
    require(type(margin) in (int, float) and math.isfinite(margin) and 0 <= margin <= 1, 'Invalid placement dominance margin')
    with open_text(root/'pancgi.genotype_evidence.jsonl.gz') as handle:
        for line in handle:
            record = json.loads(line)
            validate_record(record, loci, alleles, members, genomes, sequences, margin, by_locus)
            lid, g = record['locus_id'], record['hal_genome']
            require(lid in loci and g in genomes, 'Unknown genotype evidence key')
            require((lid, g) not in evidence_keys, 'Duplicate genotype evidence key')
            require(str(record['genotype']) == matrices['locus'][lid][g], 'Genotype evidence contradicts matrix')
            evidence_keys.add((lid, g))
            reasons[(record['genotype'], record['evidence']['reason'])] += 1
    require(len(evidence_keys) == len(loci)*len(genomes), 'Incomplete genotype evidence')
    qc = json.loads((root/'qc.json').read_text())
    require(qc['accepted_member_n'] == len(members), 'QC member count mismatch')
    require(qc['input_member_n'] == qc['accepted_member_n']+qc['excluded_member_n'], 'QC counts do not close')
    require(sum(qc['excluded_reasons'].values()) == qc['excluded_member_n'], 'QC exclusion reasons do not close')
    require(qc['genotype_evidence_n'] == len(evidence_keys), 'QC genotype evidence count mismatch')
    require(qc['genotype_reasons'] == [dict(genotype=gt, reason=reason, cells=count) for (gt,reason),count in sorted(reasons.items())], 'QC genotype reasons mismatch')
    seqs = {}
    current = None
    with open_text(root/'pancgi.alleles.fa.gz') as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith('>'):
                current = line[1:]
                require(current not in seqs, 'Duplicate FASTA identifier')
                seqs[current] = 0
            elif line:
                require(current is not None, 'Sequence before FASTA header')
                seqs[current] += len(line)
    require(set(seqs) == set(alleles), 'FASTA/catalogue keys mismatch')
    require(all(n == int(alleles[aid]['length_bp']) for aid,n in seqs.items()), 'FASTA length mismatch')
    counts = dict(locus_n=len(loci), allele_n=len(alleles), member_n=len(members), genome_n=len(genomes), sample_haplotype_n=len(samples))
    if require_complete:
        complete = json.loads((root/'completed.json').read_text())
        require(complete['status'] == 'complete', 'Incomplete run')
        require(all(complete[k] == v for k,v in counts.items()), 'Completion counts disagree')
        require(set(required+['pancgi.alleles.fa.gz','schema.json','qc.json','pancgi.genotype_evidence.jsonl.gz','run_manifest.json','validation.json']) <= set(complete['output_sha256']), 'Missing required output checksums')
        for name, sha in complete['output_sha256'].items():
            require(Path(name).name == name, 'Nonlocal checksum filename')
            require(digest(root/name) == sha, f'Output checksum mismatch: {name}')
    return dict(status='pass', schema_version=VERSION, **counts)


def main():
    p = argparse.ArgumentParser(description='Validate PanCGI result schemas, relationships, genotypes, counts and checksums.')
    p.add_argument('results')
    args = p.parse_args()
    print(json.dumps(validate(args.results), indent=2))


if __name__ == '__main__':
    main()
