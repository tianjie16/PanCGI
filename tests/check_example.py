import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
expected = json.loads(Path(sys.argv[2]).read_text())
complete = json.loads((root/'completed.json').read_text())
for key in ('genome_n','sample_haplotype_n','locus_n','allele_n','member_n'):
    if complete[key] != expected[key]:
        raise ValueError(f'Example count mismatch: {key}')
with gzip.open(root/'pancgi.members.tsv.gz','rt') as handle:
    counts = Counter(r['mechanism_class'] for r in csv.DictReader(handle, delimiter='\t') if r['role']=='sample')
if dict(counts) != expected['sample_member_mechanisms']:
    raise ValueError(f'Example mechanism mismatch: {counts}')
with gzip.open(root/'pancgi.alleles.tsv.gz', 'rt') as handle:
    allele_ids = [r['allele_id'] for r in csv.DictReader(handle, delimiter='\t')]
if allele_ids != expected['allele_ids']:
    raise ValueError(f'Example allele identifier or order mismatch: {allele_ids}')
print(json.dumps({'status':'pass','scope':expected['scope']}))
