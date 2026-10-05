import csv
import json
from pathlib import Path

from pancgi_annotations import MECHANISMS
from pancgi_contract import open_text
from pancgi_evidence import EVIDENCE_SCHEMA

VERSION = '1.3'
FIELDS = {}
PRIMARY_KEYS = {
    'genomes.tsv': ['hal_genome'],
    'pancgi.loci.tsv.gz': ['locus_id'],
    'pancgi.alleles.tsv.gz': ['allele_id'],
    'pancgi.members.tsv.gz': ['member_id'],
    'pancgi.member_sv.tsv.gz': ['member_id', 'sv_id'],
    'pancgi.locus_genotypes.tsv.gz': ['locus_id'],
    'pancgi.allele_genotypes.tsv.gz': ['allele_id']
}
FOREIGN_KEYS = {
    'pancgi.members.tsv.gz': {'locus_id':'pancgi.loci.tsv.gz:locus_id', 'allele_id':'pancgi.alleles.tsv.gz:allele_id', 'hal_genome':'genomes.tsv:hal_genome'},
    'pancgi.alleles.tsv.gz': {'locus_id':'pancgi.loci.tsv.gz:locus_id', 'representative_member_id':'pancgi.members.tsv.gz:member_id'},
    'pancgi.loci.tsv.gz': {'anchor_member_id':'pancgi.members.tsv.gz:member_id'},
    'pancgi.member_sv.tsv.gz': {'member_id':'pancgi.members.tsv.gz:member_id', 'hal_genome':'genomes.tsv:hal_genome', 'allele_id':'pancgi.alleles.tsv.gz:allele_id', 'locus_id':'pancgi.loci.tsv.gz:locus_id'},
    'pancgi.locus_genotypes.tsv.gz': {'locus_id':'pancgi.loci.tsv.gz:locus_id'},
    'pancgi.allele_genotypes.tsv.gz': {'allele_id':'pancgi.alleles.tsv.gz:allele_id'}
}


def define(names, dtype, definition, unit='', nullable=False, values=None):
    for name in names.split():
        FIELDS[name] = dict(type=dtype, definition=definition, unit=unit, nullable=nullable)
        if values is not None:
            FIELDS[name]['values'] = values


define('locus_id member_id', 'string', 'Run-scoped unique identity at the named level.')
define('allele_id', 'string', 'Representative CGI identity: {hal_genome}_{hal_sequence}:{start0}:{end0}, using its source assembly 0-based half-open interval. Each name is UTF-8 percent-encoded except ASCII letters, digits and -._~. Names retain underscores; use the representative member fields rather than splitting this ID. Collisions are rejected. IDs may change when the representative changes.')
define('representative_member_id', 'string', 'Member selected by the sequence-medoid rule; allele properties and mechanism use this member.')
define('anchor_member_id', 'string', 'Primary-reference CGI member defining this reference locus.', nullable=True)
define('hal_genome', 'string', 'Literal HAL genome name, representing one haplotype or reference.')
define('hal_sequence', 'string', 'Literal source HAL sequence name.')
define('role', 'string', 'Genome role; frequency denominators include sample roles only.', values=['sample','primary_reference','comparison_reference'])
define('input_cgi_name', 'string', 'Original BED10 name; not necessarily globally unique.', nullable=True)
define('locus_type', 'string', 'Association to a CGI of the selected primary reference.', values=['Reference','Novel'])
define('start0', 'integer', 'Source CGI inclusive 0-based start.', 'bp')
define('end0', 'integer', 'Source CGI exclusive end.', 'bp')
define('length_bp', 'integer', 'Assembly CGI sequence length; representative member for allele rows.', 'bp')
define('cpg_n', 'integer', 'Supplied caller count of CpG dinucleotides.', 'dinucleotides')
define('gc_n', 'integer', 'Supplied caller count of C plus G bases.', 'bases')
define('pct_cpg', 'number', 'Percentage of bases in CpG dinucleotides: 200*cpg_n/length_bp.', 'percent')
define('pct_gc', 'number', 'Supplied caller GC percentage: 100*gc_n/length_bp.', 'percent')
define('observed_expected', 'number', 'Supplied caller observed/expected CpG ratio.')
define('graph_coverage', 'number', 'Fraction of the source CGI covered by extracted graph steps.', 'fraction')
define('graph_steps', 'integer', 'Number of graph steps covering the source CGI.', 'steps')
define('mechanism_group', 'string', 'Qualified supplied INS/DEL evidence on this member, or on the representative for an allele.', nullable=True, values=['SV-related','non-SV'])
define('mechanism_class', 'string', 'Five exclusive SV classes or no qualified event in the supplied INS/DEL annotation.', nullable=True, values=list(MECHANISMS))
define('mechanism_status', 'string', 'SV evaluation status; reference-role mechanism fields are empty.', values=['evaluated_supplied_ins_del','not_applicable_reference_role'])
define('hal_mapping_status', 'string', 'HAL projection resolution/QC status.')
define('primary_sequence', 'string', 'Literal primary-reference HAL sequence, when assigned.', nullable=True)
define('primary_start0', 'integer', 'Inclusive 0-based primary start. Member: HAL projection. Locus: see coordinate_kind.', 'bp', True)
define('primary_end0', 'integer', 'Exclusive primary end. A novel-site cluster is not necessarily a CGI-length interval.', 'bp', True)
define('primary_strand', 'string', 'Resolved HAL relative strand on the primary reference.', nullable=True, values=['+','-'])
define('hal_query_coverage', 'number', 'HAL aligned query fraction.', 'fraction', True)
define('hal_identity', 'number', 'HAL projection identity statistic.', 'fraction', True)
define('reference_ins_id', 'string', 'Longest reference-CGI-matching INS selected for reference association, within this genome.', nullable=True)
define('longest_overlapping_ins_id', 'string', 'Globally longest overlapping INS, which need not be the reference-matching INS.', nullable=True)
define('anchor_assignment', 'string', 'Applied primary-reference association rule.')
define('all_ins_matched_reference_members', 'json_array', 'IDs of primary-reference members hit by all qualified INS anchors.')
define('coordinate_kind', 'string', 'Locus coordinate type: reference CGI interval, novel-site cluster extent or unavailable.', values=['reference_cgi_interval','novel_site_cluster','unavailable'])
define('coordinate_source', 'string', 'Primary-reference CGI or recorded novel-site grouping evidence source.', nullable=True)
define('allele_mechanism_counts', 'json_object', 'Counts of alleles per representative-member mechanism, including not_applicable_reference_role.')
define('sample_haplotype_n', 'integer', 'Number of manifest sample haplotypes, excluding both reference roles.', 'haplotypes')
define('ac', 'integer', 'Sample haplotypes with genotype 1, not member/copy count.', 'haplotypes')
define('n_called', 'integer', 'Sample haplotypes with genotype 0 or 1.', 'haplotypes')
define('n_missing', 'integer', 'Sample haplotypes with NA.', 'haplotypes')
define('call_rate', 'number', 'n_called/sample_haplotype_n.', 'fraction')
define('af', 'number', 'ac/n_called; empty when all sample haplotypes are missing.', 'fraction', True)
define('maf', 'number', 'min(af,1-af); empty when AF is undefined.', 'fraction', True)
define('in_sample_set', 'integer', '1 if at least one sample haplotype has an observed CGI; 0 for reference-only catalogue entries.', values=[0,1])
define('sv_id', 'string', 'Exact supplied SV ID, unique within hal_genome; null and None are valid literal IDs.')
define('sv_type', 'string', 'Supplied decomposed event type.', values=['INS','DEL'])
define('ref_contig', 'string', 'Literal primary-reference sequence containing the SV anchor.')
define('ref_pos1', 'integer', '1-based left padding base. Deleted bases start after this base.', 'bp')
define('ref_end1', 'integer', 'INS: equal to ref_pos1. DEL: last deleted reference base, 1-based inclusive.', 'bp')
define('asm_contig', 'string', 'Literal assembly HAL sequence containing the insertion or deletion junction.')
define('asm_start0', 'integer', '0-based inserted sequence start, or deletion-junction interbase coordinate.', 'bp')
define('asm_end0', 'integer', 'Insertion exclusive end, or the same deletion-junction coordinate as asm_start0.', 'bp')
define('sv_length', 'signed_integer', 'Original supplied SVLEN: positive for INS, negative for DEL; not derived from assembly span or reference coordinates. INS ranking uses its absolute value.', 'bp')
define('relationship', 'string', 'Individual member-event geometry; a member can have multiple relationships.', values=list(MECHANISMS[:4]))
define('overlap_bp', 'integer', 'Positive assembly overlap with INS; zero for the strictly internal DEL point.', 'bp')
define('overlap_fraction_cgi', 'number', 'overlap_bp divided by assembly CGI length.', 'fraction')
define('overlap_fraction_sv', 'number', 'overlap_bp divided by inserted assembly span; undefined for DEL junctions.', 'fraction', True)


def export_schema(directory, genome_names):
    tables = {}
    for path in sorted(Path(directory).glob('*.tsv*')):
        with open_text(path) as handle:
            columns = next(csv.reader(handle, delimiter='\t'))
        definitions = {}
        for name in columns:
            if path.name.endswith('_genotypes.tsv.gz') and name in genome_names:
                definitions[name] = dict(type='string', nullable=False, values=['0','1','NA'],
                    definition='Haplotype genotype: observed CGI=1; graph-callable absence=0; unresolved homology=NA.', unit='state')
            else:
                definitions[name] = FIELDS[name]
        tables[path.name] = dict(columns=columns, fields=definitions,
            primary_key=PRIMARY_KEYS[path.name], foreign_keys=FOREIGN_KEYS.get(path.name, {}))
    schema = dict(schema_version=VERSION, delimiter='tab', missing_value='', genotype_missing_value='NA',
                  coordinates='Coordinate conventions are declared per field.', tables=tables,
                  jsonl={'pancgi.genotype_evidence.jsonl.gz': EVIDENCE_SCHEMA})
    (Path(directory)/'schema.json').write_text(json.dumps(schema, indent=2)+'\n')
    return schema
