import json
import math
import shlex
from pathlib import Path


COMMON = 'anchor_k max_mid_anchors'
SIMILARITY = 'max_anchor_bucket exact_ro exact_size exact_ctx ro_min ro_ctx szro_ro szro_size szro_ctx weak_seq_gate ambiguity_margin kmer'
WORKERS = 'threads chunksize maxtasksperchild'
STAGES = {
    'anchor-loci-primary': ('reference_association', COMMON + ' max_medoid_members nonref_site_window_bp hal_nonref_min_coverage ' + WORKERS + ' anchor_cache_bytes anchor_max_record_bytes anchor_max_pending anchor_group_bytes anchor_weight_cache_bytes'),
    'build-features-prod': ('feature_construction', 'threads min_graph_cov flank_bp flank_max_steps hal_min_coverage hal_min_identity hal_ambig_identity_delta hal_ambig_coverage_delta hal_ambig_aligned_bp_delta sv_contig_coordinate_base'),
    'cluster-alleles-prod': ('allele_clustering', COMMON + ' identity min_len_ratio very_long_threshold parasail_match parasail_mismatch parasail_gap_open parasail_gap_extend max_allele_medoid_members ' + WORKERS),
    'cluster-loci-prod': ('locus_clustering', COMMON + ' ' + SIMILARITY),
    'make-cpgi-fasta': ('sequence_extraction', 'threads'),
    'polish-loci': ('locus_polishing', COMMON + ' ' + SIMILARITY + ' max_medoid_members ' + WORKERS),
    'strict-genotype-prod': ('genotyping', 'threads'),
    'unfold-graph': ('graph_unfolding', 'threads'),
}
ENUMS = {
    'parallel_backend': ('process', 'thread'),
    'mp_start_method': ('fork', 'spawn', 'forkserver'),
    'seq_backend': ('parasail',),
    'very_long_backend': ('external',),
    'parasail_mode': ('sg',),
    'allele_rep_strategy': ('asm_seq_medoid',),
    'allele_cluster_mode': ('allpairs_clique',),
    'compression': ('gzip', 'none'),
}
PRIVATE_FIELDS = set(('cmd features catalog out excluded hal_psl_dir in_locus in_members out_locus out_members '
    'locus_catalog locus_members out_allele out_allele_members hal contigs log_dir docker_bin docker_image hal2fasta '
    'missing_report allele_catalog allele_members out_prefix gfa gfa_inventory out_dir expected_lengths '
    'anchor_resource_report feature_store_db feature_store_rebuild feature_store_batch_size feature_store_read_chunk '
    'dump_locus_dir dump_min_members dump_pairwise_matrix dump_msa_backend dump_msa_name dump_msa_min_members '
    'dump_msa_max_members msa_mafft_bin msa_match msa_mismatch msa_gap_open msa_gap_extend flush_every_loci '
    'locus_ids_file locus_ids locus_start_index locus_end_index locus_shard_count locus_shard_index locus_shard_mode '
    'write_locus_list locus_list_only feature_load_mode very_long_external_template').split())


def public_parameters(directory):
    result = {}
    for path in sorted(Path(directory).glob('*.json')):
        value = json.loads(path.read_text())
        command = value.get('cmd')
        if command not in STAGES or path.stem != command:
            raise ValueError('Unknown parameter record')
        stage, fields = STAGES[command]
        if set(value) - set(fields.split()) - set(ENUMS) - PRIVATE_FIELDS:
            raise ValueError('Unclassified stage parameter')
        if stage in result:
            raise ValueError('Duplicate parameter record')
        selected = {}
        for field in fields.split():
            number = value[field]
            if type(number) not in (int, float) or not math.isfinite(number):
                raise ValueError('Invalid numeric parameter: ' + field)
            selected[field] = number
        for field, choices in ENUMS.items():
            if field in value:
                if value[field] not in choices:
                    raise ValueError('Unsupported public parameter: ' + field)
                selected[field] = value[field]
        if command == 'cluster-alleles-prod':
            for field in ('locus_ids_file', 'locus_ids', 'locus_start_index', 'locus_end_index', 'locus_shard_count', 'locus_shard_index', 'locus_list_only'):
                if value[field] not in ('', None, 0, False):
                    raise ValueError('Public full-run export requires complete locus selection')
            tokens = shlex.split(value['very_long_external_template'])
            if len(tokens) != 8 or Path(tokens[1]).name != 'wfa_longalign_wrapper.py' or tokens[2:7] != ['--seq1', '{seq1}', '--seq2', '{seq2}', '--log']:
                raise ValueError('Unrecognized long-alignment command')
            selected['very_long_backend'] = 'pywfa'
            selected['locus_selection'] = 'all'
        result[stage] = selected
    if set(result) != {value[0] for value in STAGES.values()}:
        raise ValueError('Incomplete stage parameter records')
    return result
