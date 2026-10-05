import pandas as pd

INTEGERS = set('''
locus_insert_anchor_member_n rep_insert_anchor_flag major_allele
has_primary_ref has_comparison_ref allele_member_n_total allele_member_n_asm allele_member_n_ref
allele_hap_n_total allele_hap_n_asm allele_hap_n_ref locus_member_n_total locus_member_n_asm locus_member_n_ref
locus_hap_n_total locus_hap_n_asm locus_hap_n_ref n_total_assemblies ref_only_flag rep_is_assembly
rep_seq_len_bp rep_cpgi_len_bp rep_len rep_cpg_n rep_gc_n rep_asm_mid0
locus_nonref_site_start1 locus_nonref_site_end1 rep_primary_start1 rep_primary_mid1 rep_primary_end1 rep_primary_anchor1
rep_sv_ins_n rep_sv_longest_site1 rep_sv_longest_start1 rep_sv_longest_end1 rep_sv_longest_len
allele_callable_asm_n allele_callable_total_n locus_callable_asm_n locus_callable_total_n
allele_gt_1_asm_n allele_gt_0_asm_n allele_gt_na_asm_n allele_gt_1_total_n allele_gt_0_total_n allele_gt_na_total_n
locus_gt_1_asm_n locus_gt_0_asm_n locus_gt_na_asm_n locus_gt_1_total_n locus_gt_0_total_n locus_gt_na_total_n
merged_cpgi_n locus_cpgi_n contains_primary_member contains_comparison_member
rep_primary_primary_start1 rep_primary_primary_end1 rep_primary_primary_site1
rep_hal_multimap_n rep_hal_primary_start0 rep_hal_primary_end0 rep_hal_primary_start1 rep_hal_primary_end1 rep_sv_all_n
'''.split())
NUMBERS = set('''
allele_asm_freq allele_freq rep_pct_gc rep_oe
locus_primary_start_median1 locus_primary_mid_median1 locus_primary_end_median1
rep_hal_query_cov rep_hal_identity
'''.split())


def merged_dataframe(rows, columns, genotype_columns):
    df = pd.DataFrame(rows, columns=columns)
    if set(genotype_columns) & (INTEGERS | NUMBERS):
        raise ValueError('Internal genotype labels collide with merged field names')
    for column in columns:
        if column in INTEGERS | NUMBERS:
            values = df[column].mask(df[column].eq(''), pd.NA)
            try:
                df[column] = pd.to_numeric(values, errors='raise').astype('Int64' if column in INTEGERS else 'Float64')
            except (ValueError, TypeError) as error:
                raise ValueError(f'Invalid merged numeric field: {column}') from error
        else:
            df[column] = df[column].astype('string')
    return df
