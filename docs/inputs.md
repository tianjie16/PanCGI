# Input formats

Provide the following files:

- A sequence-bearing GFA 1.x graph with complete `P` or `W` paths.
- The corresponding HAL file containing the selected haplotypes and references.
- `genomes.tsv`, with four tab-separated columns in this order: `hal_genome`, `role`, `cpgi_bed`, `sv_tsv`.
- `paths.tsv`, with seven tab-separated columns in this order: `gfa_record_type`, `gfa_sample`, `gfa_haplotype`, `gfa_sequence`, `gfa_path_name`, `hal_genome`, `hal_sequence`.

Use the exact column names and order, including a header. File paths in
`genomes.tsv` are absolute or relative to the directory containing that table.
Each genome represents one haplotype. The names `locus_id` and `allele_id` are
reserved for matrix identifier columns.

Select exactly one `primary_reference`, optionally one `comparison_reference`,
and at least one `sample`. The primary reference defines the reference coordinate
system and must occur in the GFA. The comparison reference contributes its CGI
sequences to the catalogue.

## Assembly and path mapping

Use the literal HAL genome and sequence names in the mapping tables and BED/SV
sequence columns. `hal_sequence` is the contig name.

| GFA record | Required GFA fields in `paths.tsv` | Fields left empty |
| --- | --- | --- |
| `W` | The original sample, haplotype and sequence fields | `gfa_path_name` |
| `P` | The complete original path name in `gfa_path_name` | `gfa_sample`, `gfa_haplotype`, `gfa_sequence` |

Apply these rules to each row, including graphs containing both record types.
Confirm that each graph path and HAL sequence identify the same haplotype,
contig, sequence version and orientation. Preserve complete path names. Each
selected HAL sequence requires one complete graph path, and each path can map
to only one HAL genome/sequence pair.

Include all paths for the selected assemblies, including sequences without CGI
calls. A selected HAL-only contig may be excluded when it has no corresponding
GFA path and no CGI, assembly-side SV or reference-side SV records. List these
contigs in a two-column TSV with headers `hal_genome` and `hal_sequence`, supply
it with `--hal-only-exclusions FILE`, and omit them from `paths.tsv`. Each
excluded contig must be listed once, and every selected genome must retain a
mapped path.

Use complete segment sequences. `W` fragments must have explicit, contiguous
coordinates beginning at zero. `P` paths must have zero overlaps, or `*` with
unambiguous matching zero-overlap links. The GFA and HAL must contain the same
assembly sequences and orientations, with matching path and sequence lengths.

## CGI BED10

Supply one headerless, tab-delimited BED10 file per genome:

| Column | Definition |
| --- | --- |
| 1 | HAL sequence name |
| 2 | CGI start, 0-based inclusive |
| 3 | CGI end, 0-based exclusive |
| 4 | Input CGI name |
| 5 | Length, equal to end minus start |
| 6 | Number of CpG dinucleotides |
| 7 | Number of C and G bases |
| 8 | Percentage CpG bases |
| 9 | Percentage C/G bases |
| 10 | Observed/expected CpG ratio |

Use these ten columns, beginning with the sequence name. PanCGI accepts existing
CGI calls. An empty BED represents a genome with no CGI calls. Intervals must
be unique, lie within their assembly sequences and have consistent lengths
and counts. Percentages must agree with counts within half of the last reported
decimal unit, up to 0.5 percentage points plus a numerical tolerance of 1e-9.
Supply the observed/expected ratio from the CGI caller and complete unmasked
assembly sequences in HAL.

## SV TSV

Supply one tab-delimited file per genome with exactly these nine columns:

```text
ref_contig	ref_pos1	ref_end1	sv_id	sv_type	asm_contig	asm_start0	asm_end0	sv_length
```

Use decomposed `INS` and `DEL` events. Convert caller-specific coordinates to
the following conventions before running PanCGI. `ref_contig` is a primary-reference
HAL sequence; `asm_contig` is a HAL sequence of the current genome. SV IDs must
be nonempty and unique within each genome.

Both reference genomes require header-only SV tables. For sample haplotypes,
a header-only table represents a zero-event callset.

| Event | Reference coordinates | Assembly coordinates |
| --- | --- | --- |
| INS | `ref_pos1 = ref_end1 = POS`, the 1-based left padding base | `[asm_start0, asm_end0)` contains the inserted sequence, excludes padding and has positive length |
| DEL | `ref_pos1 = POS`, `ref_end1 = END > POS`; deleted reference bases are `[POS, END)` in 0-based coordinates | `asm_start0 = asm_end0 = q`, the zero-length deletion junction |

Use left-padding coordinates. Supply the original signed VCF SVLEN as
`sv_length`: positive for INS and negative for DEL. SVLEN may differ from the
assembly interval length; retain each value in its corresponding field.

Assembly coordinates refer to the forward sequence stored in the selected HAL
genome. For example, a 1-based inclusive inserted interval `10987-11126` becomes
`[10986,11126)`. A deletion junction at interbase position `q` is represented as
`[q,q)`. Confirm the caller's coordinate convention before conversion.
