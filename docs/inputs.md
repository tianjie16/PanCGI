# Input Formats


- Sequence-bearing GFA 1.x with complete `P` or `W` paths.
- HAL containing every selected haplotype and reference sequence.
- `genomes.tsv`, exactly four tab-separated columns in this order: `hal_genome`, `role`, `cpgi_bed`, `sv_tsv`.
- `paths.tsv`, exactly seven tab-separated columns in this order: `gfa_record_type`, `gfa_sample`, `gfa_haplotype`, `gfa_sequence`, `gfa_path_name`, `hal_genome`, `hal_sequence`.

Headers are literal, case-sensitive, and required. Extra, renamed, reordered or missing columns are errors. Paths in `genomes.tsv` are absolute or relative to that manifest. A genome is one haplotype, not one diploid individual. There is no fixed sample count, species, chromosome name or haplotype naming convention.

The genome names `locus_id` and `allele_id` are reserved for matrix identifier columns.

Select exactly one `primary_reference`, optionally one `comparison_reference`, and at least one `sample`. The comparison reference contributes its CGI sequences but does not define the primary coordinate system or reference-locus anchors. The primary reference must also occur in the GFA.

Use HAL's literal genome and sequence names in manifests and input BED/SV sequence columns. A GFA path name may be entirely different. For each W row, set `gfa_record_type=W`, copy the literal sample, haplotype and sequence fields from the GFA, and leave `gfa_path_name` empty. For each P row, set `gfa_record_type=P`, supply the complete original path name and leave the three W fields empty. These rules apply per row, including mixed W/P files. Do not split P names or infer parent-of-origin labels. `hal_sequence` means the contig name, not its nucleotide sequence. Users confirm that both sides identify the same haplotype and contig in the same sequence version and orientation. Every selected HAL sequence requires one complete path unless explicitly declared as an unused HAL-only sequence as described below. A path cannot be assigned to two genomes or sequences.

A HAL-only contig with no GFA path and no CGI/SV input may be declared in an optional two-column TSV: `hal_genome`, `hal_sequence`. Supply it with `--hal-only-exclusions FILE` to `run` or `validate-inputs`, and omit these rows from `paths.tsv`. Eligibility requires confirming that no corresponding GFA path exists. Exclusions must be unique, selected, unmapped sequences with no CGI, assembly-side SV or reference-side SV records. Each selected genome must retain a mapped path. All analysis-scope GFA paths remain required, including those without CGI calls, because they contribute genotyping context. The validated exclusions table, input lock and public `qc.json` record the declared contigs, lengths, reason and verified zero input counts.

The inventory records path names and declared coordinates. GFA and HAL inputs must represent the same assemblies, sequence versions and orientations. GFA1 `W` fragments require explicit, contiguous coordinates from zero. GFA1 `P` paths require zero overlaps, or `*` with unambiguous matching zero-overlap links. Nonzero overlaps, gaps, jumps, topology-only graphs, GFA2 ordered groups and missing segment sequences are unsupported. During unfolding, PanCGI checks referenced nodes, computed spans and HAL sequence lengths. Whole-genome sequence equality is an input-preparation requirement; these checks cover identities, structure and lengths rather than comparing all GFA and HAL bases.

Public tables and genotype evidence use literal HAL names and run-scoped member, locus and allele IDs. The work directory retains the input identity map. Cross-run comparisons use explicit genome/sequence/interval identities and complete member sets; adding genomes or reclustering can change the IDs.

## CGI BED10

Each genome has its own headerless, tab-delimited BED10:

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

This is the ten-column caller interval format, not a UCSC browser table with an added leading `bin` field. CGI calls are supplied by the user and are not regenerated. A valid empty BED explicitly represents no CGI calls for that genome; a missing BED is an error. Duplicate intervals, inconsistent lengths, invalid counts and out-of-bounds coordinates are errors. Percentages must agree with counts within half of the last reported decimal unit, at most 0.5 percentage points, plus 1e-9 numerical tolerance. The observed/expected ratio is supplied by the caller; it is not recomputed from GC count alone. Complete unmasked assembly sequences are extracted from HAL for allele comparison.

## SV TSV

Each genome has a TSV with exactly these nine columns, in this order:

```text
ref_contig	ref_pos1	ref_end1	sv_id	sv_type	asm_contig	asm_start0	asm_end0	sv_length
```

Supply decomposed `INS` and `DEL` events using the coordinate contract below. Other event types are rejected. VCF conversion and caller-specific coordinate conversion belong in input preparation.

`ref_contig` is a primary-reference HAL sequence. `asm_contig` is a HAL sequence of the current genome. SV IDs must be nonempty and unique within that genome. A header-only TSV is an explicit zero-event callset; absence of the file is not equivalent to zero events.

`sv_length` is the mandatory original VCF SVLEN: one positive integer for INS or negative integer for DEL. INS selection uses `abs(sv_length)`; interval overlap, overlap fractions and genotype flank coordinates use the supplied assembly coordinates. SVLEN and assembly span can differ. The original signed value is retained in `pancgi.member_sv.tsv.gz`.

SV overlaps are annotated for sample haplotypes. Both reference roles require header-only SV tables.

| Event | Reference coordinates | Assembly coordinates |
| --- | --- | --- |
| INS | `ref_pos1 = ref_end1 = POS`, the 1-based left padding base | `[asm_start0, asm_end0)` is the inserted sequence, excludes padding, and has positive length |
| DEL | `ref_pos1 = POS`, `ref_end1 = END > POS`; deleted reference bases are `[POS, END)` in 0-based coordinates | `asm_start0 = asm_end0 = q`, the zero-length deletion junction |

Coordinates use a left padding base, independently of left normalization. Events requiring right padding at a sequence start are outside this contract.

DEL-junction support requires `CGI_start0 < q < CGI_end0`. Boundary-only contact is excluded. A point junction is not inflated to a one-base interval. Its assembly span and overlap are zero; the deleted reference span is `END-POS`. Neither quantity replaces the supplied SVLEN. Overlap divided by assembly SV span is undefined, not zero.

Assembly coordinates refer to the forward sequence stored in the selected HAL genome. Input preparation includes caller-specific coordinate conversion and explicit accession-to-HAL mapping. For example, a 1-based inclusive inserted interval `10987-11126` becomes `[10986,11126)`. A caller-defined junction encoded as `q+1-q+1` becomes `[q,q)` only when that interbase convention is documented by the caller. Coordinate provenance remains required because a syntactically valid one-base shift may pass range checks.

For INS, positive interval overlap is required. The mechanism priority is: `inside_ins`, then `partial_ins`, then `contains_ins`. These denote CGI contained within an insertion, crossing an insertion boundary, and spanning the complete insertion, respectively. SV overlap is annotation evidence, not a variant-alone causal claim.

Reference association checks the reference padding base of every intersecting INS. Among INS whose padding base overlaps a primary-reference CGI, the longest INS chooses the reference association. Ties retain the established descending key `(SV length, CGI overlap bp, insertion span, SV ID)`; multiple reference CGI hit by that one event retain the established reference ranking. The globally longest overlapping INS is recorded separately. A novel assignment requires no usable HAL reference-CGI match and no matched INS anchor. Reference contact is inclusive of the padding base; it does not mean the inserted bases exist in the reference CGI.
