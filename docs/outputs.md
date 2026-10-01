# Outputs

PanCGI writes its result files to `results/` within the analysis directory.

| File | Contents |
| --- | --- |
| `pancgi.loci.tsv.gz` | CGI loci, primary-reference coordinates and sample genotype statistics |
| `pancgi.alleles.tsv.gz` | CGI alleles, representative member identities, sequence properties and sample genotype statistics |
| `pancgi.members.tsv.gz` | Assembly-specific CGI calls linked to their allele, locus, source sequence and coordinates |
| `pancgi.member_sv.tsv.gz` | CGI member overlaps with supplied insertion and deletion annotations |
| `pancgi.locus_genotypes.tsv.gz` | Locus-by-haplotype CGI states |
| `pancgi.allele_genotypes.tsv.gz` | Allele-by-haplotype presence and absence |
| `pancgi.alleles.fa.gz` | Representative CGI sequences, with allele IDs as FASTA record names |
| `pancgi.genotype_evidence.jsonl.gz` | Supporting records for each locus-haplotype genotype |
| `genomes.tsv` | Genome names and their declared roles |
| `schema.json` | Field definitions, data types, units, missing values and table relationships |
| `run_manifest.json` | Software version, parameters, dependencies and input metadata |
| `qc.json` | Input, member and genotype summaries |
| `validation.json` | Result validation summary |
| `completed.json` | Result counts and file checksums |

## Table relationships

Join members to alleles by `allele_id` and alleles to loci by `locus_id`.
The `representative_member_id` links each allele to its representative in the member table.
Genotype matrix columns use the declared HAL genome names in input manifest order.

Allele IDs follow `{hal_genome}_{hal_sequence}:{start0}:{end0}`, for example
`HG002_hap1_chr1:950:1350`. The coordinates identify the representative CGI in
its source assembly. Genome and sequence names use UTF-8 percent encoding,
with ASCII letters, digits and `-._~` retained. Separate genome, sequence and
coordinate fields are provided in the member table; use these fields to recover
the components of an ID. IDs depend on the selected representative, so comparisons
between catalogues should use source identities and complete member sets.

## CGI states

| Value | Locus state | Allele state |
| --- | --- | --- |
| `1` | CGI present | Allele present |
| `0` | CGI absent at a recoverable homologous locus | Allele absent at a called locus |
| `NA` | Homologous locus unresolved | Allele state unresolved |

Other tables use empty fields for unavailable coordinates and undefined frequencies.

## Coordinates

Assembly CGI intervals are 0-based and half-open. Member primary coordinates
describe the resolved projection to the primary reference. `primary_strand`
records its relative orientation as `+` or `-`.

For loci, `coordinate_kind` identifies a primary-reference CGI interval or a
novel-site cluster summarizing insertion anchors. `coordinate_source` identifies
the source of those coordinates. SV coordinates follow the event-specific
definitions in [Input formats](inputs.md).

## Sample frequencies

- `ac`: sample haplotypes carrying the locus or allele.
- `n_called`: sample haplotypes with state `0` or `1`.
- `af`: `ac / n_called`.
- `maf`: `min(af, 1 - af)`.
- `call_rate`: `n_called / sample_haplotype_n`.
- `in_sample_set`: `1` when at least one sample haplotype carries the locus or allele.

Reference genomes are excluded from sample frequency denominators. Reference-only
entries have `in_sample_set=0`. Frequencies are empty when all sample genotypes
are missing. Apply frequency and call-rate thresholds for the intended downstream
analysis.

## Optional result check

PanCGI checks output consistency before marking a run complete. To check the
result files again after copying or transferring them, run:

```bash
pancgi validate-results analysis/results
```
