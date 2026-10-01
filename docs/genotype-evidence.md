# Genotype Evidence

`pancgi.genotype_evidence.jsonl.gz` contains one JSON record per locus and
HAL genome. Each record has `locus_id`, `hal_genome`, `genotype` and `evidence`.
Result schema 1.3 includes its machine-readable JSON Schema in `schema.json`.
The genotype is the same `0`, `1` or `NA` value as the locus matrix.

## Observed Members

A record with genotype `1` has:

```json
{"callability":"not_evaluated_observed","reason":"observed_member"}
```

Membership establishes presence directly. Graph callability is not evaluated
for that locus-haplotype. Member and allele tables identify the observed copies.

## Graph Callability

Records with genotype `0` or `NA` contain:

| Field | Meaning |
| --- | --- |
| `callability` | `callable` for genotype `0`; `unresolved` for `NA` |
| `reason` | `at_least_one_certified_source` or `no_certified_source` |
| `source_model_n` | Number of distinct evaluated sources |
| `certified_source_n` | Number of sources supporting callability |
| `supporting_sources` | Certified source objects, in source-decision order |
| `source_decisions` | One decision per source, including unresolved sources |

Each source decision records `source`, `callability`, `reason`, `best_score`,
`runner_score`, `placements` and `placement`. Reasons are
`dominant_ordered_placement`, `competing_placements` and `no_supported_placement`.
Scores and the number of competing placement groups are retained from genotyping.
`placement` is null when no supported placement exists.

## Source Identity

A member source has `type=member`, `member_id`, `hal_genome`, `hal_sequence`,
`start0` and `end0`. It links to the same interval in `pancgi.members.tsv.gz`
and is either a locus anchor or an allele representative. These coordinates
describe the source CGI; the genotyping envelope may additionally span
intersecting insertions.

A novel-locus reference context has `type=primary_reference_context`,
`locus_id`, `hal_genome`, `hal_sequence`, `start0` and `end0`. It uses the
primary-reference interval recorded for that novel locus in `pancgi.loci.tsv.gz`.
Both source types use 0-based, half-open intervals and literal HAL identities.

## Target Placement

A placement contains the target `hal_genome`, literal `hal_sequence`,
`direction` (`1` or `-1`), `estimated_start0` and `estimated_end0`.
The coordinates are the boundary estimates from the best anchor chains on
the target's forward coordinate system. They are retained as finite numbers;
equal bounds represent a point placement. They describe graph placement,
not an observed CGI interval.

`run_manifest.json` includes the mapped sequence inventory. The validator
checks source identities and coordinates, source-set completeness, target
identities, score/decision consistency, certified-source counts and agreement
with genotype matrices and observed members.
