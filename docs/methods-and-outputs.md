# Computation And Outputs


The pipeline performs path and HAL inventory, input validation, graph unfolding, CGI sequence extraction, HAL projection, feature construction, graph locus clustering, locus polishing, primary-reference association, sequence-allele clustering, multi-scale graph genotyping and output export. Selected path records are retained while building the node index. Progress records distinguish source scanning, index construction and completed paths.

### Parallel execution

`--threads N` controls feature construction and path unfolding. With one worker these stages run serially. With more workers, feature construction processes independent genomes and graph unfolding processes independent complete paths. Workers call the same per-genome and per-path scientific functions. Feature records and exclusions are merged in manifest order; the path inventory retains its sorted-ID order regardless of task completion order. The global shingle-frequency index is built once over all accepted features.

Graph source scanning, node-index construction and path spooling use a single pass. Path workers read the same completed index: dense node lengths use a shared read-only memory mapping, while sparse or arbitrary node identifiers use read-only SQLite connections. Node identifiers are assigned only by the index builder. Each path worker writes its own partial pathBED; final paths and the inventory are published only after all selected paths and the final source-state check pass. Independent path expansion is parallel; the source scan is serial.

At most N tasks are in flight in graph unfolding and feature construction; the worker count is also bounded by the number of input genomes or paths. Feature workers write private compressed shards and return counts and filenames. Shards are merged as concatenated gzip members, so readers must consume the complete gzip stream.

Feature progress at `work/features.jsonl.gz.progress.jsonl` records completed genomes and their input, accepted and excluded counts. Graph progress at `work/pathbed/progress.jsonl` records completed paths; `work/pathbed/index/path_progress/` retains per-path events. Errors stop the stage, cancel queued work and retain diagnostic logs.

Processes have private working memory and caches even when their node-length mapping is shared. Set `--threads` within the scheduler's CPU and memory budget. `--hal-threads` controls HAL concurrency separately. Scientific parameters are independent of worker counts.

GFA and HAL provenance uses file size and timestamps, recorded as `file_metadata`, with source-state checks during execution. Content hashes are computed for smaller manifests, CGI/SV tables, code and result files using SHA-256. Large GFA/HAL inputs are not hashed in a separate pass.

HAL projection uses `halLiftover --bedType 4 --outPSLWithName SOURCE.hal SOURCE_GENOME CGI.bed4 PRIMARY_GENOME output.psl`. The source genome is the literal manifest value. The four-column query uses HAL sequence names and stable CGI identities. Projection coverage defaults to 0.5 and identity to 0.0; multi-hit records remain explicit. Graph feature coverage is 0.95. Locus construction caches 1,000 bp and 32 steps of feature context; graph-genotyping context is defined separately below.

Locus construction uses graph similarity, medoid polishing and non-reference site grouping with a 100 bp grouping window. Alleles use identity threshold 0.80, minimum length ratio 0.80, all-pairs clique clustering and sequence-medoid representatives. Ordinary sequence comparison uses Parasail semiglobal `sg_stats_scan_32` with match 2, mismatch -3, gap-open 5 and gap-extension 2. Alignment scores use 32-bit precision. Sequences of at least 100,000 bp use WFA. Alignment failure or saturation stops execution. The pinned pywfa interface requires successful status and a complete CIGAR, where I consumes the first input and D the second.

### Mechanism annotation

Sample members receive one of `INS-internal`, `INS-junction`, `INS-spanning`, `DEL-junction`, `INS+DEL`, or `non-SV`. A qualified INS plus a strictly internal DEL junction gives INS+DEL; among INS relationships internal, junction, spanning priority applies. Allele classification uses its representative member, not the union of events across all members. Per-member events remain available separately. A locus exports counts of allele mechanisms rather than an arbitrary single label.

`non-SV` means no qualified event in the supplied INS/DEL annotation; it makes no causal or small-variant attribution. Reference roles use header-only SV inputs and have `not_applicable_reference_role` status with empty mechanism fields. Reference/Novel locus type and sample-set membership are independent of mechanism.

### Multi-scale graph genotyping

Observed CGI members define state `1`. An unobserved CGI is `0` when graph homology establishes a callable locus, or `NA` when graph homology is insufficient or ambiguous.

For primary-reference loci, sources are the primary CGI anchor and allele representatives with duplicate member IDs removed. Novel loci use allele representatives plus the primary-reference interval context. Insertion-associated source envelopes cover the CGI and intersecting INS intervals. Full path intervals are clipped exactly within nodes; long nodes and partial-node offsets are retained. Each source retains up to 1 Mb of available context on either side without a step-count cap.

| Parameter | Value |
| --- | --- |
| Anchor scales | 100, 1,000, 10,000 bp |
| Minimum supporting bp | 50 per side |
| Minimum support fraction | 0.05 per side |
| Placement dominance margin | 0.02 |
| Complete-link placement tolerance | 1,000 bp |
| Source/target span tolerance | 10,000 bp |
| Gap penalty / repeat exponent | 0 / 0 |
| Occurrence truncation / source competition veto | None / none |

All exact target graph occurrences are retained. Within each source, orientation-consistent left/right anchor chains are scored at all three scales and combined before placement competition. Sources are evaluated independently. One certified source is sufficient; ambiguity in an additional source does not revoke its support.

Within one locus-haplotype, any observed allele keeps `1` and establishes locus `1`. Other unobserved alleles become `0`. Without observed CGI, any certified source establishes locus `0` and all its alleles `0`. Otherwise locus and unobserved alleles remain `NA`. Input and execution errors stop the run.

## Final outputs

`results/` contains:

- `pancgi.locus_genotypes.tsv.gz`, `pancgi.allele_genotypes.tsv.gz`: identifiers followed by literal HAL genome columns, in manifest order; states are `0`, `1`, `NA`.
- `pancgi.loci.tsv.gz`, `pancgi.alleles.tsv.gz`: primary coordinates or representative member identity, CGI properties and sample-only genotype statistics.
- `pancgi.members.tsv.gz`: one row per observed called CGI copy, linked to its allele/locus and original genome, sequence, interval, annotation and reference-association evidence.
- `pancgi.member_sv.tsv.gz`: one row per intersecting member/event relationship, with both explicit coordinate systems and event IDs.
- `pancgi.alleles.fa.gz`: actual representative CGI sequences, keyed by allele ID.
- `genomes.tsv`: declared genome roles; `completed.json`: counts and output checksums.
- `schema.json`: table columns, types, missing values, units, key relationships and the genotype-evidence JSON Schema.
- `run_manifest.json`: software profile, named scientific and execution parameters, dependency versions, source and logical input hashes, genotype configuration and literal sequence inventory. Machine paths and executable command strings remain in work logs.
- `qc.json`: input/accepted/excluded member counts and exclusion/genotype reasons; `pancgi.genotype_evidence.jsonl.gz`: locus-haplotype decisions and source evidence.
- `validation.json`: relational, evidence and genotype validation before completion. Verify a delivered directory using `pancgi validate-results RESULTS`.

The [genotype-evidence specification](genotype-evidence.md) defines public source identities, target placements and the distinction between observed membership and graph callability.

Allele IDs have the form `{hal_genome}_{hal_sequence}:{start0}:{end0}`, for example `HG002_hap1_chr1:950:1350`. The names and 0-based half-open interval belong to the representative CGI's source assembly, not its primary-reference projection. Use the complete literal HAL genome name, which identifies one haplotype or reference; no additional haplotype field is inferred. Each name is percent-encoded from UTF-8, leaving only ASCII letters, digits and `-._~` unescaped. Thus `:` becomes `%3A`, `%` becomes `%25`, and a space becomes `%20`; underscores remain unchanged. Independent genome, sequence and coordinate fields are available through `representative_member_id` in the member table. Do not split the ID on underscores to recover them. Distinct alleles producing the same ID stop export rather than receiving automatic suffixes. The same ID is used in the allele catalogue, allele genotype matrix, member and member-SV tables, and FASTA headers. IDs can change when the representative changes; cross-run comparisons require explicit identities and complete member sets.

Member primary coordinates denote the resolved HAL projection. `primary_strand` is the relative orientation: equal two-sign PSL orientations yield `+`, opposite orientations yield `-`. Locus `coordinate_kind` distinguishes a primary-reference CGI interval from a novel-site cluster summarizing insertion anchors. `coordinate_source` records the evidence source. Unavailable coordinates and undefined frequencies are empty.

`ac` counts sample haplotypes with genotype 1, not member copies. `n_called` counts sample 0/1 calls; `af = ac/n_called`, `maf = min(af,1-af)`, `call_rate = n_called/sample_haplotype_n`. An all-missing frequency is blank, not zero. `in_sample_set=1` means at least one observed positive sample haplotype. Reference-only loci/alleles remain available for catalogue inspection with this flag zero; reference haplotypes are never included in sample AF/MAF denominators. PanCGI does not apply a population-analysis MAF/call-rate filter itself.

The work directory contains intermediate files and execution logs.
