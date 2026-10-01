# Resource Management

`--threads` limits concurrent computation workers; `--hal-threads` independently
limits HAL workers. Set CPU, memory and wall-time limits separately in the
scheduler or containing container.

Provide an existing writable `--scratch` directory. Set `TMPDIR` and
`SQLITE_TMPDIR` to that directory for temporary-file placement. The anchor
stage also stages SQLite tables and outputs beneath its result directory,
so place both the output and scratch trees on adequately provisioned storage.
Budget capacity for both trees, including SQLite sorting and temporary files.

Anchor metadata, assignments and grouping are disk-backed. Prepared features,
weights and metadata have bounded caches. Weighting uses the complete background.
A record or indivisible group exceeding its configured bound stops the stage.
Total memory also includes Python overhead, worker-private memory and operating-system
file cache; cache counters measure their named cache only.

Allele comparison can require pairwise work within a locus. Graph-genotyping placement work
depends on graph repetition and locus complexity. Size CPU, RAM and scratch
for the actual cohort and locus distribution.

Stage logs and receipts are written under `work/stages`. Anchor records its
execution in `work/results_internal/locus_anchored.tsv.gz.anchor.json` and
phase events in the adjacent `.progress.jsonl` file. Counters describe their
named phase. Parent-process cache statistics cover that process only.

`results/completed.json` is written after final result validation. Inputs must
remain unchanged during execution. Errors return a nonzero exit status and
retain logs and intermediate outputs. Use a new output directory to restart.
