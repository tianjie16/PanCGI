# Computational resources

Use `--threads` to set the number of analysis workers and `--hal-threads` to
set the number of concurrent HAL workers. Allocate CPU, memory and wall time
through the scheduler or container. In Docker mode, `--hal-cpus` and
`--hal-memory` set the limits for each HAL worker.

Provide an existing writable directory with `--scratch`. Set `TMPDIR` and
`SQLITE_TMPDIR` to the intended temporary-file directory. The output directory
also holds intermediate databases; provision storage for both locations.

Feature indexes on shared storage can be read from another compute node.
The source files must remain accessible and their content must match the index.
Readers verify content when the recorded filesystem state differs and reject
source or index changes during use.

CPU, memory and storage requirements depend on the number of assemblies,
graph size and CGI content. Set worker counts to match the allocated resources.
Container memory accounting includes file cache as well as process memory.
Cache budgets within the application are not whole-job memory limits.

Use a separate output directory for each analysis and keep its inputs unchanged
throughout the run. Final tables and representative sequences are written to
`results/` within the output directory.
