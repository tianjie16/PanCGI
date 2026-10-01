# Computational resources

Use `--threads` to set the number of analysis workers and `--hal-threads` to
set the number of concurrent HAL workers. Allocate CPU, memory and wall time
through the scheduler or container. In Docker mode, `--hal-cpus` and
`--hal-memory` set the limits for each HAL worker.

Provide an existing writable directory with `--scratch`. `TMPDIR` and
`SQLITE_TMPDIR` control temporary-file locations. Allow sufficient storage
for the output directory, intermediate files and temporary files.

CPU, memory and storage requirements depend on the number of assemblies,
graph size and CGI content. Set worker counts to match the allocated resources.

Use a separate output directory for each analysis and keep its inputs unchanged
throughout the run. Final tables and representative sequences are written to
`results/` within the output directory.
