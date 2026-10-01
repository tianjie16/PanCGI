# Installation

Use Linux with Python 3.12 and Bash. Python requirements are pinned in
`pyproject.toml`; `environment.yml` supplies a Conda environment definition.
Some dependencies may need a C/C++ compiler when a binary wheel is unavailable.

```bash
python -m pip install .
python -m pip check
pancgi --version
```

The wheel and source distribution contain the application, command scripts,
input documentation and synthetic example. Development tests are included
in the source distribution, not the installed application.

## HAL

Provide `halStats`, `hal2fasta` and `halLiftover` on PATH or pass their executable
paths explicitly. HAL is not bundled in the Python wheel. Native execution
uses the containing scheduler or container's resource limits. Select it
explicitly with `--hal-runtime native` for both `prepare` and `run`.

Alternatively, install Docker on the host and pull the public Cactus image:

```bash
docker pull quay.io/comparative-genomics-toolkit/cactus@sha256:646480e4b4870d57fd3a32c2724eea596079de59e8082d0067b9acfbd8cd5c5d
```

Use `--hal-runtime docker`, `--docker-image` with that image reference, and
`--hal-stats /home/cactus/bin/halStats`,
`--hal2fasta /home/cactus/bin/hal2fasta`,
`--hal-liftover /home/cactus/bin/halLiftover` as appropriate for the command.
Docker mode requires access to a running Docker daemon.

## Synthetic Example

With native HAL executables available, run from the source directory:

```bash
pancgi prepare --gfa examples/tiny/tiny.gfa --hal examples/tiny/tiny.hal \
  --hal-runtime native --out-dir example_prepared
pancgi run --prepared example_prepared \
  --gfa examples/tiny/tiny.gfa --hal examples/tiny/tiny.hal \
  --genomes examples/tiny/genomes.tsv --paths examples/tiny/paths.tsv \
  --hal-runtime native --threads 2 --hal-threads 1 --out-dir example_run
pancgi validate-results example_run/results
python tests/check_example.py example_run/results examples/tiny/expected.json
```

In Docker HAL mode, replace the native option with the runtime/image/executable options above for both
`prepare` and `run`. Every command exposes its accepted options with `--help`.
Use a fresh directory for each run; do not reuse partial results.

## Build Distributions

```bash
python -m pip install build
python -m build
```

Install the generated wheel into a clean environment with its dependencies.
The package version is reported by `pancgi --version` and recorded in outputs.
