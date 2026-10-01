# Installation

Use Linux with Python 3.12 and Bash. Python requirements are pinned in
`pyproject.toml`; `environment.yml` supplies a Conda environment definition.
Some dependencies may need a C/C++ compiler when a binary wheel is unavailable.

```bash
git clone https://github.com/tianjie16/PanCGI.git
cd PanCGI
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install .
python -m pip check
pancgi --version
```

The wheel and source distribution contain the application, command scripts,
input documentation and synthetic example. The source distribution also
includes development tests.

## Tested environment and installation time

The bundled example was tested on Ubuntu 22.04.5 LTS (x86-64, WSL2), Python
3.12.13 and HAL 2.2, using the dependency versions listed in the package.
No GPU or other specialized hardware is required. Two CPU cores and 2 GB RAM
are sufficient for the bundled example. Larger datasets require additional
memory and storage, as described in [Computational resources](resources.md).

Installing the Python wheel takes less than one minute when Python dependencies
are already available. First-time dependency downloads and HAL installation
add to the setup time.

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

## Synthetic example

With native HAL executables available, run from the source directory:

```bash
pancgi prepare --gfa examples/tiny/tiny.gfa --hal examples/tiny/tiny.hal \
  --hal-runtime native --out-dir example_prepared
pancgi run --prepared example_prepared \
  --gfa examples/tiny/tiny.gfa --hal examples/tiny/tiny.hal \
  --genomes examples/tiny/genomes.tsv --paths examples/tiny/paths.tsv \
  --hal-runtime native --threads 2 --hal-threads 1 --out-dir example_run
pancgi validate-results example_run/results
```

In Docker HAL mode, replace the native option with the runtime/image/executable options above for both
`prepare` and `run`. Every command exposes its accepted options with `--help`.
Use a fresh output directory for each run.

### Expected output and runtime

The example produces two loci, three alleles and four CGI members. Result tables
and representative sequences are written to `example_run/results/`; the main
tables are `pancgi.loci.tsv.gz`, `pancgi.alleles.tsv.gz` and
`pancgi.members.tsv.gz`. Expected counts and identifiers are supplied in
`examples/tiny/expected.json`.

Allow about one minute for the example. Preparation, analysis and output checks
took approximately 15 seconds on an Intel Core i9-14900K with two CPU cores
allocated, using the tested environment above.

## Development tests

Run the test suite from the source directory:

```bash
bash tests/run_all.sh
```

After running the bundled example, compare its results with the expected values:

```bash
python tests/check_example.py example_run/results examples/tiny/expected.json
```

## Build distributions

```bash
python -m pip install build
python -m build
```

Install the generated wheel into a clean environment with its dependencies.
The package version is reported by `pancgi --version` and recorded in outputs.
