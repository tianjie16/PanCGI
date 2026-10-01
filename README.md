# PanCGI

PanCGI builds non-redundant CpG-island loci, sequence alleles and haplotype
genotypes from CGI intervals, a sequence-bearing pangenome graph and a HAL
alignment. It reports observed membership and graph-supported presence,
absence or insufficient evidence as `1`, `0` and `NA`.

<p align="center">
  <img src="assets/PanCGI.png" alt="PanCGI locus construction, allele clustering and genotyping" width="1000">
</p>

## Installation

Linux, Python 3.12 and Bash are required. Install this source directory:

```bash
python -m pip install .
pancgi --version
pancgi --help
```

HAL executables are separate dependencies. See [installation](docs/installation.md)
for native and Docker-backed HAL execution and a runnable small example.

## Run

```bash
pancgi prepare --gfa graph.gfa --hal alignment.hal --hal-runtime native --out-dir prepared
pancgi run --prepared prepared --gfa graph.gfa --hal alignment.hal \
  --genomes genomes.tsv --paths paths.tsv --threads 4 --hal-threads 1 \
  --hal-runtime native --out-dir analysis
pancgi validate-results analysis/results
```

Review and complete the generated genome and path tables before running.
Map original GFA W/P identities to literal HAL names using the generated
tables. Use a new output directory for each run.

## Inputs And Results

- [Input formats](docs/inputs.md): genome roles, seven-column W/P mapping,
  CGI BED10, signed INS/DEL SVLEN and coordinate conventions.
- [Outputs](docs/outputs.md): result files, identifiers, coordinates,
  CGI states and sample frequencies.
- [Resource management](docs/resources.md): CPU, memory, scratch and completion.
- [Dependencies](DEPENDENCIES.md): external software and upstream projects.

CGI intervals and haplotype identities are supplied as inputs. Required input,
alignment and source-state checks must pass for a run to complete.

## Tests

```bash
bash tests/run_all.sh
```

The test suite uses synthetic data.
For contributions and issue reports, see [CONTRIBUTING.md](CONTRIBUTING.md).

## HPRC2 Resources

The [HPRC2 PanCGI resource](https://github.com/tianjie16/PanCGI/releases/tag/resource-v1.0.0)
contains 47,039 CGI alleles at 33,463 loci, assembly-specific membership,
representative sequences, CGI state matrices, evolutionary trajectories and
masked and unmasked CGI calls for 462 HPRC haplotypes.
See the [resource guide](resources/README.md) for file definitions,
coordinates and download checksums.

## License

Copyright 2026 Tianjie Liu. PanCGI is licensed under [Apache-2.0](LICENSE).
Dependencies retain their respective licenses; see [Dependencies](DEPENDENCIES.md).

## Project

Source and issue tracker: [tianjie16/PanCGI](https://github.com/tianjie16/PanCGI).
When citing the software, include its version and repository URL.

## Human Pangenome Reference Consortium

<p align="left">
  <a href="https://humanpangenome.org/">
    <img src="assets/HPRC_white_logo.png" alt="Human Pangenome Reference Consortium" width="440">
  </a>
</p>

[HPRC resources](https://github.com/human-pangenomics) ·
[Study code and resources](https://github.com/tianjie16/HPRC2_DNA_methylome_variation) ·
[Methylation framework](https://github.com/tianjie16/HPRC2_methylation_framework)
