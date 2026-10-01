# Dependencies

Python packages are installed with PanCGI. HAL executables are supplied
separately, either natively or through the Cactus container. PanCGI is licensed
under Apache-2.0; dependencies retain their respective upstream licenses.

| Dependency | Role | Upstream |
| --- | --- | --- |
| Python 3.12 | Runtime | https://www.python.org/ |
| NumPy 2.4.2 | Numerical operations | https://numpy.org/ |
| pandas 3.0.1 | Tabular data | https://pandas.pydata.org/ |
| Apache Arrow / PyArrow 23.0.1 | Parquet output | https://arrow.apache.org/ |
| Biopython 1.86 | Sequence utilities and optional MSA output | https://biopython.org/ |
| Parasail 1.3.4 | Sequence alignment | https://github.com/jeffdaily/parasail |
| pywfa 0.5.1 | Long-sequence alignment | https://github.com/kcleal/pywfa |
| HAL | Genome inventory, extraction and projection | https://github.com/ComparativeGenomicsToolkit/hal |
| Cactus container | Optional HAL runtime | https://github.com/ComparativeGenomicsToolkit/cactus |
