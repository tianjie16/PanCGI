set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH="$ROOT/src/pancgi_app${PYTHONPATH:+:$PYTHONPATH}"
"$PYTHON_BIN" "$ROOT/tests/test_graph_unfold.py"
"$PYTHON_BIN" "$ROOT/tests/test_alignment_backends.py"
"$PYTHON_BIN" "$ROOT/tests/test_hal_pipeline.py"
"$PYTHON_BIN" "$ROOT/tests/test_locus_polish_exact.py"
"$PYTHON_BIN" -m unittest discover -s "$ROOT/tests" -p 'test_*.py' -v
