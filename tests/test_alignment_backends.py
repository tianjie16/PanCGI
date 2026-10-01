from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cpgi_nr_prod as prod
import wfa_longalign_wrapper as wfa


def main() -> None:
    parasail_identity = prod.parasail_semiglobal_identity("ACGTACGT", "ACGTACGT")
    wfa_result = wfa.alignment_identity("ACGTACGT", "ACGTACGT")
    assert parasail_identity == 1.0
    assert wfa_result["identity"] == 1.0
    assert wfa_result["alignment_length"] == 8
    print("OK: PanCGI alignment backend tests passed")


if __name__ == "__main__":
    main()
