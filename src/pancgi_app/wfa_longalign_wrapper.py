from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import List, Tuple

try:
    from pywfa import WavefrontAligner
except Exception:
    WavefrontAligner = None


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def read_single_fasta(path: str) -> tuple[str, str]:
    head = ""
    seq_parts: list[str] = []
    with open(path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if head:
                    raise ValueError(f"Expected single-record FASTA, got another header in {path}")
                head = line[1:]
            else:
                seq_parts.append(line)
    if not head:
        raise ValueError(f'Missing FASTA header: {path}')
    seq = "".join(seq_parts)
    if not seq:
        raise ValueError(f"Empty FASTA sequence: {path}")
    return head, seq


def parse_cigar(cigar: str) -> List[Tuple[str, int]]:
    if not cigar or ''.join(m.group(0) for m in _CIGAR_RE.finditer(cigar)) != cigar:
        raise ValueError('Invalid WFA CIGAR')
    out: List[Tuple[str, int]] = []
    for n, op in _CIGAR_RE.findall(cigar or ""):
        if int(n) <= 0 or op not in {'M', 'I', 'D', '=', 'X'}:
            raise ValueError('Unsupported WFA global CIGAR operation')
        out.append((op, int(n)))
    return out


def replay_cigar_identity(seq1: str, seq2: str, cigar: str) -> Tuple[int, int, dict[str, int]]:
    ops = parse_cigar(cigar)
    i = j = matches = aln_len = 0
    hist: dict[str, int] = {}
    for op, n in ops:
        hist[op] = hist.get(op, 0) + n
        if op in {'=', 'X', 'M'}:
            if i + n > len(seq1) or j + n > len(seq2):
                raise RuntimeError('WFA CIGAR exceeds input sequence')
            actual = sum(a == b for a, b in zip(seq1[i:i+n], seq2[j:j+n]))
            if (op == '=' and actual != n) or (op == 'X' and actual != 0):
                raise RuntimeError('WFA CIGAR contradicts sequence bases')
            matches += actual
            i += n
            j += n
        elif op == 'I':
            i += n
        elif op == 'D':
            j += n
        aln_len += n
        if i > len(seq1) or j > len(seq2):
            raise RuntimeError('WFA CIGAR exceeds input sequence')
    if i != len(seq1) or j != len(seq2):
        raise RuntimeError('WFA CIGAR does not consume both complete input sequences')
    return matches, aln_len, hist


def alignment_identity(seq1: str, seq2: str) -> dict[str, object]:
    if WavefrontAligner is None:
        raise RuntimeError("pywfa is not installed; install pywfa to use the WFA very-long wrapper")
    aligner = WavefrontAligner()
    aligner(seq1, seq2)
    if aligner.status != 0:
        raise RuntimeError(f'WFA alignment failed with status {aligner.status}')
    cigar = aligner.cigarstring
    matches, aln_len, op_hist = replay_cigar_identity(seq1, seq2, cigar)
    if aln_len <= 0:
        raise RuntimeError('Empty WFA alignment')
    ident = float(matches) / float(aln_len)
    return {
        "identity": ident,
        "matches": matches,
        "alignment_length": aln_len,
        "score": getattr(aligner, "score", None),
        "status": getattr(aligner, "status", None),
        "cigar": cigar,
        "op_hist": op_hist,
    }


def main() -> None:
    p = argparse.ArgumentParser(description="Minimal WFA2-lib wrapper for cpgi very-long external backend")
    p.add_argument("--seq1", required=True)
    p.add_argument("--seq2", required=True)
    p.add_argument("--log", default="")
    args = p.parse_args()

    h1, s1 = read_single_fasta(args.seq1)
    h2, s2 = read_single_fasta(args.seq2)
    out = alignment_identity(s1, s2)
    if args.log:
        with open(args.log, "at") as fh:
            fh.write(
                "\t".join(
                    [
                        h1,
                        h2,
                        str(len(s1)),
                        str(len(s2)),
                        str(out.get("identity")),
                        str(out.get("alignment_length")),
                    ]
                )
                + "\n"
            )
    print(json.dumps(out, ensure_ascii=False))


if __name__ == "__main__":
    main()
