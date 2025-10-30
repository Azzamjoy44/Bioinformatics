#!/usr/bin/env python3
"""
fasta_symbol_percentages.py

Read a FASTA file (DNA/RNA/protein) and display relative percentages
for each symbol in the sequence alphabet.

Usage:
  python fasta_symbol_percentages.py path/to/file.fasta
  python fasta_symbol_percentages.py file.fasta --alphabet auto    # default
  python fasta_symbol_percentages.py file.fasta --alphabet dna
  python fasta_symbol_percentages.py file.fasta --alphabet rna
  python fasta_symbol_percentages.py file.fasta --alphabet protein
  python fasta_symbol_percentages.py file.fasta --include-gaps --include-stop
"""

import argparse
import gzip
import io
import sys
from collections import Counter
from typing import Iterable, Tuple, Set

DNA_IUPAC = set("ACGTRYSWKMBDHVN")  # includes ambiguous + N
RNA_IUPAC = set("ACGURYSWKMBDHVN")  # includes ambiguous + N
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")  # 20 common AAs
PROTEIN_EXTENDED = PROTEIN_ALPHABET | set("BXZJUO")  # include common extras

def open_maybe_gzip(path: str):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")

def parse_fasta(handle: Iterable[str]) -> Iterable[str]:
    """Yield sequence strings (without headers) from a FASTA stream."""
    seq_chunks = []
    for line in handle:
        if not line:
            continue
        if line.startswith(">") or line.startswith(";"):  # headers/comments
            if seq_chunks:
                yield "".join(seq_chunks)
                seq_chunks = []
            continue
        seq_chunks.append(line.strip())
    if seq_chunks:
        yield "".join(seq_chunks)

def detect_alphabet(letters_present: Set[str]) -> Tuple[str, Set[str]]:
    """Guess alphabet based on observed letters (A–Z only)."""
    if not letters_present:
        return "auto", set()
    if letters_present <= DNA_IUPAC:
        return "dna", DNA_IUPAC
    if letters_present <= RNA_IUPAC:
        return "rna", RNA_IUPAC
    # If it looks like protein (has symbols outside DNA/RNA sets), assume protein
    return "protein", PROTEIN_EXTENDED

def tally_symbols(fasta_path: str,
                  alphabet: str = "auto",
                  include_gaps: bool = False,
                  include_stop: bool = False) -> Tuple[Counter, int, str, Set[str], Counter]:
    """
    Returns:
      counts_in_alphabet: Counter for symbols considered part of the alphabet
      total_in_alphabet: total symbols counted toward percentages
      resolved_alpha_name: 'dna'|'rna'|'protein'|'auto'
      resolved_alpha_set: the alphabet set used (if applicable)
      counts_other: Counter for other symbols (not in alphabet), for diagnostics
    """
    with open_maybe_gzip(fasta_path) as fh:
        seq_iter = parse_fasta(fh)
        counts_all = Counter()

        allowed_nonletters = set()
        if include_gaps:
            allowed_nonletters |= set("-.")
        if include_stop:
            allowed_nonletters |= set("*")

        for seq in seq_iter:
            # keep only letters and optionally allowed non-letters
            for ch in seq:
                if ch.isalpha():
                    counts_all[ch.upper()] += 1
                elif ch in allowed_nonletters:
                    counts_all[ch] += 1
                # else ignore digits/whitespace/others by default

    # Split letters vs non-letters (like '-' or '*') for alphabet detection
    letters_present = {k for k in counts_all if k.isalpha() and counts_all[k] > 0}

    if alphabet == "auto":
        resolved_name, resolved_set = detect_alphabet(letters_present)
    elif alphabet == "dna":
        resolved_name, resolved_set = "dna", DNA_IUPAC
    elif alphabet == "rna":
        resolved_name, resolved_set = "rna", RNA_IUPAC
    elif alphabet == "protein":
        resolved_name, resolved_set = "protein", PROTEIN_EXTENDED
    else:
        raise ValueError(f"Unknown alphabet option: {alphabet}")

    counts_in_alpha = Counter()
    counts_other = Counter()

    # Count only characters that belong to the chosen alphabet (letters)
    for sym, cnt in counts_all.items():
        if sym.isalpha():
            if sym in resolved_set or resolved_name == "auto":
                counts_in_alpha[sym] += cnt
            else:
                counts_other[sym] += cnt
        else:
            # Non-letter symbols are never part of the biological alphabet,
            # but we may still want to report them in 'Other'
            if cnt > 0:
                counts_other[sym] += cnt

    total_in_alpha = sum(counts_in_alpha.values())
    return counts_in_alpha, total_in_alpha, resolved_name, resolved_set, counts_other

def main():
    ap = argparse.ArgumentParser(description="Compute symbol percentages from a FASTA file.")
    ap.add_argument("fasta", help="Path to FASTA file (use '-' for STDIN). Supports .gz")
    ap.add_argument("--alphabet", choices=["auto", "dna", "rna", "protein"],
                    default="auto", help="Alphabet to use (auto-detect by default)")
    ap.add_argument("--include-gaps", action="store_true",
                    help="Also tally gap characters '-' and '.' in 'Other'")
    ap.add_argument("--include-stop", action="store_true",
                    help="Also tally stop '*' in 'Other'")
    ap.add_argument("--sort", choices=["symbol", "percent", "count"], default="percent",
                    help="How to sort the output rows")
    args = ap.parse_args()

    counts_in_alpha, total, resolved_name, resolved_set, counts_other = tally_symbols(
        args.fasta,
        alphabet=args.alphabet,
        include_gaps=args.include_gaps,
        include_stop=args.include_stop,
    )

    if total == 0:
        print("No sequence symbols found. Check that the FASTA contains sequence lines.")
        return

    # Prepare rows
    rows = []
    for sym, cnt in counts_in_alpha.items():
        pct = (cnt / total) * 100.0
        rows.append((sym, cnt, pct))

    if args.sort == "symbol":
        rows.sort(key=lambda x: x[0])
    elif args.sort == "count":
        rows.sort(key=lambda x: x[1], reverse=True)
    else:  # percent
        rows.sort(key=lambda x: x[2], reverse=True)

    print(f"# Alphabet: {resolved_name.upper()} (auto-detected)" if args.alphabet == "auto"
          else f"# Alphabet: {args.alphabet.upper()}")
    if resolved_set:
        print(f"# Symbols counted toward percentages: {''.join(sorted(resolved_set))}")
    print(f"# Total symbols counted: {total}\n")

    # Print table
    print(f"{'Symbol':<8}{'Count':>12}{'Percent':>12}")
    print("-" * 32)
    for sym, cnt, pct in rows:
        print(f"{sym:<8}{cnt:>12}{pct:>11.2f}%")

    # Report 'other' symbols, if any
    other_total = sum(counts_other.values())
    if other_total:
        print("\n# Characters outside the chosen alphabet (not included in %):")
        for sym in sorted(counts_other):
            print(f"  {sym!r}: {counts_other[sym]}")

if __name__ == "__main__":
    main()
