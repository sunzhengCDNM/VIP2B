#!/usr/bin/env python3
"""
convert_marisa.py  —  one-time conversion of VIP2B marisa_trie databases to TSV.

Run this ONCE per database before using vip2b-abundance (Rust).

Usage:
    python3 scripts/convert_marisa.py  \\
        database/8Enzyme/AlfI_BcgI.species.uniq.marisa  \\
        database/8Enzyme/AlfI_BcgI.species.uniq.tsv

The TSV format (two tab-separated columns, no header):
    tag40  <TAB>  gcf8

  tag40  = 40-character zero-padded tag sequence (marisa key, unchanged)
  gcf8   = 8-character GCF/genome identifier (one row per gcf per tag)

Multiple rows with the same tag40 are written when a tag maps to multiple genomes.

The script also reads the companion *.stat.xls if present — no conversion needed
since vip2b-abundance reads that file directly.
"""

import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',  help='Input .marisa file (RecordTrie format)')
    parser.add_argument('output', help='Output .tsv file')
    parser.add_argument('--format', default='8c',
        help='marisa_trie RecordTrie format string (default: 8c)')
    args = parser.parse_args()

    try:
        import marisa_trie
    except ImportError:
        sys.exit("ERROR: marisa_trie not installed.  Run: pip install marisa-trie")

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: Input file not found: {args.input}")

    print(f"INFO: Loading {args.input} ...", flush=True)
    trie = marisa_trie.RecordTrie(args.format).mmap(args.input)

    # Count keys for progress reporting
    all_keys = list(trie.keys())
    n_keys   = len(all_keys)
    print(f"INFO: {n_keys} unique tags in the trie", flush=True)

    written = 0
    with open(args.output, 'w', buffering=1 << 20) as OUT:
        for i, key in enumerate(all_keys):
            # Each record in RecordTrie('8c') is a tuple of 8 one-byte objects
            for record in trie[key]:
                gcf_id = ''.join(b.decode('utf-8') for b in record)
                OUT.write(f"{key}\t{gcf_id}\n")
                written += 1

            if (i + 1) % 500_000 == 0:
                pct = (i + 1) / n_keys * 100
                print(f"  {i+1}/{n_keys} tags ({pct:.1f}%)  rows written: {written}",
                      flush=True)

    print(f"INFO: Done.  {written} rows written to {args.output}", flush=True)


if __name__ == '__main__':
    main()
