#!/usr/bin/env python3
"""
quan_from_tags.py — Quantitative abundance from pre-computed unique tag counts.

Called by the Rust vip2b binary in Phase 5.  Rust reads the digest .fa.gz,
deduplicates tags into a unique tag→count map, and pipes lines of the form

    tag40 <TAB> count

to this script's stdin.  This script loads the per-sample marisa DB natively,
does the C++ marisa lookups for only the unique tags (not all 21M raw reads),
aggregates to taxonomy level, and writes the standard .xls output.

No marisa→TSV conversion is ever performed.

Usage:
    python3 quan_from_tags.py \\
        -d mkdb/sample/sample  \\   # DB prefix (.marisa + .stat.xls)
        -c classify.txt.gz     \\   # abfh_classify file
        -t Species             \\   # taxonomy level
        -o quan/               \\   # output parent dir (sample subdir created)
        -s sample_id           \\   # sample name
        -ct 0.6                     # coverage threshold (0 = keep all)
"""

import sys, os, argparse, gzip, math
import marisa_trie

# ── Taxonomy level → column index (1-based, slice to col+1) ──────────────────
LEVEL_COL = {'Class': 3, 'Order': 4, 'Family': 5, 'Genus': 6, 'Species': 7}

def open_maybe_gz(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def load_stat(stat_path):
    """Read {prefix}.stat.xls → (gcf_theo_dic, tax_theo_dic)."""
    gcf_theo, tax_theo = {}, {}
    with open(stat_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            tmp = line.split('\t')
            if len(tmp) < 4:
                continue
            gcf, tax = tmp[0], tmp[1]
            try:
                gcf_theo[gcf] = float(tmp[2])
                tax_theo[tax] = float(tmp[3])
            except ValueError:
                pass
    return gcf_theo, tax_theo

def load_classify(path, level):
    """Read classify file → {gcf: taxonomy_string}."""
    col = LEVEL_COL.get(level, 7)
    gcf_tax = {}
    with open_maybe_gz(path) as f:
        for line in f:
            tmp = line.rstrip('\n').split('\t')
            if not tmp or tmp[0].startswith('#'):
                continue
            gcf = tmp[0]
            # Match Python original: ','.join(tmp[1 : col+1])
            gcf_tax.setdefault(gcf, ','.join(tmp[1:col + 1]))
    return gcf_tax

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d',  dest='database',  required=True,
                        help='per-sample mkdb prefix (without extension)')
    parser.add_argument('-c',  dest='classify',  required=True,
                        help='abfh_classify_with_speciename.txt.gz')
    parser.add_argument('-t',  dest='level',     default='Species',
                        choices=['Class', 'Order', 'Family', 'Genus', 'Species'],
                        help='taxonomy level')
    parser.add_argument('-o',  dest='output',    required=True,
                        help='output parent directory')
    parser.add_argument('-s',  dest='sample',    required=True,
                        help='sample ID')
    parser.add_argument('-ct', dest='threshold', type=float, default=0.0,
                        help='coverage threshold (default 0 = keep all)')
    args = parser.parse_args()

    # ── Load per-sample marisa DB (mmap — fast, no full copy into RAM) ────────
    marisa_path = args.database + '.marisa'
    stat_path   = args.database + '.stat.xls'

    if not os.path.exists(marisa_path):
        sys.stderr.write(f"ERROR: {marisa_path} not found\n")
        sys.exit(1)

    tag_trie = marisa_trie.RecordTrie('8c').mmap(marisa_path)
    gcf_theo, tax_theo = load_stat(stat_path)
    gcf_tax  = load_classify(args.classify, args.level)

    # ── Read unique tag counts from Rust via stdin ────────────────────────────
    # Format: tag40 <TAB> count  (one line per unique zero-padded tag)
    # gcf_tag_dic: gcf → { tag40: count }  (max-count deduplication)
    gcf_tag_dic = {}

    for line in sys.stdin:
        line = line.rstrip('\n')
        if not line:
            continue
        parts = line.split('\t', 1)
        if len(parts) != 2:
            continue
        tag40 = parts[0]
        try:
            count = int(parts[1])
        except ValueError:
            continue

        try:
            gcf_records = tag_trie[tag40]
        except KeyError:
            continue

        for rec in gcf_records:
            gcf = rec[0].decode('utf-8', errors='replace').rstrip()
            entry = gcf_tag_dic.setdefault(gcf, {})
            if tag40 not in entry or count > entry[tag40]:
                entry[tag40] = count

    # ── Aggregate GCF → taxonomy level ───────────────────────────────────────
    tax_tag_dic  = {}   # tax → {tag40: max_count}
    for gcf, tag_counts in gcf_tag_dic.items():
        tax = gcf_tax.get(gcf)
        if tax is None:
            continue
        t = tax_tag_dic.setdefault(tax, {})
        for tag, cnt in tag_counts.items():
            if tag not in t or cnt > t[tag]:
                t[tag] = cnt

    # ── Write output ──────────────────────────────────────────────────────────
    outdir = os.path.join(args.output, args.sample)
    os.makedirs(outdir, exist_ok=True)

    xls_path = os.path.join(outdir, f"{args.sample}.xls")
    with open(xls_path, 'w') as OUT:
        OUT.write(
            'Taxonomy\tTheoretical_Tag_Num\tSequenced_Tag_Num\tPercent\t'
            'Sequenced_Reads_Num\tSequenced_Reads_Num/Theoretical_Tag_Num\t'
            'Sequenced_Reads_Num/Sequenced_Tag_Num\t'
            'Sequenced_Tag_Num(depth>1)\tG_Score\n'
        )
        for tax in sorted(tax_tag_dic.keys()):
            tag_counts = tax_tag_dic[tax]
            seq_tag_num  = len(tag_counts)
            seq_read_num = sum(tag_counts.values())
            theo         = tax_theo.get(tax, 1.0)
            percent      = round(seq_tag_num / theo, 4)
            if percent < args.threshold:
                continue
            g_score = round(math.sqrt(seq_tag_num * seq_read_num), 4)
            OUT.write(
                f"{tax}\t{theo}\t{seq_tag_num}\t{percent}\t{seq_read_num}\t"
                f"{round(seq_read_num / theo, 4)}\t"
                f"{round(seq_read_num / seq_tag_num, 4)}\t"
                f"NA\t{g_score}\n"
            )

    sys.stdout.write(f"INFO [{args.sample}]: quan_from_tags done → {xls_path}\n")

if __name__ == '__main__':
    main()
