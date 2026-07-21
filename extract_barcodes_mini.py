#!/usr/bin/env python3
"""
extract_barcodes.py

Extract barcodes flanked by known constant sequences from FASTQ reads,
count occurrences, and report a CSV summary.

The barcode is defined as whatever sequence lies between:
    LEFT flank  : ttgcagagctca
    RIGHT flank : aatacagctccc

Both the forward read sequence and its reverse complement are scanned,
so barcodes are found regardless of which strand/orientation they were
sequenced on.

Usage
-----
    python extract_barcodes.py -i reads.fastq -o barcode_counts.csv
    python extract_barcodes.py -i R1.fastq.gz -i2 R2.fastq.gz -o out.csv
    python extract_barcodes.py -i reads.fastq -o out.csv \
        --left TTGCAGAGCTCA --right AATACAGCTCCC \
        --min-len 10 --max-len 30 --mismatches 1

Notes
-----
- Reads can be plain text or gzip-compressed (.gz), autodetected by extension.
- You can pass one file (-i) or a pair of paired-end files (-i / -i2); both
  are searched identically (fw + rv complement), so R1/R2 don't need to be
  distinguished -- a barcode read from either mate/orientation is captured.
- Fuzzy flank matching (allowing N mismatches) requires the third-party
  `regex` package. If it isn't installed, the script automatically falls
  back to exact matching and prints a note. Install it with:
      pip install regex --break-system-packages
"""

import argparse
import csv
import gzip
import sys
from collections import Counter

# Try to use the third-party `regex` module for fuzzy (mismatch-tolerant)
# matching. Fall back to the standard `re` module (exact matches only).
try:
    import regex as re_engine
    HAVE_FUZZY = True
except ImportError:
    import re as re_engine
    HAVE_FUZZY = False

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def fastq_reader(path):
    """Yield (header, sequence, plus, quality) tuples from a FASTQ file."""
    with open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not (header and seq and plus and qual):
                break
            yield header.rstrip("\n"), seq, plus, qual


def build_pattern(left, right, min_len, max_len, mismatches):
    """
    Build a compiled regex pattern that captures the barcode between
    the left and right flanking sequences.
    """
    left = left.upper()
    right = right.upper()
    barcode_wildcard = f"[ACGTN]{{{min_len},{max_len}}}"

    if mismatches > 0:
        if not HAVE_FUZZY:
            print(
                "WARNING: --mismatches > 0 requested but the 'regex' package "
                "is not installed. Falling back to exact flank matching.\n"
                "Install fuzzy support with: pip install regex --break-system-packages",
                file=sys.stderr,
            )
            mismatches = 0

    if mismatches > 0:
        # NOTE: these MUST be non-capturing groups (?:...). If they were
        # capturing groups, the match would contain 3 groups (left, barcode,
        # right) instead of 1, and code that grabs "the last group" as the
        # barcode would silently grab the right-flank match instead -- which
        # is exactly the bug that produced barcode_seq values identical/near-
        # identical to the right flank sequence.
        left_pat = f"(?:{left}){{e<={mismatches}}}"
        right_pat = f"(?:{right}){{e<={mismatches}}}"
    else:
        left_pat = left
        right_pat = right

    pattern = f"{left_pat}(?P<barcode>{barcode_wildcard}){right_pat}"

    if mismatches > 0:
        # IMPORTANT: without BESTMATCH, the `regex` engine returns the first
        # fuzzy match satisfying the edit-distance tolerance, not the one with
        # the fewest errors. Since the flanks can be a near-match to sequence
        # just inside the barcode (e.g. a single-base deletion lets "AATACAG..."
        # fuzzy-match "ATACAG..."), this can silently steal bases from the
        # barcode into the flank and return truncated/shifted "barcodes".
        # BESTMATCH forces it to pick the true minimal-edit-distance match.
        return re_engine.compile(pattern, re_engine.BESTMATCH)
    return re_engine.compile(pattern)


def scan_read(seq, pattern):
    """Return all barcode matches found in a sequence (non-overlapping)."""
    seq = seq.upper()
    return [m.group("barcode") for m in pattern.finditer(seq)]


def is_nnk(seq):
    """
    Check whether a sequence conforms to an NNK codon pattern: each codon
    is N-N-K, where N = any base and K = G or T (keto base) at the wobble
    (3rd) position. Requires len(seq) to be a multiple of 3; any 'N'
    (ambiguous base call) at a wobble position also fails the check, since
    it can't be confirmed as G/T.
    """
    seq = seq.upper()
    if len(seq) % 3 != 0:
        return False
    for i in range(0, len(seq), 3):
        wobble = seq[i + 2]
        if wobble not in ("G", "T"):
            return False
    return True



STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_stop_codons(seq):
    """
    Return the list of in-frame stop codons found in a sequence (3 bases at
    a time from the start), in order. Empty list if none, or if len(seq)
    isn't a multiple of 3.
    """
    seq = seq.upper()
    if len(seq) % 3 != 0:
        return []
    return [seq[i:i + 3] for i in range(0, len(seq), 3) if seq[i:i + 3] in STOP_CODONS]


def has_stop_codon(seq):
    """Check whether any in-frame codon is a stop codon (TAA, TAG, TGA)."""
    return bool(find_stop_codons(seq))


def gc_content(seq):
    """GC content of a sequence, as a percentage rounded to 1 decimal place."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in ("G", "C"))
    return round(100 * gc / len(seq), 1)


def max_homopolymer_run(seq):
    """
    Return (run_length, repeated_sequence) for the longest run of a single
    repeated base in the sequence, e.g. "AAAAA" -> (5, "AAAAA").
    """
    seq = seq.upper()
    if not seq:
        return 0, ""
    max_run = 1
    max_start = 0
    cur_run = 1
    cur_start = 0
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur_run += 1
        else:
            cur_run = 1
            cur_start = i
        if cur_run > max_run:
            max_run = cur_run
            max_start = cur_start
    return max_run, seq[max_start:max_start + max_run]


def process_file(path, pattern, counter, stats):
    for header, seq, plus, qual in fastq_reader(path):
        stats["total_reads"] += 1
        found_this_read = False

        for orientation_seq in (seq, revcomp(seq)):
            matches = scan_read(orientation_seq, pattern)
            for bc in matches:
                counter[bc] += 1
                found_this_read = True

        if found_this_read:
            stats["reads_with_barcode"] += 1


def main():
    parser = argparse.ArgumentParser(
        description="Extract flanked barcodes from FASTQ reads and count them."
    )
    parser.add_argument(
        "-i", "--input", required=True, action="append",
        help="Input FASTQ file (.fastq or .fastq.gz). Repeatable for multiple files "
             "(e.g. -i R1.fastq.gz -i R2.fastq.gz).",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output CSV file path.",
    )
    parser.add_argument(
        "--left", default="ttgcagagctca",
        help="Left flanking sequence (default: ttgcagagctca).",
    )
    parser.add_argument(
        "--right", default="aatacagctccc",
        help="Right flanking sequence (default: aatacagctccc).",
    )
    parser.add_argument(
        "--min-len", type=int, default=1,
        help="Minimum expected barcode length (default: 1). Ignored if --length is set.",
    )
    parser.add_argument(
        "--max-len", type=int, default=200,
        help="Maximum expected barcode length (default: 200). Ignored if --length is set.",
    )
    parser.add_argument(
        "--length", type=int, default=None,
        help="Exact known/designed barcode length (e.g. 21 for an NNKx7 library). "
             "Strongly recommended whenever the barcode length is fixed by design: "
             "it sets --min-len/--max-len to this value and prevents the search from "
             "latching onto short, spurious matches near the flanks (a real failure "
             "mode when --mismatches > 0). Overrides --min-len/--max-len.",
    )
    parser.add_argument(
        "--length-tolerance", type=int, default=0,
        help="Allow barcode length to vary by +/- this many bases around --length "
             "(useful if the library may contain occasional indels; default: 0, "
             "i.e. exact length only). Only used together with --length.",
    )
    parser.add_argument(
        "--mismatches", type=int, default=0,
        help="Allowed mismatches in each flanking sequence (requires the "
             "'regex' package; default: 0, i.e. exact match).",
    )
    parser.add_argument(
        "--check_nnk", type=int, default=None,
        help="If set, adds 'NNK_pattern', 'Stop_codon', and 'Stop_codon_seq' columns "
             "to the output CSV. NNK_pattern checks each unique barcode against the "
             "NNK degenerate-codon pattern (N-N-K per codon, where K = G or T). "
             "Stop_codon (Yes/No) checks whether any in-frame codon is a stop codon "
             "(TAA/TAG/TGA); Stop_codon_seq lists the actual stop codon(s) found, "
             "comma-separated (empty if none). All checks are only applied to "
             "barcodes whose length equals this value (e.g. --check_nnk 21 checks "
             "only 21nt barcodes); barcodes of any other length are marked 'NA' in "
             "these columns, since the checks aren't meaningful for them.",
    )
    args = parser.parse_args()

    min_len, max_len = args.min_len, args.max_len
    if args.length is not None:
        min_len = args.length - args.length_tolerance
        max_len = args.length + args.length_tolerance
        print(
            f"Using fixed barcode length {args.length} "
            f"(+/-{args.length_tolerance}) -> searching range [{min_len}, {max_len}]",
            file=sys.stderr,
        )

    pattern = build_pattern(
        args.left, args.right, min_len, max_len, args.mismatches
    )

    counter = Counter()
    stats = {"total_reads": 0, "reads_with_barcode": 0}

    for path in args.input:
        print(f"Processing {path} ...", file=sys.stderr)
        process_file(path, pattern, counter, stats)

    with open(args.output, "w", newline="") as out_fh:
        writer = csv.writer(out_fh)
        base_header = [
            "barcode_seq", "length", "read_counts",
            "GC_content", "Max_homopolymer_run", "Max_homopolymer_seq",
        ]
        if args.check_nnk is not None:
            writer.writerow(base_header + ["NNK_pattern", "Stop_codon", "Stop_codon_seq"])
            for bc, count in counter.most_common():
                homo_len, homo_seq = max_homopolymer_run(bc)
                row = [bc, len(bc), count, gc_content(bc), homo_len, homo_seq]
                if len(bc) == args.check_nnk:
                    nnk_flag = "Yes" if is_nnk(bc) else "No"
                    stops = find_stop_codons(bc)
                    stop_flag = "Yes" if stops else "No"
                    stop_seq = ",".join(stops) if stops else ""
                else:
                    nnk_flag = "NA"
                    stop_flag = "NA"
                    stop_seq = "NA"
                writer.writerow(row + [nnk_flag, stop_flag, stop_seq])
        else:
            writer.writerow(base_header)
            for bc, count in counter.most_common():
                homo_len, homo_seq = max_homopolymer_run(bc)
                writer.writerow([bc, len(bc), count, gc_content(bc), homo_len, homo_seq])

    print(
        f"\nDone.\n"
        f"  Total reads processed : {stats['total_reads']}\n"
        f"  Reads with >=1 barcode: {stats['reads_with_barcode']}\n"
        f"  Unique barcodes found : {len(counter)}\n"
        f"  Total barcode calls   : {sum(counter.values())}\n"
        f"  Output written to     : {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
