#!/usr/bin/env python3

import sys
import os

def usage():
    print("Usage: python wig2bed_merged.py <wig_file> <max_occurrences> <out_bed>")
    sys.exit(1)

# -------------------------------
# Argument parsing
# -------------------------------
if len(sys.argv) != 4:
    usage()

wig_file = sys.argv[1]
try:
    max_occurrences = float(sys.argv[2])
except ValueError:
    print("Error: max_occurrences must be a number")
    sys.exit(1)

out_bed = sys.argv[3]

if not os.path.isfile(wig_file):
    print(f"Error: WIG file not found: {wig_file}")
    sys.exit(1)

# -------------------------------
# Main processing
# -------------------------------
with open(wig_file) as f, open(out_bed, "w") as out:

    chrom = None
    merge_start = None
    merge_end = None
    merge_score_sum = 0
    merge_count = 0

    def flush_interval(chrom, merge_start, merge_end, merge_score_sum, merge_count):
        """Write a merged interval to BED if valid"""
        if merge_start is not None and merge_count > 0 and chrom:
            avg_score = merge_score_sum / merge_count
            out.write(f"{chrom}\t{merge_start}\t{merge_end}\t{avg_score:.6f}\n")
        return None, None, 0, 0

    for line in f:
        line = line.strip()
        if not line:
            continue

        # New chromosome
        if line.startswith("variableStep"):
            merge_start, merge_end, merge_score_sum, merge_count = flush_interval(
                chrom, merge_start, merge_end, merge_score_sum, merge_count
            )
            parts = line.split()
            chrom_parts = [p.split("=")[1] for p in parts if p.startswith("chrom=")]
            chrom = chrom_parts[0] if chrom_parts else None
            continue

        # Lines with position and score
        try:
            pos_str, score_str = line.split()
            pos = int(pos_str) - 1  # BED is 0-based
            score = float(score_str)
        except ValueError:
            print(f"Skipping invalid line: {line}", file=sys.stderr)
            continue

        # Calculate occurrences
        occurrences = 1 / score if score > 0 else float('inf')

        if occurrences <= max_occurrences:
            # Merge consecutive bases
            if merge_start is None:
                merge_start = pos
                merge_end = pos + 1
                merge_score_sum = score
                merge_count = 1
            elif pos == merge_end:
                merge_end = pos + 1
                merge_score_sum += score
                merge_count += 1
            else:
                merge_start, merge_end, merge_score_sum, merge_count = flush_interval(
                    chrom, merge_start, merge_end, merge_score_sum, merge_count
                )
                merge_start = pos
                merge_end = pos + 1
                merge_score_sum = score
                merge_count = 1
        else:
            merge_start, merge_end, merge_score_sum, merge_count = flush_interval(
                chrom, merge_start, merge_end, merge_score_sum, merge_count
            )

    # Flush last interval if any
    merge_start, merge_end, merge_score_sum, merge_count = flush_interval(
        chrom, merge_start, merge_end, merge_score_sum, merge_count
    )
