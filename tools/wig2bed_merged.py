#!/usr/bin/env python3

import sys

if len(sys.argv) != 4:
    print("Usage: python3 wig2bed_merged.py <wig_file> <max_occurrences> <out_bed>")
    sys.exit(1)

wig_file = sys.argv[1]
max_occurrences = float(sys.argv[2])
out_bed = sys.argv[3]

def flush_interval(chrom, start, end, score_sum, count, out_handle):
    """Write merged interval to BED"""
    if start is not None and count > 0:
        avg_score = score_sum / count
        out_handle.write("{}\t{}\t{}\t{:.6f}\n".format(chrom, start, end, avg_score))
    return None, None, 0, 0

with open(wig_file) as f, open(out_bed, "w") as out:
    chrom = ""
    merge_start = None
    merge_end = None
    merge_score_sum = 0
    merge_count = 0

    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith("variableStep"):
            # New chromosome: flush previous interval
            merge_start, merge_end, merge_score_sum, merge_count = flush_interval(
                chrom, merge_start, merge_end, merge_score_sum, merge_count, out
            )
            parts = line.split()
            chrom = [p.split("=")[1] for p in parts if p.startswith("chrom=")][0]
            continue

        # Each line: "position score"
        try:
            pos_str, score_str = line.split()
        except ValueError:
            print("Skipping invalid line: {}".format(line), file=sys.stderr)
            continue

        pos = int(pos_str) - 1  # BED is 0-based
        score = float(score_str)

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
                    chrom, merge_start, merge_end, merge_score_sum, merge_count, out
                )
                merge_start = pos
                merge_end = pos + 1
                merge_score_sum = score
                merge_count = 1
        else:
            merge_start, merge_end, merge_score_sum, merge_count = flush_interval(
                chrom, merge_start, merge_end, merge_score_sum, merge_count, out
            )

    # Flush last interval
    flush_interval(chrom, merge_start, merge_end, merge_score_sum, merge_count, out)
