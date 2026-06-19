#!/usr/bin/env bash
set -o nounset
set -o pipefail

cd "$1" 

if command -v module >/dev/null 2>&1; then
    module load samtools/1.16.1 || true
fi

ls

# md5 a file, skipping comment/header lines (Picard "# Started on:", command lines, etc.).
md5_nohdr() {
    echo "$1 $(grep -v '^#' "$1" | md5sum | cut -d' ' -f1)"
}

# Picard metrics (aggregate + per-interval): strip headers, then hash the tables.
for f in *.duplicate_metrics *.dupmarked.metrics \
         *.wgs_metrics.txt *.raw_wgs_metrics.txt \
         *.alignment_summary_metrics \
         *.gc_bias.summary_metrics *.gc_bias.detail_metrics \
         *.quality_distribution_metrics; do
    [ -e "$f" ] || continue
    md5_nohdr "$f"
done

# samtools stats: the "# The command line was:" line holds a temp path; strip comments.
for f in *.samtools_stats.txt; do
    [ -e "$f" ] || continue
    md5_nohdr "$f"
done

# Read-length distribution: plain counts, no header; hash directly.
for f in *.read_length_distribution.txt; do
    [ -e "$f" ] || continue
    echo "$f $(md5sum < "$f" | cut -d' ' -f1)"
done

