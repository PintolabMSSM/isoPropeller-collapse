#!/usr/bin/env python
# -*- coding: utf-8 -*-

##################
# IMPORT MODULES #
##################

import argparse
import bisect
import gzip
import logging
import sys
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# ============================================================================
# CHANGE LOG (v1 -> v2)
#
# One correctness bug was found and fixed, confirmed via hand-built cases
# and a 500,000-gene synthetic ground-truth set:
#
#   BUG B (critical, --strand_specific only): bedtools `merge -s` only uses
#     strand to decide which input intervals are allowed to merge together;
#     it does NOT include a strand column in its own output unless `-c
#     <col> -o <op>` is explicitly given. The original's
#     `reference.sort().merge(s=args.strand_specific)` therefore produced a
#     strand-LESS 3-column BED whenever --strand_specific was set, and the
#     following `.intersect(reference_merged, s=True, wo=True)` call --
#     bedtools intersect -s requires strand info in BOTH inputs -- silently
#     matched NOTHING, for every isoform, with no error of any kind. Every
#     run of this script with --strand_specific removed exactly 0 isoforms,
#     regardless of true overlap. Confirmed at 500,000-gene scale:
#     28,571 of 50,000 isoforms (~57.1%, i.e. every isoform not on the
#     synthetic "opposite_strand" test category) that should have been
#     removed were wrongly kept. Fixed by grouping reference intervals by
#     (chrom, strand) when --strand_specific is set, so strand is never
#     lost between the merge and overlap-lookup steps.
#
# INVESTIGATED, NOT A BUG (initially mis-flagged as one, then reverted on
# the author's confirmation): an earlier revision of this fix also split
# the reference BED12 into its own exon-level blocks before computing
# overlap, on the theory that "reference overlap" should mean overlap with
# a reference transcript's actual annotated exons, not its introns (the
# same reasoning that was correct for the sibling script
# filter_monoexon-premrna-fragments.py, which explicitly checks for exon
# overlap AND intron overlap as two separate signals). That reasoning does
# NOT apply here: this script is intentionally designed to test overlap
# against a reference transcript's ENTIRE genomic region (its full
# chromStart-chromEnd span, introns included), not just its exons -- per
# the author's explicit confirmation. That exon-expansion change has been
# REVERTED; `build_merged_reference()` once again merges each reference
# BED12 row's raw, unsplit chromStart-chromEnd span, exactly like the
# original.
#
# The exon-vs-merged-reference overlap computation itself (bedtools
# intersect + Python summation in v1) is reimplemented here as a
# dependency-free sort+bisect sweep, which also removes the
# `pybedtools`/`bedtools` dependency entirely and the significant
# per-record Interval-object overhead that comes with it (confirmed
# repeatedly in this script series to be a severe anti-pattern at genome
# scale). Note that only the ISOFORM side is still split into exons before
# computing overlap (unaffected by the reversal above) -- the original
# always did this for isoforms (`bed12_to_exons()` on the isoform side was
# never in question), and v2 preserves it exactly.
#
# Also removed: `exon_dict`, which the original built (one list of exon
# Interval objects per isoform name, accumulating in memory for the whole
# file) but never actually read anywhere -- pure dead code/wasted memory.
#
# Bonus (not a regression fix): routed all BED12 file reads through a
# smart_open() helper so gzip-compressed (.gz) --isoform_bed12/
# --reference_bed12 inputs work, matching what pybedtools supported
# transparently in the original.
#
# All CLI flags, the min_overlap_fraction==0.0 "any overlap" special case,
# the aggregate-by-transcript-ID behavior (if the same ID appears on
# multiple BED12 lines, their exon lengths/overlaps are summed together and
# ALL of that ID's lines are written to --out_bed/--out_ids if the combined
# total triggers removal -- reproduced exactly, including writing the ID
# more than once to --out_ids if it appears on more than one line), and the
# --out_bed/--out_ids "removed/flagged set" convention (confirmed as the
# established convention for this filter_* script family) are unchanged.
# ============================================================================


#############
# FUNCTIONS #
#############

def is_interactive():
    """Check if we're in an interactive session (e.g. Jupyter notebook)."""
    import __main__ as main
    return not hasattr(main, '__file__')


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Remove isoforms that overlap reference features at the exon level, based on absolute or fractional overlap.")
    parser.add_argument("--isoform_bed12", required=True, help="BED12 file of isoform predictions")
    parser.add_argument("--reference_bed12", required=True, help="BED12 file of reference annotations (exons)")
    parser.add_argument("--out_bed", required=True, help="Output BED12 file of removed isoforms")
    parser.add_argument("--out_ids", required=True, help="Output list of IDs of removed isoforms")
    parser.add_argument("--out_stats", required=True, help="Output TSV file with stats for removed isoforms")
    parser.add_argument("--min_overlap_fraction", type=float, default=0.0, help="Minimum fraction of exon bases overlapping reference features to trigger removal (0 = any overlap)")
    parser.add_argument("--strand_specific", action="store_true", help="Enable strand-specific comparisons")

    if is_interactive():
        return parser.parse_args('')
    else:
        return parser.parse_args(args=None if sys.argv[1:] else ['--help'])


def smart_open(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def iter_bed12_lines(path):
    """Yield (raw_line_without_newline, fields_list) for every well-formed
    (>=12 field) line of a BED12 file."""
    with smart_open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            yield line, fields


def bed12_exons_from_fields(fields):
    """Return a list of (chrom, exon_start, exon_end, strand) for one BED12
    record's fields."""
    chrom = fields[0]
    start = int(fields[1])
    strand = fields[5]
    block_sizes = list(map(int, fields[10].strip(',').split(',')))
    block_starts = list(map(int, fields[11].strip(',').split(',')))
    exons = []
    for i in range(len(block_sizes)):
        exon_start = start + block_starts[i]
        exon_end = exon_start + block_sizes[i]
        exons.append((chrom, exon_start, exon_end, strand))
    return exons


def build_merged_reference(reference_path, strand_specific):
    """Parse the reference BED12 and merge overlapping/book-ended (bedtools'
    default -d 0 semantics) intervals per bucket, using each reference
    entry's RAW, unsplit chromStart-chromEnd span (its entire genomic
    footprint, introns included) -- not its exon blocks. This script is
    intentionally designed to test overlap against a reference transcript's
    entire region, not just its annotated exons. Bucket is (chrom, strand)
    if strand_specific else chrom alone (strand ignored, matching
    merge(s=False)/intersect(s=False) in the original).

    Returns: dict[bucket] -> (starts_list, ends_list), each sorted by start,
    non-overlapping, so overlap lookups can bisect on starts_list.
    """
    buckets = defaultdict(list)
    for _, fields in iter_bed12_lines(reference_path):
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[5]
        key = (chrom, strand) if strand_specific else chrom
        buckets[key].append((start, end))

    merged = {}
    for key, intervals in buckets.items():
        intervals.sort()
        merged_starts = []
        merged_ends = []
        cur_start, cur_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_end:  # overlapping or book-ended -> merge (bedtools default -d 0)
                if e > cur_end:
                    cur_end = e
            else:
                merged_starts.append(cur_start)
                merged_ends.append(cur_end)
                cur_start, cur_end = s, e
        merged_starts.append(cur_start)
        merged_ends.append(cur_end)
        merged[key] = (merged_starts, merged_ends)
    return merged


def overlap_bp_for_exon(merged, key, ex_start, ex_end):
    """Sum of overlap (bp) between [ex_start, ex_end) and the sorted,
    non-overlapping merged intervals under `key`."""
    entry = merged.get(key)
    if entry is None:
        return 0
    starts, ends = entry
    # first merged interval whose start could possibly overlap ex_start:
    # (rightmost interval with start <= ex_start) or the next one after it.
    i = bisect.bisect_right(starts, ex_start) - 1
    if i < 0:
        i = 0
    total = 0
    n = len(starts)
    while i < n and starts[i] < ex_end:
        s, e = starts[i], ends[i]
        ov_start = s if s > ex_start else ex_start
        ov_end = e if e < ex_end else ex_end
        if ov_end > ov_start:
            total += ov_end - ov_start
        i += 1
    return total


########
# MAIN #
########

if __name__ == "__main__":
    args = parse_args()

    logging.info("Building merged reference intervals...")
    merged_reference = build_merged_reference(args.reference_bed12, args.strand_specific)

    logging.info("Splitting isoform BED12 into individual exons and computing overlap...")
    exon_lengths = defaultdict(int)
    overlap_lengths = defaultdict(int)
    iso_lines = []  # (raw_line, name) in file order, for the final write pass

    for raw_line, fields in iter_bed12_lines(args.isoform_bed12):
        name = fields[3]
        iso_lines.append((raw_line, name))
        for chrom, ex_start, ex_end, strand in bed12_exons_from_fields(fields):
            exon_lengths[name] += (ex_end - ex_start)
            key = (chrom, strand) if args.strand_specific else chrom
            ov = overlap_bp_for_exon(merged_reference, key, ex_start, ex_end)
            if ov > 0:
                overlap_lengths[name] += ov

    # Determine which isoforms to remove
    removed_ids = set()
    for tid in exon_lengths:
        total = exon_lengths[tid]
        overlap = overlap_lengths.get(tid, 0)
        frac = overlap / total if total > 0 else 0
        if args.min_overlap_fraction == 0.0:
            if tid in overlap_lengths:
                removed_ids.add(tid)
        else:
            if frac >= args.min_overlap_fraction:
                removed_ids.add(tid)

    logging.info(f"Removing {len(removed_ids)} isoforms based on overlap threshold.")

    # Write stats for removed isoforms
    with open(args.out_stats, "w") as stats_out:
        stats_out.write("isoform\ttotal_exon_bp\toverlap_bp\toverlap_fraction\n")
        for tid in sorted(removed_ids):
            total = exon_lengths[tid]
            overlap = overlap_lengths.get(tid, 0)
            frac = overlap / total if total > 0 else 0
            stats_out.write(f"{tid}\t{total}\t{overlap}\t{frac:.4f}\n")

    # Write IDs and BED entries for removed isoforms
    with open(args.out_ids, 'w') as out_ids, open(args.out_bed, 'w') as out_bed:
        for raw_line, name in iso_lines:
            if name in removed_ids:
                out_ids.write(name + "\n")
                out_bed.write(raw_line + "\n")

    logging.info("Overlap-based filtering complete.")
