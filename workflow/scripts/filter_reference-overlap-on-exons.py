#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pybedtools
import logging
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def parse_args():
    parser = argparse.ArgumentParser(description="Remove isoforms that overlap reference features at the exon level, based on absolute or fractional overlap.")
    parser.add_argument("--isoform_bed12", required=True, help="BED12 file of isoform predictions")
    parser.add_argument("--reference_bed12", required=True, help="BED12 file of reference annotations (exons)")
    parser.add_argument("--out_bed", required=True, help="Output BED12 file of removed isoforms")
    parser.add_argument("--out_ids", required=True, help="Output list of IDs of removed isoforms")
    parser.add_argument("--out_stats", required=True, help="Output TSV file with stats for removed isoforms")
    parser.add_argument("--min_overlap_fraction", type=float, default=0.0, help="Minimum fraction of exon bases overlapping reference features to trigger removal (0 = any overlap)")
    parser.add_argument("--strand_specific", action="store_true", help="Enable strand-specific comparisons")
    return parser.parse_args()

def bed12_to_exons(bed12_entry):
    """Expand a BED12 entry into a list of exon-level BED intervals."""
    chrom = bed12_entry.chrom
    start = int(bed12_entry.start)
    strand = bed12_entry.strand
    name = bed12_entry.name

    block_sizes = list(map(int, bed12_entry.fields[10].strip(',').split(',')))
    block_starts = list(map(int, bed12_entry.fields[11].strip(',').split(',')))

    exons = []
    for i in range(len(block_sizes)):
        exon_start = start + block_starts[i]
        exon_end = exon_start + block_sizes[i]
        exons.append(pybedtools.create_interval_from_list([
            chrom, str(exon_start), str(exon_end), name, '0', strand
        ]))
    return exons

if __name__ == "__main__":
    args = parse_args()

    isoforms = pybedtools.BedTool(args.isoform_bed12)
    reference = pybedtools.BedTool(args.reference_bed12)
    reference_merged = reference.sort().merge(s=args.strand_specific)

    logging.info("Splitting isoform BED12 into individual exons...")
    exon_dict = defaultdict(list)
    all_exons = []

    for iso in isoforms:
        exons = bed12_to_exons(iso)
        exon_dict[iso.name].extend(exons)
        all_exons.extend(exons)

    exons_bed = pybedtools.BedTool(all_exons).saveas()

    overlaps = exons_bed.intersect(reference_merged, s=args.strand_specific, wo=True)

    # Calculate cumulative exon lengths and overlaps
    exon_lengths = defaultdict(int)
    overlap_lengths = defaultdict(int)

    for exon in exons_bed:
        exon_lengths[exon.name] += exon.length

    for feature in overlaps:
        tid = feature.name
        overlap_bp = int(feature[-1])
        overlap_lengths[tid] += overlap_bp

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
        for iso in isoforms:
            if iso.name in removed_ids:
                out_ids.write(iso.name + "\n")
                out_bed.write(str(iso))

    logging.info("Overlap-based filtering complete.")
