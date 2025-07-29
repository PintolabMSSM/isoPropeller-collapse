#!/usr/bin/env python
# -*- coding: utf-8 -*-

##################
# IMPORT MODULES #
##################

import argparse
import pybedtools
import logging
import sys
import os
from collections import defaultdict

logging.basicConfig(
    level   = logging.INFO,
    format  = '%(asctime)s - %(levelname)s - %(message)s',
    datefmt = '%Y-%m-%d %H:%M:%S'
)


#############
# FUNCTIONS #
#############


def is_interactive():
    """Check if we're in an interactive session (e.g. Jupyter notebook)."""
    import __main__ as main
    return not hasattr(main, '__file__')


def parse_args():
    parser = argparse.ArgumentParser(description="Remove isoforms whose splice chains overlap reference isoforms on the antisense strand.")
    parser.add_argument("--isoform_bed12", required=True, help="BED12 file of isoform predictions")
    parser.add_argument("--reference_bed12", required=True, help="BED12 file of reference annotations")
    parser.add_argument("--out_bed", required=True, help="Output BED12 file of isoforms removed due to antisense splice chain overlap")
    parser.add_argument("--out_ids", required=True, help="Output file listing transcript IDs removed due to antisense splice chain overlap")

    if is_interactive():
        args = parser.parse_args('')
    else:
        args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    return args


def extract_splice_junctions(bed):
    """Extract splice junctions as tuples: (chrom, donor, acceptor, strand, transcript_id)."""
    splice_junctions = defaultdict(list)
    for entry in bed:
        chrom = entry.chrom
        strand = entry.strand
        tid = entry.name
        start = int(entry.start)
        block_sizes = list(map(int, entry.fields[10].strip(',').split(',')))
        block_starts = list(map(int, entry.fields[11].strip(',').split(',')))

        for i in range(len(block_starts) - 1):
            donor = start + block_starts[i] + block_sizes[i]
            acceptor = start + block_starts[i + 1]
            splice_junctions[tid].append((chrom, donor, acceptor, strand))
    return splice_junctions


def get_splice_chain_set(splice_junctions):
    """Convert list of junctions to set of (chrom, donor, acceptor) tuples ignoring strand."""
    chain_set = set()
    for sj in splice_junctions:
        chain_set.add((sj[0], sj[1], sj[2]))
    return chain_set


########
# MAIN #
########

if __name__ == "__main__":
    args = parse_args()

    isoform_bed = pybedtools.BedTool(args.isoform_bed12)
    reference_bed = pybedtools.BedTool(args.reference_bed12)

    logging.info("Extracting splice junctions...")
    isoform_junctions = extract_splice_junctions(isoform_bed)
    reference_junctions = extract_splice_junctions(reference_bed)

    # Build antisense splice chain set
    antisense_chain_set = set()
    for sjs in reference_junctions.values():
        if not sjs:
            continue
        strand = sjs[0][3]
        if strand == "+":
            antisense = "-"
        elif strand == "-":
            antisense = "+"
        else:
            continue
        flipped = set((c, d, a) for (c, d, a, s) in sjs if s == antisense)
        antisense_chain_set.update(flipped)

    removed_ids = set()
    with open(args.out_bed, "w") as out_bed, open(args.out_ids, "w") as out_ids:
        for entry in isoform_bed:
            tid = entry.name
            if tid not in isoform_junctions:
                continue  # skip isoforms without junctions
            sj_set = get_splice_chain_set(isoform_junctions[tid])
            if sj_set & antisense_chain_set:
                out_bed.write(str(entry))
                removed_ids.add(tid)
                out_ids.write(tid + "\\n")

    logging.info(f"Removed {len(removed_ids)} isoforms with antisense splice chain overlap.")
