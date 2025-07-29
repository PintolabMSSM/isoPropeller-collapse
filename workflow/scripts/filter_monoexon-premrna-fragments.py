#!/usr/bin/env python
# -*- coding: utf-8 -*-

##################
# IMPORT MODULES #
##################

import argparse
import pybedtools
import sys
import os
import gc
import logging
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
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Remove mono-exonic isoforms that likely represent pre-mRNA fragments (overlap reference exon AND intron, unless matching a reference mono-exon)."
    )
    parser.add_argument("--isoform_bed12",       required=True, help="BED12 file of isoform predictions")
    parser.add_argument("--reference_bed12",     required=True, help="Reference annotations in BED12 format")
    parser.add_argument("--out_bed",             required=True, help="Output BED12 of filtered isoforms")
    parser.add_argument("--out_ids",             required=True, help="Transcript IDs flagged as likely pre-mRNA")
    parser.add_argument("--min_intron_overlap",  type=int, default=10, help="Minimum base pairs of overlap with a reference intron to flag (default: 10)")
    
    if is_interactive():
        return parser.parse_args('')
    else:
        return parser.parse_args(args=None if sys.argv[1:] else ['--help'])


def extract_introns_from_bed12(bed):
    """
    Extract intron regions from BED12 entries.
    Assumes that entries have blockCount, blockSizes, and blockStarts fields.
    Returns a BedTool object with intron intervals.
    """
    intron_lines = []
    for entry in bed:
        block_count = int(entry[9])
        block_sizes = list(map(int, entry[10].strip(',').split(',')))
        block_starts = list(map(int, entry[11].strip(',').split(',')))

        if block_count < 2:
            continue

        for i in range(block_count - 1):
            intron_start = int(entry.start) + block_starts[i] + block_sizes[i]
            intron_end   = int(entry.start) + block_starts[i + 1]
            if intron_end - intron_start > 0:
                intron_lines.append(
                    f"{entry.chrom}\t{intron_start}\t{intron_end}\t{entry.name}\t.\t{entry.strand}"
                )
    return pybedtools.BedTool(intron_lines)


def extract_monoexonic_references(bed):
    """
    Identify mono-exonic reference transcripts from a BED12 file.
    Returns a BedTool object of exon intervals corresponding to mono-exonic transcripts.
    """
    mono_ref_lines = []
    for f in bed:
        block_count = int(f[9])
        if block_count == 1:
            start = f.start + int(f.fields[11].split(',')[0])
            end = start + int(f.fields[10].split(',')[0])
            mono_ref_lines.append(f"{f.chrom}\t{start}\t{end}\t{f.name}\t.\t{f.strand}")
    logging.info(f"Identified {len(mono_ref_lines)} mono-exonic reference transcripts.")
    return pybedtools.BedTool(mono_ref_lines)


def filter_premrna_like_monoexons(isoform_bed, reference_bed, min_intron_overlap):
    """
    Filter mono-exonic isoforms that overlap both reference exons and introns,
    but do not overlap known mono-exonic reference transcripts.

    Parameters:
        isoform_bed (BedTool): BED12 of predicted isoforms
        reference_bed (BedTool): BED12 reference annotations
        min_intron_overlap (int): Minimum overlap with introns to flag

    Returns:
        set: Transcript IDs flagged as likely pre-mRNA
    """
    reference_exons = reference_bed
    reference_introns = extract_introns_from_bed12(reference_bed)
    mono_ref_bed = extract_monoexonic_references(reference_bed)

    mono_lines = []
    mono_ids = set()
    for f in isoform_bed:
        block_count = int(f[9])
        if block_count == 1:
            start = f.start + int(f.fields[11].split(',')[0])
            end = start + int(f.fields[10].split(',')[0])
            mono_ids.add(f.name)
            mono_lines.append(f"{f.chrom}\t{start}\t{end}\t{f.name}\t.\t{f.strand}")
    mono_bed = pybedtools.BedTool(mono_lines)
    logging.info(f"Identified {len(mono_ids)} mono-exonic isoforms.")

    exon_hits = mono_bed.intersect(reference_exons, u=True, s=True)
    exon_hit_ids = {line.name for line in exon_hits}
    logging.info(f"Mono-exons overlapping reference exons: {len(exon_hit_ids)}")

    intron_hit_ids = set()
    if len(reference_introns) > 0:
        intron_hits = mono_bed.intersect(reference_introns, wo=True, s=True)
        for feature in intron_hits:
            overlap_bp = int(feature[-1])
            if overlap_bp >= min_intron_overlap:
                intron_hit_ids.add(feature.name)
        logging.info(f"Mono-exons overlapping ≥{min_intron_overlap} bp of introns: {len(intron_hit_ids)}")

    mono_ref_overlap_ids = {f.name for f in mono_bed.intersect(mono_ref_bed, u=True, s=True)}
    logging.info(f"Mono-exons overlapping mono-exonic reference transcripts: {len(mono_ref_overlap_ids)}")

    flagged_ids = (exon_hit_ids & intron_hit_ids) - mono_ref_overlap_ids
    logging.info(f"Flagged as likely pre-mRNA (after mono-exon ref exclusion): {len(flagged_ids)}")

    return flagged_ids


########
# MAIN #
########

if __name__ == '__main__':
    args = parse_args()

    logging.info("Filtering mono-exonic isoforms for potential pre-mRNA fragments.")

    isoform_bed = pybedtools.BedTool(args.isoform_bed12)
    reference_bed = pybedtools.BedTool(args.reference_bed12)

    flagged_ids = filter_premrna_like_monoexons(
        isoform_bed=isoform_bed,
        reference_bed=reference_bed,
        min_intron_overlap=args.min_intron_overlap
    )

    logging.info(f"Writing {len(flagged_ids)} transcript IDs to {args.out_ids}")
    with open(args.out_ids, 'w') as out_id_file:
        for tid in sorted(flagged_ids):
            out_id_file.write(tid + '\n')

    logging.info(f"Writing filtered BED12 to {args.out_bed}")
    with open(args.out_bed, 'w') as outfile:
        for f in isoform_bed:
            if f.name in flagged_ids:
                outfile.write(str(f))

    gc.collect()
    logging.info("Processing complete.")
