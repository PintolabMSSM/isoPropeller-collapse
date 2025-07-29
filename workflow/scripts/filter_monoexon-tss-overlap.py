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
import gc

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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
        description="Filter mono-exonic BED12 isoforms whose TSS is not within a given distance of any reference TSS peak."
    )
    parser.add_argument("--isoform_bed12",     required=True,  help="BED12 file of isoforms")
    parser.add_argument("--isoform_tss_bed",   required=True,  help="BED6 file of isoform TSSs (column 4 = transcript ID)")
    parser.add_argument("--reftss_bed",        required=True,  help="BED file of reference TSS peaks (e.g., CAGE)")
    parser.add_argument("--genome_index",      required=True,  help="FAI index file for genome (used for slop)")
    parser.add_argument("--out_bed",           required=True,  help="Output BED12 file of filtered isoforms")
    parser.add_argument("--out_ids",           required=True,  help="Transcript IDs that were removed")
    parser.add_argument("--max_distance",      type=int, default=10, help="Max distance from isoform TSS to reference TSS (default: 10 bp)")

    if is_interactive():
        return parser.parse_args('')
    else:
        return parser.parse_args(args=None if sys.argv[1:] else ['--help'])


def identify_monoexonic_ids(bed12_path):
    """
    Identify mono-exonic isoforms from a BED12 file.

    Args:
        bed12_path (str): Path to the input BED12 file.

    Returns:
        set: Transcript IDs (column 4) with only one exon (blockCount == 1).
    """
    mono_ids = set()
    with open(bed12_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            block_count = int(fields[9])
            if block_count == 1:
                tid = fields[3]
                mono_ids.add(tid)
    logging.info(f"Identified {len(mono_ids)} mono-exonic isoforms.")
    return mono_ids


def filter_tss_support(tss_bed_path, reftss_bed_path, genome_index, mono_ids, max_distance):
    """
    Filter mono-exonic isoforms whose TSS is not supported by proximity to a reference TSS.

    Args:
        tss_bed_path (str): Path to BED6 file of isoform TSSs (column 4 is transcript ID).
        reftss_bed_path (str): Path to BED file of reference TSS peaks (e.g., CAGE).
        genome_index (str): Path to FAI file (used to slop reference TSSs).
        mono_ids (set): Transcript IDs identified as mono-exonic.
        max_distance (int): Maximum distance allowed between isoform TSS and reference TSS.

    Returns:
        set: Transcript IDs of mono-exonic isoforms removed due to insufficient TSS support.
    """
    all_tss = pybedtools.BedTool(tss_bed_path)

    # Get transcript IDs present in the TSS BED
    tss_ids = {f.name for f in all_tss}
    matched_ids = mono_ids & tss_ids
    missing_ids = mono_ids - tss_ids

    logging.info(f"Matched {len(matched_ids)} mono-exonic isoforms with their respective TSS ranges.")
    if missing_ids:
        logging.warning(f"{len(missing_ids)} mono-exonic isoforms were missing from the TSS BED and skipped from filtering.")

    # Filter to valid chromosomes and matched IDs
    with open(genome_index) as f:
        valid_chroms = {line.split('\t')[0] for line in f}

    tss_bed = pybedtools.BedTool([f for f in all_tss if f.name in matched_ids and f.chrom in valid_chroms])
    reftss = pybedtools.BedTool([f for f in pybedtools.BedTool(reftss_bed_path) if f.chrom in valid_chroms])
    reftss_slop = reftss.slop(b=max_distance, g=genome_index)

    # Compare without strand: find isoform TSS not overlapping refTSS
    unsupported = tss_bed.intersect(reftss_slop, v=True)
    removed_ids = {f.name for f in unsupported}

    logging.info(f"{len(removed_ids)} mono-exonic isoforms had TSS ranges but no nearby reference TSS.")
    return removed_ids


########
# MAIN #
########

if __name__ == '__main__':
    args = parse_args()

    logging.info("Filtering mono-exonic BED12 isoforms based on TSS proximity to reference TSS peaks.")

    mono_ids = identify_monoexonic_ids(args.isoform_bed12)
    removed_ids = filter_tss_support(
        tss_bed_path     = args.isoform_tss_bed,
        reftss_bed_path  = args.reftss_bed,
        genome_index     = args.genome_index,
        mono_ids         = mono_ids,
        max_distance     = args.max_distance
    )

    logging.info(f"Writing {len(removed_ids)} removed transcript IDs to {args.out_ids}")
    with open(args.out_ids, 'w') as out_f:
        for tid in sorted(removed_ids):
            out_f.write(tid + '\n')

    logging.info(f"Writing filtered BED12 to {args.out_bed}")
    with open(args.isoform_bed12) as infile, open(args.out_bed, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) >= 4 and fields[3] in removed_ids:
                outfile.write(line)

    gc.collect()
    logging.info("Processing complete.")
