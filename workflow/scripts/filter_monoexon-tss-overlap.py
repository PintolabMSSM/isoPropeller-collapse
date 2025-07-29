#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# filter_monoexon-tss-overlap
#

##################
# IMPORT MODULES #
##################

import argparse
import pybedtools
import re
import sys
import os
import gc
import logging
from collections import defaultdict

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)  ## Suppress futurewarnings

# Configure logging to include the date and time in the log messages
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
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter mono-exonic isoforms whose 5' end is not within a certain distance of any reference TSS peak."
    )
    parser.add_argument("--isoform_gtf",     required=True,  help="Input GTF file of isoforms")
    parser.add_argument("--reftss_bed",      required=True,  help="BED file with reference TSS (e.g., CAGE peaks)")
    parser.add_argument("--out_gtf",         required=True,  help="Output GTF file of filtered isoforms")
    parser.add_argument("--out_ids",         required=True,  help="Output file of filtered transcript IDs")
    parser.add_argument("--max_distance",    type=int, default=100, help="Maximum distance from TSS to reference TSS peak (default: 100 bp)")
    parser.add_argument("--genome_index",    required=True,  help="FAI index file of the reference genome")

    if is_interactive():
        args = parser.parse_args('')
    else:
        args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    return args


def extract_transcript_id(attribute_field):
    """Extract transcript_id from GTF attribute field."""
    match = re.search(r'transcript_id\s+"([^"]+)"', attribute_field)
    return match.group(1) if match else None


def filter_monoexons_without_cage(isoforms_bed, tss_bed_path, max_distance, genome_index):
    """Return set of mono-exonic transcript IDs not supported by nearby reference TSS."""
    with open(genome_index) as f:
        valid_chroms = {line.split('\t')[0] for line in f}

    raw_tss = pybedtools.BedTool(tss_bed_path)
    filtered_tss = [f for f in raw_tss if f.chrom in valid_chroms]
    logging.info(f"Filtered TSS peaks to valid chromosomes: {len(filtered_tss)} retained.")

    tss_peaks = pybedtools.BedTool(filtered_tss).slop(b=max_distance, g=genome_index)

    transcript_exons = defaultdict(list)
    for feature in isoforms_bed:
        if feature[2] == 'exon':
            tid = extract_transcript_id(feature[8])
            if tid:
                transcript_exons[tid].append(feature)

    mono_tss_lines = []
    mono_ids = set()
    for tid, exons in transcript_exons.items():
        if len(exons) == 1:
            mono_ids.add(tid)
            exon = exons[0]
            tss = exon.start if exon.strand == '+' else exon.end - 1
            mono_tss_lines.append(f"{exon.chrom}\t{tss}\t{tss+1}\t{tid}")

    logging.info(f"Identified {len(mono_ids)} mono-exonic transcripts.")

    tss_bed = pybedtools.BedTool(mono_tss_lines)
    unsupported = tss_bed.intersect(tss_peaks, v=True)
    removed_ids = {line.name for line in unsupported}
    logging.info(f"{len(removed_ids)} transcripts have unsupported TSS (beyond {max_distance} bp).")

    return removed_ids


########
# MAIN #
########

if __name__ == '__main__':
    args = parse_args()

    logging.info("Loading isoform GTF as BED.")
    isoforms_bed = pybedtools.BedTool(args.isoform_gtf)

    removed_ids = filter_monoexons_without_cage(
        isoforms_bed,
        args.reftss_bed,
        args.max_distance,
        args.genome_index
    )

    logging.info(f"Writing {len(removed_ids)} filtered transcript IDs to {args.out_ids}.")
    with open(args.out_ids, "w") as tid_file:
        for tid in sorted(removed_ids):
            tid_file.write(tid + "\n")

    logging.info(f"Filtering GTF and writing output to {args.out_gtf}.")
    with open(args.isoform_gtf) as infile, open(args.out_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            tid = extract_transcript_id(fields[8])
            if tid in removed_ids:
                outfile.write(line)

    gc.collect()
    logging.info("Processing complete.")
