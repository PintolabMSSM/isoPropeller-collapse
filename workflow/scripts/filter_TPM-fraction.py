#!/usr/bin/env python
# -*- coding: utf-8 -*-

##################
# IMPORT MODULES #
##################

import argparse
import logging
import csv
from pybedtools import BedTool

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

#############
# FUNCTIONS #
#############

def is_interactive():
    """Check if we're in an interactive session (e.g. Jupyter notebook)."""
    import __main__ as main
    return not hasattr(main, '__file__')


def parse_args():
    parser = argparse.ArgumentParser(description="Filter isoforms based on expression across samples using a CPM threshold.")
    parser.add_argument("--count_matrix", required=True, help="TSV file of raw counts. First column is isoform ID; header row required.")
    parser.add_argument("--isoform_bed12", required=True, help="BED12 file of isoforms (ID in 4th column).")
    parser.add_argument("--min_fraction_samples", type=float, default=0.05, help="Minimum fraction of samples in which an isoform must be expressed (CPM >= threshold). Default: 0.05 (5%%).")
    parser.add_argument("--min_tpm", type=float, default=0.0, help="Minimum CPM to count a sample as 'expressed'. Default: 0.")
    parser.add_argument("--out_ids", required=True, help="Output file with IDs of filtered-out isoforms.")
    parser.add_argument("--out_bed", required=True, help="Output BED file of filtered-out isoforms.")
    
    if is_interactive():
        return parser.parse_args('')
    else:
        return parser.parse_args(args=None if sys.argv[1:] else ['--help'])


def compute_column_sums_and_count(filepath):
    sums = []
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        sample_count = len(header) - 1
        sums = [0] * sample_count
        for row in reader:
            for i in range(sample_count):
                sums[i] += int(row[i + 1])
    return sums, sample_count


def filter_low_expression_ids(filepath, col_sums, sample_count, min_fraction, min_tpm):
    low_expr_ids = set()
    threshold_count = int(min_fraction * sample_count)
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            tid = row[0]
            counts = list(map(int, row[1:]))
            cpm = [(1e6 * c / col_sums[i]) if col_sums[i] > 0 else 0 for i, c in enumerate(counts)]
            print(cpm)
            expressed = sum(1 for val in cpm if val >= min_tpm)
            if expressed < threshold_count:
                low_expr_ids.add(tid)
    return low_expr_ids


def write_filtered_bed(isoform_bed, filtered_ids, output_bed):
    bed = BedTool(isoform_bed)
    with open(output_bed, 'w') as out_bed:
        for entry in bed:
            if entry.name in filtered_ids:
                out_bed.write(str(entry))


########
# MAIN #
########

if __name__ == "__main__":
    args = parse_args()

    logging.info("Calculating column sums for CPM normalization...")
    col_sums, sample_count = compute_column_sums_and_count(args.count_matrix)

    print(col_sums)

    logging.info("Filtering isoforms based on expression threshold...")
    filtered_ids = filter_low_expression_ids(
        args.count_matrix, col_sums, sample_count,
        args.min_fraction_samples, args.min_tpm
    )

    logging.info(f"Identified {len(filtered_ids)} isoforms to remove.")
    with open(args.out_ids, 'w') as out_id_file:
        for tid in sorted(filtered_ids):
            out_id_file.write(tid + '\n')

    logging.info("Writing BED file of removed isoforms...")
    write_filtered_bed(args.isoform_bed12, filtered_ids, args.out_bed)

    logging.info("Expression filtering complete.")
