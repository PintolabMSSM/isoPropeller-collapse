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
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

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
    parser = argparse.ArgumentParser(description="Identify isoforms with non-canonical splice junctions from BED12.")
    parser.add_argument("--isoform_bed12",  required=True, help="BED12 file of isoforms")
    parser.add_argument("--genome_fasta",   required=True, help="Reference genome FASTA file")
    parser.add_argument("--out_ids",        required=True, help="Output file for isoform IDs with non-canonical junctions")
    parser.add_argument("--out_motifs",     required=True, help="Output TSV listing non-canonical motifs per isoform")
    parser.add_argument("--out_bed",      required=True, help="Output BED12 file of isoforms with non-canonical junctions")
    parser.add_argument("--allowed_motifs", nargs='*', default=["GTAG", "GCAG", "ATAC"],
                        help="Allowed splice site motifs (default: GTAG GCAG ATAC)")

    if is_interactive():
        args = parser.parse_args('')
    else:
        args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    return args

def load_genome(fasta_path):
    """Load genome FASTA into dictionary of chrom -> sequence."""
    logging.info(f"Loading genome FASTA: {fasta_path}")
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    return genome


def get_junctions_from_bed12(bed_entry):
    """
    Given a BED12 entry, return list of introns as tuples:
    (chrom, donor_site, acceptor_site, strand)
    """
    chrom = bed_entry.chrom
    start = int(bed_entry.start)
    strand = bed_entry.strand
    block_sizes = list(map(int, bed_entry.fields[10].strip(",").split(",")))
    block_starts = list(map(int, bed_entry.fields[11].strip(",").split(",")))

    junctions = []
    for i in range(len(block_starts) - 1):
        exon_end = start + block_starts[i] + block_sizes[i]
        next_exon_start = start + block_starts[i + 1]
        junctions.append((chrom, exon_end, next_exon_start, strand))
    return junctions


def fetch_splice_motif(chrom, donor_pos, acceptor_pos, strand, genome):
    """
    Fetch canonical splice site motif (2bp donor + 2bp acceptor), strand-aware.
    """
    if chrom not in genome:
        return "NNNN"

    seq = genome[chrom].seq
    try:
        donor = seq[donor_pos:donor_pos+2]
        acceptor = seq[acceptor_pos-2:acceptor_pos]
        motif = donor + acceptor
        if strand == "-":
            motif = motif.reverse_complement()
        return str(motif).upper()
    except Exception:
        return "NNNN"


########
# MAIN #
########

if __name__ == '__main__':
    args = parse_args()
    genome = load_genome(args.genome_fasta)
    bed = pybedtools.BedTool(args.isoform_bed12)

    allowed_motifs = set(args.allowed_motifs)
    logging.info(f"Allowed splice site motifs: {allowed_motifs}")

    noncanonical_ids = set()
    noncanonical_motifs_per_id = defaultdict(set)

    for entry in bed:
        if len(entry.fields) < 12:
            continue
        tid = entry.name
        strand = entry.strand
        junctions = get_junctions_from_bed12(entry)

        for (chrom, donor_pos, acceptor_pos, strand) in junctions:
            motif = fetch_splice_motif(chrom, donor_pos, acceptor_pos, strand, genome)
            if motif not in allowed_motifs:
                noncanonical_ids.add(tid)
                noncanonical_motifs_per_id[tid].add(motif)

    logging.info(f"Found {len(noncanonical_ids)} isoforms with non-canonical splice sites.")

    # Write list of noncanonical transcript IDs
    with open(args.out_ids, 'w') as out_f:
        for tid in sorted(noncanonical_ids):
            out_f.write(tid + "\n")

    # Write TSV of motifs per transcript
    with open(args.out_motifs, 'w') as outf:
        for tid in sorted(noncanonical_motifs_per_id):
            motifs = ",".join(sorted(noncanonical_motifs_per_id[tid]))
            outf.write(f"{tid}\t{motifs}\n")

    # Write filtered BED12 of noncanonical isoforms
    logging.info(f"Writing BED12 entries for non-canonical isoforms to {args.out_bed}")
    with open(args.out_bed, 'w') as outbed:
        for entry in bed:
            if entry.name in noncanonical_ids:
                outbed.write(str(entry))

    gc.collect()
    logging.info("Processing complete.")
