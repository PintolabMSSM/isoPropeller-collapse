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

logging.basicConfig(
    level   = logging.INFO,
    format  = '%(asctime)s - %(levelname)s - %(message)s',
    datefmt = '%Y-%m-%d %H:%M:%S'
)

PATSEQLEN = 10  # core length of sequence flanks (excluding wiggle)


#############
# FUNCTIONS #
#############

def is_interactive():
    """Check if we're in an interactive session (e.g. Jupyter notebook)."""
    import __main__ as main
    return not hasattr(main, '__file__')


def parse_args():
    """Parse command-line arguments for RT switching detection."""
    parser = argparse.ArgumentParser(description="Detect RT switching artifacts using SQANTI3 logic from BED12.")
    parser.add_argument("--isoform_bed12", required=True, help="BED12 file of isoforms")
    parser.add_argument("--genome_fasta", required=True, help="Reference genome FASTA file")
    parser.add_argument("--out_rts_tsv", required=True, help="Output TSV with detailed RTS match info")
    parser.add_argument("--out_ids", required=True, help="Output file of isoform IDs with RT switching")
    parser.add_argument("--out_bed", required=True, help="Output BED12 file of isoforms with RT switching")
    parser.add_argument("--min_match", type=int, default=8, choices=range(4, 11), help="Minimum match length (default: 8)")
    parser.add_argument("--wiggle", type=int, default=1, choices=range(0, 4), help="Wiggle bases allowed at junction (default: 1)")
    parser.add_argument("--allow_mismatch", action='store_true', help="Allow 1 mismatch in matching region")

    if is_interactive():
        args = parser.parse_args('')
    else:
        args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    return args
    

def load_genome(fasta_path):
    """Load the reference genome FASTA into a dictionary of chromosome -> SeqRecord."""
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))


def get_junctions_from_bed12(entry):
    """Extract junctions from a BED12 entry as a list of tuples: (chrom, donor, acceptor, strand, junction_number)."""
    chrom = entry.chrom
    start = int(entry.start)
    strand = entry.strand
    block_sizes = list(map(int, entry.fields[10].strip(',').split(',')))
    block_starts = list(map(int, entry.fields[11].strip(',').split(',')))

    junctions = []
    for i in range(len(block_starts) - 1):
        exon_end = start + block_starts[i] + block_sizes[i]
        next_exon_start = start + block_starts[i + 1]
        junctions.append((chrom, exon_end, next_exon_start, strand, i + 1))
    return junctions


def seq_match(exseq, inseq, allow_mismatch):
    """Compare two sequences and return whether they match exactly or with one mismatch.

    Args:
        exseq (str): Exon sequence
        inseq (str): Intron sequence
        allow_mismatch (bool): Whether to allow up to one mismatch

    Returns:
        (bool, int): Tuple indicating match status and mismatch count
    """
    if len(exseq) != len(inseq):
        return False, None
    if exseq == inseq:
        return True, 0
    if allow_mismatch:
        mismatches = sum(1 for a, b in zip(exseq, inseq) if a != b)
        return mismatches <= 1, mismatches
    return False, None


def check_for_repeat_pat(seq_exon, seq_intron, min_match, allow_mismatch):
    """Detect repeat motifs between exon and intron sequences indicating RT switching.

    Args:
        seq_exon (str): Sequence upstream of donor site
        seq_intron (str): Sequence downstream of acceptor site
        min_match (int): Minimum number of bases that must match
        allow_mismatch (bool): Whether to allow one mismatch in extended region

    Returns:
        # Return start indices to enable boundary-aware filtering:
        (bool, int, str, int, int, int)
        -> (match flag, length, pattern, mismatch count, exon_start_idx, intron_start_idx)
    """
    seedsize = min_match // 2
    n = len(seq_exon)
    for i in range(n - seedsize + 1):
        seed = seq_exon[i:i + seedsize]
        offset = 0
        while True:
            j = seq_intron.find(seed, offset)
            if j == -1:
                break
            offset = j + 1
            k = seedsize
            while i + k < n and j + k < len(seq_intron) and seq_exon[i + k] == seq_intron[j + k]:
                k += 1
            m = min_match - k
            if i + k + m <= n and j + k + m <= len(seq_intron):
                flag, mismatch = seq_match(seq_exon[i + k:i + k + m], seq_intron[j + k:j + k + m], allow_mismatch)
                if flag:
                    # Include start indices for right-extension case
                    return True, k + m, seq_exon[i:i + k + m], mismatch, i, j
            if i - m >= 0 and j - m >= 0:
                flag, mismatch = seq_match(seq_exon[i - m:i], seq_intron[j - m:j], allow_mismatch)
                if flag:
                    # Include start indices for left-extension case
                    return True, k + m, seq_exon[i - m:i + k], mismatch, i - m, j - m
    # Extend the return tuple with None placeholders for indices
    return False, None, None, None, None, None


def detect_rt_switching(entry, genome, min_match, wiggle, allow_mismatch):
    """Perform RT switching detection for all junctions in a BED12 entry.

    Args:
        entry: BED12 entry
        genome (dict): Reference genome sequences
        min_match (int): Minimum match length
        wiggle (int): Bases allowed to shift search region
        allow_mismatch (bool): Allow one mismatch if True

    Returns:
        list of dict: Each dict contains info for one detected RT junction
    """
    chrom = entry.chrom
    if chrom not in genome:
        return []
    results = []
    junctions = get_junctions_from_bed12(entry)
    seq = genome[chrom].seq

    for (chrom, donor, acceptor, strand, jnum) in junctions:
        cnt = PATSEQLEN + 2 * wiggle
        if strand == "+":
            ex_start = donor - cnt + wiggle - 1
            in_start = acceptor - cnt + wiggle
            seq_exon = str(seq[ex_start:ex_start + cnt]).upper()
            seq_intron = str(seq[in_start:in_start + cnt]).upper()
        else:
            in_start = donor - wiggle - 1
            ex_start = acceptor - wiggle
            seq_intron = str(seq[in_start:in_start + cnt].reverse_complement()).upper()
            seq_exon = str(seq[ex_start:ex_start + cnt].reverse_complement()).upper()
        if not seq_exon or not seq_intron:
            continue

        # Run matcher and capture start indices for boundary-aware filtering
        flag, matchLen, matchPat, mismatch, i_ex, j_in = check_for_repeat_pat(
            seq_exon, seq_intron, min_match, allow_mismatch
        )

        if flag:
            # Perform boundary-aware post-filtering:
            # Define the maximum core index that a match is allowed to reach (exclusive end index).
            # The last 'wiggle' bases on the right of each window are overhang and must not be used.
            core_right = len(seq_exon) - wiggle  # same length for both windows

            right_edge_ex = i_ex + matchLen   # exclusive end in exon window
            right_edge_in = j_in + matchLen   # exclusive end in intron window

            # Accept if both ends stay within the core (can end exactly at the core boundary).
            if right_edge_ex <= core_right and right_edge_in <= core_right:
                results.append({
                    'isoform': entry.name,
                    'junction_number': jnum,
                    'chrom': chrom,
                    'strand': strand,
                    'donor_pos': donor,
                    'acceptor_pos': acceptor,
                    'match_len': matchLen,
                    'match_pat': matchPat,
                    'mismatch': mismatch
                })
            else:
                # Provide verbose diagnostics on filtered hits
                logging.debug(
                    f"Filtered overhang match at {entry.name} j{jnum}: "
                    f"ex_right={right_edge_ex}, in_right={right_edge_in}, core_right={core_right}"
                )
                pass

    return results

########
# MAIN #
########

if __name__ == "__main__":
    args = parse_args()
    genome = load_genome(args.genome_fasta)
    flagged_ids = set()
    flagged_entries = []
    bed = pybedtools.BedTool(args.isoform_bed12)

    fields = ['isoform', 'junction_number', 'chrom', 'strand', 'donor_pos', 'acceptor_pos', 'match_len', 'match_pat', 'mismatch']
    with open(args.out_rts_tsv, "w") as out_f:
        out_f.write("\t".join(fields) + "\n")
        for entry in bed:
            if len(entry.fields) < 12:
                continue
            results = detect_rt_switching(entry, genome, args.min_match, args.wiggle, args.allow_mismatch)
            if results:
                flagged_ids.add(entry.name)
                flagged_entries.append(entry)
                for rec in results:
                    out_f.write("\t".join(str(rec[f]) for f in fields) + "\n")

    # Write isoform IDs with detected RT switching
    with open(args.out_ids, "w") as id_out:
        for tid in sorted(flagged_ids):
            id_out.write(tid + "\n")

    # Write BED12 entries of isoforms with RT switching
    with open(args.out_bed, "w") as bed_out:
        for entry in flagged_entries:
            bed_out.write(str(entry))

    gc.collect()
    logging.info("RT switching detection complete.")
