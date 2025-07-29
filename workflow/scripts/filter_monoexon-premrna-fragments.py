#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  monoexon-premrna-filter
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
    import __main__ as main
    return not hasattr(main, '__file__')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove mono-exonic isoforms that likely represent pre-mRNA fragments (overlap reference exon AND intron)."
    )
    parser.add_argument("--isoform_gtf",         required=True,  help="GTF file of isoform predictions")
    parser.add_argument("--reference",           required=True,  help="Reference annotations in GTF or BED format")
    parser.add_argument("--genome_index",        required=True,  help="FAI index file for genome (used for intron handling)")
    parser.add_argument("--out_gtf",             required=True,  help="Output GTF of filtered isoforms")
    parser.add_argument("--out_ids",             required=True,  help="Transcript IDs flagged as likely pre-mRNA")
    parser.add_argument("--min_intron_overlap",  type=int, default=10, help="Minimum base pairs of overlap with a reference intron to flag (default: 10)")
    
    if is_interactive():
        return parser.parse_args('')
    else:
        return parser.parse_args(args=None if sys.argv[1:] else ['--help'])


def extract_transcript_id(attribute_field):
    match = re.search(r'transcript_id\s+"([^"]+)"', attribute_field)
    return match.group(1) if match else None


def extract_introns_from_gtf(reference_gtf):
    exons_by_transcript = defaultdict(list)
    for f in reference_gtf:
        if f[2] == 'exon':
            tid = extract_transcript_id(f[8])
            if tid:
                exons_by_transcript[tid].append(f)

    intron_lines = []
    for tid, exons in exons_by_transcript.items():
        if len(exons) < 2:
            continue
        sorted_exons = sorted(exons, key=lambda x: x.start)
        for i in range(len(sorted_exons) - 1):
            intron_start = sorted_exons[i].end
            intron_end = sorted_exons[i+1].start
            if intron_end - intron_start > 0:
                intron_lines.append(
                    f"{sorted_exons[i].chrom}\t{intron_start}\t{intron_end}\t{tid}\t.\t{sorted_exons[i].strand}"
                )
    return pybedtools.BedTool(intron_lines)


def extract_introns_from_bed12(bed):
    """Extract introns from BED12 format (blockStarts and blockSizes)."""
    intron_lines = []
    for entry in bed:
        if len(entry.fields) < 12:
            continue
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


def load_reference(reference_path):
    """Load reference exons and introns depending on file extension and format."""
    _, ext = os.path.splitext(reference_path.lower())
    bed = pybedtools.BedTool(reference_path)

    if ext in [".gtf", ".gff"]:
        exons = bed.filter(lambda f: f[2] == 'exon').saveas()
        introns = extract_introns_from_gtf(bed)
        logging.info(f"Loaded GTF with {len(exons)} exons and {len(introns)} derived introns.")
        return exons, introns

    elif ext in [".bed", ".bed6", ".bed12"]:
        exons = bed.saveas()
        has_blocks = all(len(f.fields) >= 12 for f in bed[:10])  # check first few lines
        if has_blocks:
            introns = extract_introns_from_bed12(bed)
            logging.info(f"Loaded BED12 with {len(exons)} exons and {len(introns)} derived introns.")
        else:
            introns = pybedtools.BedTool([])
            logging.warning("BED file lacks block structure (not BED12); intron overlap check skipped.")
        return exons, introns

    else:
        raise ValueError(f"Unsupported reference file extension: {ext}")


def filter_premrna_like_monoexons(isoform_gtf_path, reference_path, genome_index, min_intron_overlap):
    isoform_bed = pybedtools.BedTool(isoform_gtf_path)
    reference_exons, reference_introns = load_reference(reference_path)

    # Group isoform exons by transcript
    transcript_exons = defaultdict(list)
    for feature in isoform_bed:
        if feature[2] == 'exon':
            tid = extract_transcript_id(feature[8])
            if tid:
                transcript_exons[tid].append(feature)

    # Collect mono-exonic transcripts
    mono_lines = []
    mono_ids = set()
    for tid, exons in transcript_exons.items():
        if len(exons) == 1:
            mono_ids.add(tid)
            exon = exons[0]
            mono_lines.append(f"{exon.chrom}\t{exon.start}\t{exon.end}\t{tid}\t.\t{exon.strand}")
    mono_bed = pybedtools.BedTool(mono_lines)
    logging.info(f"Identified {len(mono_ids)} mono-exonic transcripts.")

    # Intersect mono-exons with reference exons
    exon_hits = mono_bed.intersect(reference_exons, u=True, s=True)
    exon_hit_ids = {line.name for line in exon_hits}
    logging.info(f"Mono-exons overlapping reference exons: {len(exon_hit_ids)}")

    # Intersect mono-exons with introns (if available)
    intron_hit_ids = set()
    if len(reference_introns) > 0:
        intron_hits = mono_bed.intersect(reference_introns, wo=True, s=True)
        for feature in intron_hits:
            overlap_bp = int(feature[-1])
            if overlap_bp >= min_intron_overlap:
                intron_hit_ids.add(feature.name)
        logging.info(f"Mono-exons overlapping ≥{min_intron_overlap} bp of introns: {len(intron_hit_ids)}")
    else:
        logging.warning("No introns loaded — skipping intron overlap filtering.")

    flagged_ids = exon_hit_ids & intron_hit_ids if len(reference_introns) > 0 else set()
    logging.info(f"Flagged as likely pre-mRNA: {len(flagged_ids)} mono-exons.")

    return flagged_ids


########
# MAIN #
########

if __name__ == '__main__':
    args = parse_args()

    logging.info("Filtering mono-exonic isoforms for potential pre-mRNA fragments.")

    flagged_ids = filter_premrna_like_monoexons(
        isoform_gtf_path    = args.isoform_gtf,
        reference_path      = args.reference,
        genome_index        = args.genome_index,
        min_intron_overlap  = args.min_intron_overlap
    )

    logging.info(f"Writing {len(flagged_ids)} transcript IDs to {args.out_ids}")
    with open(args.out_ids, 'w') as out_id_file:
        for tid in sorted(flagged_ids):
            out_id_file.write(tid + '\n')

    logging.info(f"Writing filtered GTF to {args.out_gtf}")
    with open(args.isoform_gtf) as infile, open(args.out_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            tid = extract_transcript_id(fields[8])
            if tid not in flagged_ids:
                outfile.write(line)

    gc.collect()
    logging.info("Processing complete.")
