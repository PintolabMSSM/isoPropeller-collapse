#!/usr/bin/env python
# -*- coding: utf-8 -*-

##################
# IMPORT MODULES #
##################

import argparse
import gzip
import logging
import sys
import gc

logging.basicConfig(
    level   = logging.INFO,
    format  = '%(asctime)s - %(levelname)s - %(message)s',
    datefmt = '%Y-%m-%d %H:%M:%S'
)

PATSEQLEN = 10  # core length of sequence flanks (excluding wiggle)

##########################################################################
# CHANGE LOG (v2)
#
# 1. CRITICAL correctness bug fixed: the "+"-strand ex_start formula (and,
#    symmetrically, the "-"-strand in_start formula) had a spurious extra
#    "-1" that shifted the donor-based search window 1bp away from the
#    splice junction. Confirmed via a strand-symmetric "mirror" test
#    (reverse-complementing an entire synthetic chromosome to build a
#    biologically IDENTICAL "-"-strand transcript from a "+"-strand one):
#    the unmodified original MISSED the textbook "zero-gap" RT-switching
#    signature (an exact min_match-bp repeat immediately flanking the
#    junction, no gap) on the "+" strand at every tested --wiggle value
#    (0-3), while correctly detecting the identical biological case on the
#    "-" strand. Fixed by removing the extra "-1" from both formulas.
# 2. CRITICAL correctness bug fixed: check_for_repeat_pat() capped its
#    reported match length at exactly --min_match whenever the natural
#    exact-match extension already reached or exceeded that length on its
#    own, discarding the TRUE (possibly much longer) extension length. This
#    defeated the function's own boundary-aware "overhang" filter in
#    detect_rt_switching() (whose job is to reject a match whose real extent
#    reaches into the "wiggle" padding zone, since that signals an
#    ambiguous/shiftable junction rather than a clean RTS signature):
#    confirmed directly with a hand-built case where the true exact-match
#    run extended 2bp past core_right, yet the capped matchLen made the
#    overhang filter wrongly ACCEPT it. Fixed by reporting the true
#    extension length (not the min_match-capped one) whenever that
#    extension alone already satisfies min_match.
# 3. Performance: genome loading replaced Bio.SeqIO (a Bio.SeqRecord/Bio.Seq
#    object per chromosome) with a small dependency-free FASTA parser into a
#    plain dict of chrom -> uppercase str, and per-junction sequence
#    extraction/reverse-complementing replaced Bio.Seq slicing/
#    .reverse_complement() with plain string slicing and a
#    str.translate()-based IUPAC-aware reverse-complement table. The `Bio`
#    (Biopython) dependency is dropped entirely.
# 4. Performance: BED12 parsing replaced pybedtools' per-line Interval
#    object construction with a single plain-text pass (smart_open-wrapped
#    for gzip support). Since no real `bedtools` operation is ever used in
#    this script, the `pybedtools` dependency is dropped entirely.
# 5. Bonus capability (not a regression fix): the original's Bio.SeqIO
#    genome loader could not read a gzip-compressed (.fa.gz) genome FASTA
#    directly. The new loader is routed through smart_open(), so a
#    gzip-compressed --genome_fasta now works where it previously required
#    manual decompression first.
#
# All CLI flags, output file formats (--out_rts_tsv, --out_ids, --out_bed),
# and the RT-switching detection logic/semantics (including the
# seed-and-extend matching algorithm, the --min_match/--wiggle/
# --allow_mismatch behavior, and the duplicate-free isoform-file-order
# convention for --out_bed) are otherwise unchanged -- only the two bugs
# above were corrected, and the implementation was optimized for speed.
##########################################################################


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


def smart_open(path, mode="rt"):
    """Open a file transparently, decompressing if it ends in .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# IUPAC-ambiguity-aware complement table, built with str.translate() for
# speed, matching Bio.Seq.reverse_complement()'s supported alphabet
# (uppercase and lowercase nucleotide codes).
_COMPLEMENT_PAIRS = {
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
_COMPLEMENT_MAP = {}
for _b, _c in _COMPLEMENT_PAIRS.items():
    _COMPLEMENT_MAP[_b] = _c
    _COMPLEMENT_MAP[_b.lower()] = _c.lower()
_COMPLEMENT_TABLE = str.maketrans(_COMPLEMENT_MAP)


def reverse_complement(s):
    return s.translate(_COMPLEMENT_TABLE)[::-1]


def load_genome(fasta_path):
    """Load the reference genome FASTA into a dict of chrom -> uppercase str.

    A small dependency-free parser, in place of Bio.SeqIO -- gzip-aware via
    smart_open() (a new capability vs. the original, which could not read a
    .fa.gz genome directly).
    """
    logging.info(f"Loading genome FASTA: {fasta_path}")
    # A per-line streaming parse (rather than reading the whole file into one
    # string and splitting it) was deliberately kept here: a bulk
    # read-then-split approach was measured faster for CPU time but held
    # several large intermediate copies of the whole genome in memory at
    # once, roughly doubling peak RSS for a large genome -- a bad trade for
    # a ~0.2s difference in load time.
    genome = {}
    chrom = None
    chunks = []
    with smart_open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if chrom is not None:
                    genome[chrom] = "".join(chunks).upper()
                chrom = line[1:].split()[0].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if chrom is not None:
            genome[chrom] = "".join(chunks).upper()
    return genome


def _iter_bed12_lines(path):
    """Yield (raw_line, fields) for each BED12 line with >=12 tab fields."""
    with smart_open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            yield line, fields


def get_junctions_from_fields(fields):
    """Given BED12 fields, return list of introns as tuples:
    (chrom, donor_site, acceptor_site, strand, junction_number)
    """
    chrom = fields[0]
    start = int(fields[1])
    strand = fields[5]
    block_sizes = list(map(int, fields[10].strip(',').split(',')))
    block_starts = list(map(int, fields[11].strip(',').split(',')))

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
                    # BUGFIX: when the natural exact extension (k) already
                    # meets/exceeds min_match, m<=0 and the TRUE matched run
                    # is k itself, not the min_match-capped k+m -- report the
                    # true length so the boundary-aware overhang filter
                    # cannot be fooled by an under-reported match extent.
                    true_len = k if m <= 0 else k + m
                    return True, true_len, seq_exon[i:i + true_len], mismatch, i, j
            if i - m >= 0 and j - m >= 0:
                flag, mismatch = seq_match(seq_exon[i - m:i], seq_intron[j - m:j], allow_mismatch)
                if flag:
                    # Include start indices for left-extension case
                    return True, k + m, seq_exon[i - m:i + k], mismatch, i - m, j - m
    # Extend the return tuple with None placeholders for indices
    return False, None, None, None, None, None


def detect_rt_switching(chrom, name, fields, genome, min_match, wiggle, allow_mismatch):
    """Perform RT switching detection for all junctions in a BED12 entry.

    Args:
        chrom (str): chromosome
        name (str): transcript ID
        fields (list): BED12 fields for this entry
        genome (dict): Reference genome sequences (chrom -> uppercase str)
        min_match (int): Minimum match length
        wiggle (int): Bases allowed to shift search region
        allow_mismatch (bool): Allow one mismatch if True

    Returns:
        list of dict: Each dict contains info for one detected RT junction
    """
    if chrom not in genome:
        return []
    results = []
    junctions = get_junctions_from_fields(fields)
    seq = genome[chrom]

    for (chrom, donor, acceptor, strand, jnum) in junctions:
        cnt = PATSEQLEN + 2 * wiggle
        if strand == "+":
            # BUGFIX: removed a spurious extra "-1" that shifted this window
            # 1bp away from the junction (donor is already the first
            # intronic base in 0-based half-open BED coordinates, so
            # seq[donor-cnt+wiggle : donor+wiggle] already ends exactly at
            # the junction for wiggle=0, with no further adjustment needed).
            ex_start = donor - cnt + wiggle
            in_start = acceptor - cnt + wiggle
            # Plain string slicing mirrors Bio.Seq's slicing semantics
            # exactly (truncates at the sequence end, and -- since donor,
            # acceptor, and wiggle are all always >=0, a negative start here
            # can only arise together with a non-negative stop, which
            # Python's slicing always resolves to an empty string, not a
            # wraparound read from the far end of the chromosome -- the same
            # structurally-safe behavior the original relied on).
            seq_exon = seq[ex_start:ex_start + cnt]
            seq_intron = seq[in_start:in_start + cnt]
        else:
            # BUGFIX: same spurious extra "-1" removed (mirror of the "+"
            # strand fix above -- this is the donor-based window on the "-"
            # strand side).
            in_start = donor - wiggle
            ex_start = acceptor - wiggle
            seq_intron = reverse_complement(seq[in_start:in_start + cnt])
            seq_exon = reverse_complement(seq[ex_start:ex_start + cnt])
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
                    'isoform': name,
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
                    f"Filtered overhang match at {name} j{jnum}: "
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
    flagged_lines = []  # (raw_line) in isoform-file order, matching original convention

    fields_cols = ['isoform', 'junction_number', 'chrom', 'strand', 'donor_pos', 'acceptor_pos', 'match_len', 'match_pat', 'mismatch']
    with open(args.out_rts_tsv, "w") as out_f:
        out_f.write("\t".join(fields_cols) + "\n")
        for raw_line, bfields in _iter_bed12_lines(args.isoform_bed12):
            chrom = bfields[0]
            name = bfields[3]
            results = detect_rt_switching(chrom, name, bfields, genome, args.min_match, args.wiggle, args.allow_mismatch)
            if results:
                flagged_ids.add(name)
                flagged_lines.append(raw_line)
                for rec in results:
                    out_f.write("\t".join(str(rec[f]) for f in fields_cols) + "\n")

    # Write isoform IDs with detected RT switching
    with open(args.out_ids, "w") as id_out:
        for tid in sorted(flagged_ids):
            id_out.write(tid + "\n")

    # Write BED12 entries of isoforms with RT switching (isoform-file order,
    # matching the original's convention)
    with open(args.out_bed, "w") as bed_out:
        for raw_line in flagged_lines:
            bed_out.write(raw_line + "\n")

    gc.collect()
    logging.info("RT switching detection complete.")
