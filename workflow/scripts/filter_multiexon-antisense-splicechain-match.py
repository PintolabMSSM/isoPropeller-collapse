#!/usr/bin/env python
# -*- coding: utf-8 -*-

##################
# IMPORT MODULES #
##################

# Change log:
#   v2  2026-10  Correctness fixes + speed optimization.
#
#     1) CRITICAL BUG, confirmed to make the original script a complete
#        no-op: when building `antisense_chain_set`, the original computed,
#        for each reference transcript, `strand` (that transcript's own
#        strand) and `antisense` (the opposite strand), then filtered that
#        SAME transcript's own junction list with
#        `set((c, d, a) for (c, d, a, s) in sjs if s == antisense)`.
#        Every junction tuple `(c, d, a, s)` for a given reference
#        transcript has `s == strand` (its own strand, the same for every
#        junction of that one transcript) -- never `s == antisense`, by
#        construction (`antisense` is defined as the opposite of `strand`).
#        So this filter condition can never be true for any junction of any
#        reference transcript, `antisense_chain_set` is unconditionally
#        empty after the loop regardless of input, `sj_set & antisense_
#        chain_set` is therefore always the empty set for every isoform,
#        and `removed_ids` is always empty. **Confirmed directly**: a
#        hand-built test with two isoforms whose splice junctions exactly
#        match a reference transcript's junction coordinates but on the
#        opposite strand (the textbook case this script exists to catch)
#        -- the unmodified original reports "Removed 0 isoforms" and
#        writes empty `--out_bed`/`--out_ids` files, every time, regardless
#        of how much genuine antisense overlap is present in the input.
#        Fixed by keying the accumulated reference junction coordinates by
#        the STRAND they would be antisense *to* (i.e. `antisense_by_
#        strand["-"]` holds every "+"-strand reference junction's
#        coordinates, and vice versa), and, when checking an isoform,
#        looking up its own strand in that dict -- `sj_set & antisense_by_
#        strand.get(isoform_strand, set())` -- instead of comparing against
#        one single (and, in the original, unconditionally empty) flat set.
#        This also required reading each isoform's own strand (available
#        in `isoform_junctions[tid]`, which the original already builds but
#        never consulted for strand) rather than ignoring it.
#
#     2) SECOND BUG (dormant in the original only because bug #1 means
#        `removed_ids` is always empty, so this line is never actually
#        reached): `out_ids.write(tid + "\\n")` -- note the DOUBLE
#        backslash in the source. In a Python string literal, `"\\n"` is
#        an escaped backslash (`\`) followed by the literal letter `n`,
#        i.e. the 2-character string `\n` -- NOT a newline character
#        (which needs a single backslash in source, `"\n"`). Had bug #1
#        not made `removed_ids` always empty, every removed transcript ID
#        would have been written to `--out_ids` back-to-back on a single
#        physical line with literal `\n` text between them instead of real
#        line breaks, breaking any downstream tool expecting one ID per
#        line (`wc -l`, `grep -f`, etc). **Confirmed directly** with an
#        isolated reproduction of the same string literal. Fixed to a
#        proper single-backslash `"\n"`, matching the convention used
#        everywhere else in this pipeline.
#
#     3) Performance: `extract_splice_junctions()` iterated both the
#        isoform and reference BED12 files through pybedtools' per-line
#        Interval object construction merely to pull out a handful of
#        tab-separated fields (the same anti-pattern identified and fixed
#        in the sibling `filter_monoexon-premrna-fragments.py` and
#        `filter_monoexon-tss-overlap.py` scripts), and `main()` then
#        iterated the isoform file a SECOND time, via pybedtools Intervals
#        again, just to decide what to write out -- on top of the full
#        pass already done to build `isoform_junctions`. This script never
#        actually calls a real `bedtools` subcommand (no `.intersect()`,
#        `.slop()`, `.merge()`, etc. anywhere) -- pybedtools was being used
#        purely, and expensively, as a BED12 field parser. Rewritten as a
#        single plain-text pass per file (reference file first, to build
#        the by-strand antisense coordinate sets; then one pass over the
#        isoform file that extracts junctions, decides removal, and writes
#        both outputs for each isoform as it goes, in the same single
#        pass and in the isoform file's original row order, matching the
#        original's own non-sorted output-order convention for this
#        script). The `pybedtools` import/dependency is no longer needed
#        at all.
#
#     Caught during this rewrite, before it shipped: switching
#     `--isoform_bed12`/`--reference_bed12` from pybedtools' native
#     Interval-based iteration to plain `open()` (fix #3 above) silently
#     broke gzip support for both inputs, which the original supported
#     transparently via pybedtools (confirmed directly: the unmodified
#     original handles `.gz` isoform/reference files with identical -- if
#     always-empty -- output to plain-text ones; the initial plain-`open()`
#     draft of this rewrite raised `UnicodeDecodeError: 'utf-8' codec can't
#     decode byte 0x8b...` on the same files). Fixed by adding a
#     `smart_open()` helper (matching the same pattern used in the sibling
#     `consolidate_splice_chain_fragments.py`/
#     `prune_rare_isoforms_from_clusters.py`/
#     `filter_monoexon-premrna-fragments.py`/
#     `filter_monoexon-tss-overlap.py` fixes) and routing both BED12 reads
#     through it.
#
#     All CLI flags and output file formats (one BED12 line per removed
#     isoform in `--out_bed`, one transcript ID per line in `--out_ids`,
#     written in isoform-file order) are unchanged.
#
##################

import argparse
import gzip
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


def smart_open(path):
    """
    Open a text file for reading, transparently decompressing it if its name
    ends in `.gz`. Plain `open()` does not auto-detect gzip, unlike
    pybedtools' native Interval-based iteration -- used here so that
    `--isoform_bed12`/`--reference_bed12` still support gzip-compressed
    inputs the way the original script did via pybedtools.
    """
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def _iter_bed12_fields(path):
    """
    Yield each BED12 line's fields as a plain list of strings, parsed
    directly with str.split() rather than through pybedtools' per-line
    Interval object construction (this script never calls a real
    `bedtools` subcommand, so pybedtools was only ever being used here as a
    slow BED12 field parser).
    """
    with smart_open(path) as fh:
        for line in fh:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            yield fields, line


def _junctions_from_fields(fields):
    """
    Given one BED12 entry's fields, return (chrom, strand, start, junctions)
    where `junctions` is a list of (donor, acceptor) genomic-coordinate
    tuples, one per splice junction (empty for a mono-exonic/single-block
    entry).
    """
    chrom = fields[0]
    start = int(fields[1])
    strand = fields[5]
    block_sizes = list(map(int, fields[10].strip(',').split(',')))
    block_starts = list(map(int, fields[11].strip(',').split(',')))

    junctions = []
    for i in range(len(block_starts) - 1):
        donor = start + block_starts[i] + block_sizes[i]
        acceptor = start + block_starts[i + 1]
        junctions.append((donor, acceptor))
    return chrom, strand, junctions


def build_antisense_chain_sets(reference_bed12_path):
    """
    Scan the reference BED12 file once and build a dict keyed by the
    strand a reference junction would be ANTISENSE *to*:
        antisense_by_strand["+"] -> coordinates of every "-"-strand
                                     reference junction
        antisense_by_strand["-"] -> coordinates of every "+"-strand
                                     reference junction
    so that an isoform on strand X can be checked against `antisense_by_
    strand[X]` to test whether any of its own junctions coincide with a
    reference junction annotated on the opposite strand.

    Returns:
        dict: {"+": set of (chrom, donor, acceptor), "-": set of (...)}
    """
    antisense_by_strand = {"+": set(), "-": set()}
    opposite = {"+": "-", "-": "+"}

    for fields, _ in _iter_bed12_fields(reference_bed12_path):
        chrom, strand, junctions = _junctions_from_fields(fields)
        if not junctions or strand not in opposite:
            continue
        target_strand = opposite[strand]
        antisense_by_strand[target_strand].update(
            (chrom, donor, acceptor) for donor, acceptor in junctions
        )
    return antisense_by_strand


########
# MAIN #
########

if __name__ == "__main__":
    args = parse_args()

    logging.info("Extracting splice junctions...")
    antisense_by_strand = build_antisense_chain_sets(args.reference_bed12)

    removed_ids = set()
    with open(args.out_bed, "w") as out_bed, open(args.out_ids, "w") as out_ids:
        for fields, line in _iter_bed12_fields(args.isoform_bed12):
            tid = fields[3]
            chrom, strand, junctions = _junctions_from_fields(fields)
            if not junctions:
                continue  # skip isoforms without junctions (mono-exonic)

            sj_set = {(chrom, donor, acceptor) for donor, acceptor in junctions}
            antisense_targets = antisense_by_strand.get(strand, set())

            if sj_set & antisense_targets:
                out_bed.write(line if line.endswith("\n") else line + "\n")
                removed_ids.add(tid)
                out_ids.write(tid + "\n")

    logging.info(f"Removed {len(removed_ids)} isoforms with antisense splice chain overlap.")
