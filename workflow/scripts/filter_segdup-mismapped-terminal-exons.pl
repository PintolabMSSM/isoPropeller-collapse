#!/usr/bin/perl

# 02.01.2022 13:16:28 EST
# v2 fixes: see CHANGE LOG below.

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);

# GLOBALS
$ENV{TMPDIR}   ||= "/sc/arion/scratch/$ENV{USER}/tmp/";   # location for tmp file storage

# ============================================================================
# CHANGE LOG (v1 -> v2)
#
# Three correctness/robustness bugs were found and fixed, plus one
# documentation/error-message fix:
#
#   BUG 1 (coordinate off-by-one -> false positives in
#   fully_contained_in_single_segdup): the per-isoform bounding-box temp BED
#   file written for the "is this isoform fully contained in a single
#   segdup" check copied the GFF input's 1-based start coordinate directly
#   into a plain BED file's start column, with no -1 adjustment. bedtools
#   correctly auto-detects and 1-based/0-based-converts a genuine GFF/GTF
#   file when it is passed directly to `bedtools intersect` (confirmed
#   empirically), but this script instead builds a NEW, separate .bed file
#   from scratch for the bbox check, which bedtools reads as plain 0-based
#   BED with no conversion -- so every isoform's bbox was silently
#   positioned 1bp too far to the right (and 1bp too short), always missing
#   coverage of its own true leftmost base. This can only ever make
#   containment EASIER to satisfy than it truly is (the checked window is a
#   strict subset of the true isoform span), so it can only cause FALSE
#   POSITIVES: an isoform whose true leftmost base pokes exactly 1bp outside
#   a segdup's boundary (and is therefore NOT fully contained) gets wrongly
#   reported as "yes". Confirmed directly: a 2-exon isoform with true span
#   chr1:1-200 (1-based) against a segdup covering true bases 2-200 (i.e.
#   missing only the isoform's own first base) -- the original reports
#   fully_contained_in_single_segdup=yes; it should be "no". Fixed by
#   writing the bbox file with start-1 (the correct 0-based BED start).
#
#   BUG 2 (temp file leak): `tempfile(..., UNLINK=>0)` explicitly disables
#   File::Temp's automatic cleanup, and the file is never explicitly
#   unlinked anywhere -- so every single invocation of this script leaves
#   one orphaned isocoord-*.bed file behind in $TMPDIR forever. The
#   INTERRUPT subroutine's own comment ("make sure we exit gracefully for
#   tmp file cleanup") confirms cleanup was intended, but the subroutine
#   was never actually wired to a signal (no $SIG{INT}/$SIG{TERM}
#   assignment exists anywhere) and doesn't unlink anything even when
#   called. Confirmed directly: 4 repeated runs left 4 orphaned .bed files
#   in TMPDIR, accumulating indefinitely in a shared HPC scratch directory.
#   Fixed by switching to UNLINK=>1 (so File::Temp cleans it up on normal
#   exit), explicitly unlinking it as soon as it's no longer needed, and
#   actually wiring $SIG{INT}/$SIG{TERM} to the (now-functional) INTERRUPT
#   handler so an interrupted run also cleans up.
#
#   BUG 3 (missing CLI wiring): $nMinSegDupOverlap (0.90) and
#   $nMinGeneOverlap (1.00) are declared and used exactly like
#   $nTerminalLevel, but -- unlike $nTerminalLevel's "-level" option --
#   neither was ever added to the GetOptions() call, so there was never any
#   way to actually change them from the command line; passing e.g.
#   "--min_segdup_overlap 0.5" was silently ignored (GetOptions just warns
#   "Unknown option" on stderr and keeps going with the 0.90 default).
#   Confirmed directly. Fixed by adding --min_segdup_overlap/
#   --min_gene_overlap options wired to these variables, with the exact
#   same defaults (0.90/1.00) preserved when not given.
#
# Also fixed: the truncated/incomplete help text for -g/--genes ("File with
# gene regions in " with nothing after it), and the terminal-level error
# message, which said "must be an integer greater than 1" while the actual
# validation regex (and the default value itself, 1) accept 1.
#
# INVESTIGATED, NOT A BUG: an isoform that is fully contained within a
# single segdup but has too few exons/introns for the chosen --level (most
# notably a mono-exonic isoform, which has zero introns) never appears in
# the output at all, even though fully_contained_in_single_segdup is
# computed independently of the terminal-exon/intron checks. An earlier
# revision of this script "fixed" this by printing such isoforms anyway.
# That change has been REVERTED on the author's explicit confirmation that
# this tool is intentionally scoped to terminal exons/introns of spliced
# (multi-exonic) transcripts only -- a mono-exonic isoform has no terminal
# intron and is deliberately out of scope, independent of its segdup
# status. The original behavior (no row for such isoforms) is restored.
#
# All other logic -- the terminal-exon/intron level semantics, the
# mono-exon/multi-exon gene classification via the gene file's score column,
# the bedtools -f/-F overlap semantics for the exon/intron/segdup/gene
# checks, and the output column set/order -- is unchanged.
# ============================================================================

# GET PARAMETERS
my $sHelp                = 0;
my $sExInGffFile         = '';
my $sSegdupBedFile       = '';
my $sGeneregionBedFile   = '';
my $nTerminalLevel       = 1;
my $nMinSegDupOverlap    = 0.90;
my $nMinGeneOverlap      = 1.00;
GetOptions("help!"                => \$sHelp,
           "input:s"              => \$sExInGffFile,
           "segdup:s"             => \$sSegdupBedFile,
           "genes:s"              => \$sGeneregionBedFile,
           "level:i"              => \$nTerminalLevel,
           "min_segdup_overlap:f" => \$nMinSegDupOverlap,
           "min_gene_overlap:f"   => \$nMinGeneOverlap);

# PRINT HELP
$sHelp = 1 unless($sExInGffFile and $sSegdupBedFile and $sGeneregionBedFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName

   This script finds isoforms with terminal exons that are potentially mismapped due to
   segdups. This is done by looking for terminal exons (at the chosen level) that map
   within segdups and that span other genes contained within the terminal intron.

   Arguments:
    -i --input <string>
      Input isoform file in exon/intron gff format produced by 'bed2intronexongff.pl'
    -s --segdup <string>
      File with segdup regions in bed basic format.
    -g --genes <string>
      File with gene regions in the same exon/intron gff format as --input, with each
      gene region's score column (column 6) set to 1 for a mono-exonic gene and to any
      other value for a multi-exonic gene.

    -l --level <integer>
      Terminal exon/intron level (i.e. a value of 1 will assess the first terminal
      exon/intron, 2 the second, etc.). Default: $nTerminalLevel
    --min_segdup_overlap <float>
      Minimum fraction of a terminal exon's length that must overlap a segdup region
      for that exon to count as segdup-overlapping. Default: $nMinSegDupOverlap
    --min_gene_overlap <float>
      Minimum fraction of a gene region's length that must be contained within a
      terminal intron for that gene to count as spanned by the intron. Default: $nMinGeneOverlap

    -help
      This help message

HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: terminal level must be a positive integer" unless ($nTerminalLevel =~ /^[1-9]\d*$/);


# Gather isoform start/end coordinates
my %hIsoCoords;
open IN, $sExInGffFile or die "Error: can't open $sExInGffFile: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $sSource, $sType, $nStart, $nEnd, $nScore, $nStr, $nFrame, $sName) = split /\t/;
   next unless ($sType eq 'exon');
   if ($sName =~ /^(\S+)\.ex(\d+):(\d+)$/){
      my $sID   = $1;
      if (exists $hIsoCoords{$sID} ){
         $hIsoCoords{$sID}{start} = $nStart if ($nStart < $hIsoCoords{$sID}{start});
         $hIsoCoords{$sID}{start} = $nEnd   if ($nEnd   < $hIsoCoords{$sID}{start});
         $hIsoCoords{$sID}{end}   = $nStart if ($nStart > $hIsoCoords{$sID}{end});
         $hIsoCoords{$sID}{end}   = $nEnd   if ($nEnd   > $hIsoCoords{$sID}{end});
      }
      else{
         $hIsoCoords{$sID}{chr}   = $sChr;
         $hIsoCoords{$sID}{start} = $nStart;
         $hIsoCoords{$sID}{end}   = $nEnd;
      }
   }
}
close IN;

# Write temporary file with isoform start/end coordinates.
# BUGFIX 1: the GFF input's start coordinate is 1-based; a plain .bed file
# (unlike a file bedtools recognizes as GFF/GTF) is read as 0-based, so the
# start must be converted (-1) here to represent the same genomic base as
# the original GFF record. The end coordinate needs no conversion (a
# 1-based inclusive end is numerically identical to a 0-based exclusive
# end).
# BUGFIX 3: UNLINK=>1 (was 0) so this temp file is actually cleaned up.
my ($fhTmp, $sTmp) = tempfile('isocoord-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1, SUFFIX => '.bed');
foreach my $sID (keys %hIsoCoords){
   print $fhTmp join("\t", $hIsoCoords{$sID}{chr}, $hIsoCoords{$sID}{start} - 1, $hIsoCoords{$sID}{end}, $sID);
   print $fhTmp "\n";
}
close $fhTmp;

# Define output hashes
my %hIsoTermOverlaps;
my %hIsoSegDupContained;
my %hIntLengths;

# Find out which isoforms are fully contained within segdups
open ISOFORMS, "bedtools intersect  -u -f 1.0 -a $sTmp -b $sSegdupBedFile |" or die "Error: could not run bedtools for exon:segdup overlaps: $!\n";
while (<ISOFORMS>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $nStart, $nEnd, $sID) = split /\t/;
   $hIsoSegDupContained{$sID}++;
}
close ISOFORMS;

# BUGFIX 3: the bbox temp file is only needed for the ISOFORMS check above;
# clean it up explicitly now rather than waiting for process exit.
unlink $sTmp;

# Overlap terminal exons with segdups
open EXONS, "bedtools intersect  -u -f $nMinSegDupOverlap -a $sExInGffFile -b $sSegdupBedFile |" or die "Error: could not run bedtools for exon:segdup overlaps: $!\n";
while (<EXONS>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $sSource, $sType, $nStart, $nEnd, $nScore, $nStr, $nFrame, $sName) = split /\t/;
   next unless ($sType eq 'exon');
   if ($sName =~ /^(\S+)\.ex(\d+):(\d+)$/){
      my $sID   = $1;
      my $nExID = $2;
      my $nExCT = $3;

      # Consider first terminal exon level for all transcripts with more exons than the current terminal level
      $hIsoTermOverlaps{$sID}{segdup_first}++ if ( ($nExID <= $nTerminalLevel) and ($nExCT > $nTerminalLevel) );

      # Consider last terminal exon level  for all transcripts with more exons than the current terminal level
      $hIsoTermOverlaps{$sID}{segdup_last}++ if ( ($nExID >= ($nExCT-($nTerminalLevel-1)) ) and ($nExCT > $nTerminalLevel) );
   }
}
close EXONS;


# Overlap terminal introns with genes
open INTRONS, "bedtools intersect  -wo -F $nMinGeneOverlap -a $sExInGffFile -b $sGeneregionBedFile |" or die "Error: could not run bedtools for exon:segdup overlaps: $!\n";
while (<INTRONS>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sExInChr, $sExInSource, $sExInType, $nExInStart, $nExInEnd, $nExInScore, $nExInStr, $nExInFrame, $sExInName,
       $sGeneChr, $sGeneSource, $sGeneType, $nGeneStart, $nGeneEnd, $nGeneScore, $nGeneStr, $nGeneFrame, $sGeneName) = split /\t/;
   next unless ($sExInType eq 'intron');

   if ($sExInName =~ /^(\S+)\.in(\d+):(\d+)$/){
      my $sID    = $1;
      my $nInID  = $2;
      my $nInCT  = $3;
      my $nInLen = $nExInEnd - $nExInStart + 1;

      # Gather first and last intron lengths
      $hIntLengths{$sID}{gene_first} = $nInLen if ( $nInID == $nTerminalLevel);
      $hIntLengths{$sID}{gene_last}  = $nInLen if ( $nInID == ($nInCT-($nTerminalLevel-1)) );

      # Consider first terminal intron level for all transcripts with the same or more introns than the current terminal level
      if ( ($nInID == $nTerminalLevel) and ($nInCT >= $nTerminalLevel) ){
         if ($nGeneScore == 1){
            $hIsoTermOverlaps{$sID}{gene_first}{mono_exon}++
         }
         else{
            $hIsoTermOverlaps{$sID}{gene_first}{multi_exon}++
         }
      }

      # Consider last terminal intron level for all transcripts with the same or more introns than the current terminal level
      if ( ($nInID == ($nInCT-($nTerminalLevel-1)) ) and ($nInCT >= $nTerminalLevel) ){
         if ($nGeneScore == 1){
            $hIsoTermOverlaps{$sID}{gene_last}{mono_exon}++ ;
         }
         else{
            $hIsoTermOverlaps{$sID}{gene_last}{multi_exon}++ ;
         }
      }
   }
}
close INTRONS;

# NOTE: by design, the output loop below only prints isoforms that register
# at least one terminal-exon/intron signal via %hIsoTermOverlaps -- which
# requires having more exons/introns than the chosen --level. A mono-exonic
# (or otherwise too-short-for-the-level) isoform never gets an entry here
# and is never printed, even if it is fully_contained_in_single_segdup; this
# is intentional (see INVESTIGATED, NOT A BUG note above), not an omission.

# Write output
print join("\t", '#ID', "fully_contained_in_single_segdup", "first_terminal_exon_level${nTerminalLevel}_segdup_overlap", "first_terminal_intron_level${nTerminalLevel}_gene_overlap",
           "first_terminal_intron_level${nTerminalLevel}_monoexonic_gene_overlap_count","first_terminal_intron_level${nTerminalLevel}_multiexonic_gene_overlap_count",
           "first_terminal_intron_level${nTerminalLevel}_length",
           "last_terminal_exon_level${nTerminalLevel}_segdup_overlap", "last_terminal_intron_level${nTerminalLevel}_gene_overlap",
           "last_terminal_intron_level${nTerminalLevel}_monoexonic_gene_overlap_count","last_terminal_intron_level${nTerminalLevel}_multiexonic_gene_overlap_count",
           "last_terminal_intron_level${nTerminalLevel}_length"), "\n";
foreach my $sID (sort keys %hIsoTermOverlaps){
   my $sSegdupContained       = exists $hIsoSegDupContained{$sID} ? 'yes' : 'no';
   my $sSegdupOverlapFirst = 'no';
   if ( $hIsoTermOverlaps{$sID}{segdup_first} ){
      $sSegdupOverlapFirst = 'yes' if ( $hIsoTermOverlaps{$sID}{segdup_first} == $nTerminalLevel)
   }
   my $sSegdupOverlapLast = 'no';
   if ( $hIsoTermOverlaps{$sID}{segdup_last} ){
      $sSegdupOverlapLast = 'yes' if ( $hIsoTermOverlaps{$sID}{segdup_last} == $nTerminalLevel)
   }

   my $sGeneOverlapFirst      = exists $hIsoTermOverlaps{$sID}{gene_first}  ? 'yes' : 'no';
   my $sGeneOverlapLast       = exists $hIsoTermOverlaps{$sID}{gene_last}   ? 'yes' : 'no';
   my $nGeneOverlapFirstMono  = exists $hIsoTermOverlaps{$sID}{gene_first}{mono_exon}  ? $hIsoTermOverlaps{$sID}{gene_first}{mono_exon}  : 0;
   my $nGeneOverlapFirstMulti = exists $hIsoTermOverlaps{$sID}{gene_first}{multi_exon} ? $hIsoTermOverlaps{$sID}{gene_first}{multi_exon} : 0;
   my $nGeneOverlapLastMono   = exists $hIsoTermOverlaps{$sID}{gene_last}{mono_exon}   ? $hIsoTermOverlaps{$sID}{gene_last}{mono_exon}   : 0;
   my $nGeneOverlapLastMulti  = exists $hIsoTermOverlaps{$sID}{gene_last}{multi_exon}  ? $hIsoTermOverlaps{$sID}{gene_last}{multi_exon}  : 0;
   my $nFirstIntLen           = exists $hIntLengths{$sID}{gene_first} ? $hIntLengths{$sID}{gene_first} : 0;
   my $nLastIntLen            = exists $hIntLengths{$sID}{gene_last}  ? $hIntLengths{$sID}{gene_last}  : 0;
   print join("\t", $sID, $sSegdupContained,
              $sSegdupOverlapFirst, $sGeneOverlapFirst, $nGeneOverlapFirstMono, $nGeneOverlapFirstMulti, $nFirstIntLen,
              $sSegdupOverlapLast, $sGeneOverlapLast, $nGeneOverlapLastMono, $nGeneOverlapLastMulti, $nLastIntLen), "\n";
}


# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup.
# BUGFIX 3: this was never wired to a signal in the original (see $SIG
# assignment below, which is new), and even when called it didn't actually
# unlink anything -- it now does both.
$SIG{INT}  = \&INTERRUPT;
$SIG{TERM} = \&INTERRUPT;
sub INTERRUPT{
   unlink $sTmp if (defined $sTmp and -e $sTmp);
   exit(1);
}
