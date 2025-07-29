#!/usr/bin/perl

# 02.01.2022 13:16:28 EST

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);

# GLOBALS
$ENV{TMPDIR}   ||= "/sc/arion/scratch/$ENV{USER}/tmp/";   # location for tmp file storage

# GET PARAMETERS
my $sHelp                = 0;
my $sExInGffFile         = '';
my $sSegdupBedFile       = '';
my $sGeneregionBedFile   = '';
my $nTerminalLevel       = 1;
my $nMinSegDupOverlap    = 0.90;
my $nMinGeneOverlap      = 1.00;
GetOptions("help!"    => \$sHelp,
           "input:s"  => \$sExInGffFile,
           "segdup:s" => \$sSegdupBedFile,
           "genes:s"  => \$sGeneregionBedFile,
           "level:i"  => \$nTerminalLevel);

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
      File with gene regions in 
     
    -l --level <integer>
      Terminal exon/intron level (i.e. a value of 1 will assess the first terminal 
      exon/intron, 2 the second, etc.). Default: $nTerminalLevel
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: terminal level must be an integer greater than 1" unless ($nTerminalLevel =~ /^[1-9]\d*$/);


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

# Write temporary file with isoform start/end coordinates
my ($fhTmp, $sTmp) = tempfile('isocoord-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>0, SUFFIX => '.bed');
foreach my $sID (keys %hIsoCoords){
   print $fhTmp join("\t", $hIsoCoords{$sID}{chr}, $hIsoCoords{$sID}{start}, $hIsoCoords{$sID}{end}, $sID);
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
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
