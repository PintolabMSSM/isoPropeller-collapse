#!/bin/sh

# 29.07.2025 10:13:11 EDT


cd /sc/arion/scratch/pintod02/isopropeller-collapse-test/04_isoPropeller-merge

# Arguments that need to be defined in the main config.yaml
REFGEN="/sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa"
REFGTF="/sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/gencode.v41.annotation_no-mi-sn-sno-sca-RNAs.gtf"
REFTSS="/sc/arion/work/pintod02/opt/isoseq_pipeline/data/combined_cage_Pitt-Fantom5-119FrontalLob_refTSS3.3_refseq_gencode_extended.merged.sorted_chr.bed"
RMSKBED="/sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome_rmsk_nr.bed"
PARBED="/hpc/users/pintod02/opt/isoseq_pipeline/data/GRCh38.p13_PARs.bed"
SEGDUP="/sc/arion/work/pintod02/opt/isoseq_pipeline/data/SegDups_hg38.merged.bed"

FILT_MONOEXON_TSS_MAXDIST_BP=10
FILT_MONOEXON_MIN_INTRON_OVLP_BP=10
FILT_RMSK_MIN_OVLP_FRACT=0.9
FILT_PAR_MIN_OVLP_FRACT=0
FILT_TPM_MIN_COUNT=10
FILT_TPM_MIN_FRACT=0.5

# Arguments that are passed from previous tasks
SAMPLE="ISOP_all"

ISOGTF="ISOP_all.gtf"
ISOBED="ISOP_all.bed"
ISOTSS="ISOP_all_tss.bed"
ISOEXP="ISOP_all_exp.txt"

REFBED="/sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/gencode.v41.annotation.bed"
REFFAI="/sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa.fai"


# Monoexon tss filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-tss-overlap.py \
   --isoform_bed12       ${ISOBED} \
   --isoform_tss_bed     ${ISOTSS} \
   --reftss_bed          ${REFTSS} \
   --genome_index        ${REFFAI} \
   --max_distance        ${FILT_MONOEXON_TSS_MAXDIST_BP} \
   --out_bed             isoqc_fail_${SAMPLE}_monoexon-no-reftss-overlap.bed \
   --out_ids             isoqc_fail_${SAMPLE}_monoexon-no-reftss-overlap.ids 

# Monoexon pre-mRNA filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-premrna-fragments.py \
   --isoform_bed12       ${ISOBED} \
   --reference_bed12     ${REFBED} \
   --min_intron_overlap  ${FILT_MONOEXON_MIN_INTRON_OVLP_BP} \
   --out_bed             isoqc_fail_${SAMPLE}_monoexon-likely-premrnas.bed \
   --out_ids             isoqc_fail_${SAMPLE}_monoexon-likely-premrnas.ids

# Noncanonical splice junction filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_multiexon-noncanonical-splices.py \
   --isoform_bed12       ${ISOBED} \
   --genome_fasta        ${REFGEN} \
   --out_bed             isoqc_fail_${SAMPLE}_multiexonic-noncanonical-splices.bed \
   --out_ids             isoqc_fail_${SAMPLE}_multiexonic-noncanonical-splices.ids \
   --out_motifs          isoqc_fail_${SAMPLE}_multiexonic-noncanonical-splices.motifs.txt

# Template switching filter (based on SQANTI5 approach)
~/opt/isoPropeller-collapse/workflow/scripts/filter_multiexon-rt-switching.py \
   --isoform_bed12       ${ISOBED} \
   --genome_fasta        ${REFGEN} \
   --out_rts_tsv         isoqc_fail_${SAMPLE}_multiexonic-rt-switching.repeats.txt \
   --out_ids             isoqc_fail_${SAMPLE}_multiexonic-rt-switching.ids \
   --out_bed             isoqc_fail_${SAMPLE}_multiexonic-rt-switching.bed

# Antisense perfect splice match filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_multiexon-antisense-splicechain-match.py \
   --isoform_bed12       ${ISOBED} \
   --reference_bed12     ${REFBED} \
   --out_ids             isoqc_fail_${SAMPLE}_multiexonic-antisense-splicechain-match.ids \
   --out_bed             isoqc_fail_${SAMPLE}_multiexonic-antisense-splicechain-match.bed

# Reference region overlap filter, contained in repeats
~/opt/isoPropeller-collapse/workflow/scripts/filter_reference-overlap-on-exons.py \
   --isoform_bed12        ${ISOBED} \
   --reference_bed12      ${RMSKBED} \
   --min_overlap_fraction ${FILT_RMSK_MIN_OVLP_FRACT} \
   --out_ids              isoqc_fail_${SAMPLE}_repeatmasker-overlap.ids \
   --out_bed              isoqc_fail_${SAMPLE}_repeatmasker-overlap.bed \
   --out_stats            isoqc_fail_${SAMPLE}_repeatmasker-overlap.stats.txt

# Reference region overlap filter, contained in repeats
~/opt/isoPropeller-collapse/workflow/scripts/filter_reference-overlap-on-exons.py \
   --isoform_bed12        ${ISOBED} \
   --reference_bed12      ${PARBED} \
   --min_overlap_fraction ${FILT_PAR_MIN_OVLP_FRACT} \
   --out_ids              isoqc_fail_${SAMPLE}_PAR-overlap.ids \
   --out_bed              isoqc_fail_${SAMPLE}_PAR-overlap.bed \
   --out_stats            isoqc_fail_${SAMPLE}_PAR-overlap.stats.txt 

# Min-TPM filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_TPM-fraction.py \
   --count_matrix         ${ISOEXP} \
   --isoform_bed12        ${ISOBED} \
   --min_tpm              ${FILT_TPM_MIN_COUNT} \
   --min_fraction_samples ${FILT_TPM_MIN_FRACT} \
   --out_ids              isoqc_fail_${SAMPLE}_min-TPM.ids \
   --out_bed              isoqc_fail_${SAMPLE}_min-TPM.bed


# Mismatched terminal exon in segdups filter
bed2intronexongff.pl -v 1 ${ISOBED} > isoqc_temp_${SAMPLE}_corrected.intronexon.gff
gtf-get-gene-regions.pl   ${REFGTF} > isoqc_temp_${SAMPLE}_reference-gene-regions.gtf
for LEVEL in 1 2 3 4
do
   ~/opt/isoPropeller-collapse/workflow/scripts/filter_segdup-mismapped-terminal-exons.pl \
      -i isoqc_temp_${SAMPLE}_corrected.intronexon.gff \
      -g isoqc_temp_${SAMPLE}_reference-gene-regions.gtf \
      -s ${SEGDUP} \
      -l ${LEVEL} \
      > isoqc_temp_${SAMPLE}_terminal-exons-in-segdup_${LEVEL}.txt
   
   awkt '($2=="no" && $3=="yes" && $6>0 && $7>100000) || ($2=="no" && $8=="yes" && $11>0 && $12>100000) {print $1}' \
      isoqc_temp_${SAMPLE}_terminal-exons-in-segdup_${LEVEL}.txt \
      > isoqc_temp_${SAMPLE}_mismapped-terminal-exon-in-segdup_${LEVEL}.txt
done
cat isoqc_temp_${SAMPLE}_mismapped-terminal-exon-in-segdup_*.txt | sort | uniq > isoqc_fail_${SAMPLE}_mismapped-terminal-exon-in-segdup.ids
intersect-by-ids -ff ${ISOBED} -fc 4 -if isoqc_fail_${SAMPLE}_mismapped-terminal-exon-in-segdup.txt > isoqc_fail_${SAMPLE}_mismapped-terminal-exon-in-segdup.bed
rm -f isoqc_temp_${SAMPLE}_*

