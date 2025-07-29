#!/bin/sh

# 29.07.2025 10:13:11 EDT


cd /sc/arion/scratch/pintod02/isopropeller-collapse-test/04_isoPropeller-merge

# Monoexon tss filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-tss-overlap.py \
   --isoform_bed12       ISOP_all.bed \
   --isoform_tss_bed     ISOP_all_tss.bed \
   --reftss_bed          /sc/arion/work/pintod02/opt/isoseq_pipeline/data/combined_cage_Pitt-Fantom5-119FrontalLob_refTSS3.3_refseq_gencode_extended.merged.sorted_chr.bed \
   --out_bed             filter-test_monoexon-no-reftss-overlap.bed \
   --out_ids             filter-test_monoexon-no-reftss-overlap.ids \
   --max_distance        10 \
   --genome_index        /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa.fai

# Monoexon pre-mRNA filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-premrna-fragments.py \
   --isoform_bed12       ISOP_all.bed \
   --reference_bed12     /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/gencode.v41.annotation.bed \
   --out_bed             filter-test_monoexon-likely-premrnas.bed \
   --out_ids             filter-test_monoexon-likely-premrnas.ids \
   --min_intron_overlap  10

# Noncanonical splice junction filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_multiexon-noncanonical-splices.py \
   --isoform_bed12       ISOP_all.bed \
   --genome_fasta        /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa \
   --out_bed             filter-test_multiexonic-noncanonical-splices.bed \
   --out_ids             filter-test_multiexonic-noncanonical-splices.ids \
   --out_motifs          filter-test_multiexonic-noncanonical-splices.motifs.txt

# Template switching filter (based on SQANTI5 approach)
~/opt/isoPropeller-collapse/workflow/scripts/filter_multiexon-rt-switching.py \
   --isoform_bed12       ISOP_all.bed \
   --genome_fasta        /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa \
   --out_rts_tsv         filter-test_multiexonic-rt-switching.repeats.txt \
   --out_ids             filter-test_multiexonic-rt-switching.ids \
   --out_bed             filter-test_multiexonic-rt-switching.bedl

# Antisense perfect splice match filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_multiexon-antisense-splicechain-match.py \
   --isoform_bed12       ISOP_all.bed \
   --reference_bed12     /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/gencode.v41.annotation.bed \
   --out_ids             filter-test_multiexonic-antisense-splicechain-match.ids \
   --out_bed             filter-test_multiexonic-antisense-splicechain-match.bed
   
# Reference region overlap filter, contained in repeats
~/opt/isoPropeller-collapse/workflow/scripts/filter_reference-overlap-on-exons.py \
   --isoform_bed12        ISOP_all.bed \
   --reference_bed12      /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome_rmsk_nr.bed \
   --out_ids              filter-test_repeatmasker-overlap.ids \
   --out_bed              filter-test_repeatmasker-overlap.bed \
   --out_stats            filter-test_repeatmasker-overlap.stats.txt \
   --min_overlap_fraction 0.9

# Reference region overlap filter, contained in repeats
~/opt/isoPropeller-collapse/workflow/scripts/filter_reference-overlap-on-exons.py \
   --isoform_bed12        ISOP_all.bed \
   --reference_bed12      /hpc/users/pintod02/opt/isoseq_pipeline/data/GRCh38.p13_PARs.bed \
   --out_ids              filter-test_PAR-overlap.ids \
   --out_bed              filter-test_PAR-overlap.bed \
   --out_stats            filter-test_PAR-overlap.stats.txt \
   --min_overlap_fraction 0

