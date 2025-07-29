#!/bin/sh

# 29.07.2025 10:13:11 EDT


cd /sc/arion/scratch/pintod02/isopropeller-collapse-test/04_isoPropeller-merge

# Monoexon tss filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-tss-overlap.py \
   --isoform_bed12       ISOP_all.bed \
   --isoform_tss_bed     ISOP_all_tss.bed \
   --reftss_bed          /sc/arion/work/pintod02/opt/isoseq_pipeline/data/combined_cage_Pitt-Fantom5-119FrontalLob_refTSS3.3_refseq_gencode_extended.merged.sorted_chr.bed \
   --out_bed             filter-test-monoexon-tss.bed \
   --out_ids             filter-test-monoexon-tss.ids \
   --max_distance        10 \
   --genome_index        /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa.fai

# Monoexon pre-mRNA filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-premrna-fragments.py \
   --isoform_bed12       ISOP_all.bed \
   --reference_bed12     /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/gencode.v41.annotation.bed \
   --out_bed             filter-test-monoexon-premrna.bed \
   --out_ids             filter-test-monoexon-premrna.ids \
   --min_intron_overlap  10

# Noncanonical splice junction filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_multiexon-noncanonical-splices.py \
   --isoform_bed12       ISOP_all.bed \
   --genome_fasta        /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa \
   --out_bed             filter-test-multiexonic-noncanonical-splices.bed \
   --out_ids             filter-test-multiexonic-noncanonical-splices.ids \
   --out_motifs          filter-test-multiexonic-noncanonical-splices.motifs.txt


# Antisense perfect splice match filter
remove_isoforms_with_antisense_perfect_splicematch: True


# Reference region overlap filter
# We allow specifying any overlap or a fraction of overlap relative to the isoform
remove_isoforms_overlapping_PAR_regions:            True  # Here we do any overlap
remove_isoforms_fully_contained_in_repeats:         True  # Here we specify the degree of overlap (e.g. >90%)

# Template switching filter
remove_isoforms_with_template_switching_artifacts:  True

