#!/bin/sh

# 29.07.2025 10:13:11 EDT


cd /sc/arion/scratch/pintod02/isopropeller-collapse-test/04_isoPropeller-merge

# Monoexon tss filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-tss-overlap.py \
   --isoform_gtf   ISOP_all.gtf \
   --reftss_bed    /sc/arion/work/pintod02/opt/isoseq_pipeline/data/combined_cage_Pitt-Fantom5-119FrontalLob_refTSS3.3_refseq_gencode_extended.merged.sorted_chr.bed \
   --out_gtf       filter-test-monoexon-tss.gtf \
   --out_ids       filter-test-monoexon-tss.ids \
   --max_distance  50 \
   --genome_index  /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa.fai

# Monoexon pre-mRNA filter
~/opt/isoPropeller-collapse/workflow/scripts/filter_monoexon-premrna-fragments.py \
   --isoform_gtf         ISOP_all.gtf \
   --reference           /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/gencode.v41.annotation.bed \
   --genome_index        /sc/arion/projects/pintod02c/reference-databases/hg38-v41-ERCC/GRCh38.primary_assembly.genome.fa.fai \
   --out_gtf             filter-test-monoexon-premrna.gtf \
   --out_ids             filter-test-monoexon-premrna.ids \
   --min_intron_overlap  10
