## isoPropeller-collapse
**isoPropeller-collapse** is a Snakemake pipeline designed for the non-redundant reconstruction of transcriptomes from long-read sequencing data (e.g., PacBio HiFi/Iso-Seq). It provides the discovery engine of the **isoPropeller** suite, transforming long read sequencing data into a high-confidence, filtered, and consolidated isoform reference.

The pipeline automates the long read sequencing workflow to collapse reads into isoforms, including:

- **High-Accuracy Mapping:** Utilizing [minimap2](https://github.com/lh3/minimap2) for splice-aware alignment.
- **Reference-Based Correction:** Parallelized error correction via [TranscriptClean](https://github.com/mortazavilab/TranscriptClean) to refine splice junctions and indels.
- **De Novo Isoform Clustering:** Grouping reads into unique transcript models based on junction chains using our [isoPropeller](https://github.com/PintolabMSSM/isoPropeller) tool.
- **Modular Quality Control:** A multi-stage filtering suite that removes technical artifacts like template switching, internal priming, and non-canonical splicing.
- **Expression Redistribution:** Consolidation of fragmented reads into full-length parent transcripts with proportional count preservation.
- **Fusion Discovery:** Identification of chimeric transcripts and genomic breakpoints using the PacBio [pbfusion](https://github.com/PacificBiosciences/pbfusion) tool.

By the end of the collapse process, the pipeline produces a refined, non-redundant GTF and expression matrix, formatted for downstream functional characterization in the [isoPropeller-annotate](https://github.com/PintolabMSSM/isoPropeller-annotate) pipeline of the isoPropeller suite.



## Installation


#### Prepare the snakemake conda environment

Installation of the required external software packages is largely handled by the pipeline itself, however a conda environment named `snakemake` needs to be present in your environment. We recommend miniconda, which is a free minimal installer for [conda](https://docs.conda.io/en/latest/miniconda.html). Follow the instructions below to start the miniconda installer on Linux. When asked whether the conda environment should automatically be initialized, select 'yes'. Note that Snakemake requires the channel_priority to be set to strict. The post-installation commands to apply this setting are included in the post-installation selection below.

```bash
# Start miniconda installation
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh

# Post-installation commands to enforce strict channel_priority (required for Snakemake)
conda config --set auto_activate_base false
conda config --set channel_priority strict
```

After installation and configuration of conda/miniconda, the following 'conda create' command can be used to set up the required 'snakemake' environment.

```bash
conda create -c conda-forge -c bioconda -n snakemake 'snakemake>=9.8' snakemake-executor-plugin-lsf snakemake-executor-plugin-cluster-generic 'tabulate>=0.8'
```



## Running isoPropeller collapse

Once the snakemake environment has been created, all that is needed to run the pipeline is an analysis folder containing a single input file listing the the paths to long read sequencing data in bam or fastq format and the sample names. For fastq input, the csv file should have two columns and include a header with the names `path_to_fastq` and `sample`, with the paths to the fastq(.gz) and sample name, respectively. Multiple fastq files can be specified per sample. These will be merged together automatically in the isoPropeller-collapse analysis.

```bash
path_to_fastq,sample
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep1.fastq.gz,WTC11_CapTrap_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep2.fastq.gz,WTC11_CapTrap_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep3.fastq.gz,WTC11_CapTrap_rep3
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep1.fastq.gz,WTC11_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep2.fastq.gz,WTC11_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep3.fastq.gz,WTC11_rep3
```

The structure is similar for bam inputs, except that the header will have the name `path_to_bam` and `sample`, e.g.:

```
path_to_bam,sample
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep1.bam,WTC11_CapTrap_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep2.bam,WTC11_CapTrap_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep3.bam,WTC11_CapTrap_rep3
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep1.bam,WTC11_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep2.bam,WTC11_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep3.bam,WTC11_rep3
```

Once the input file is ready, run start the pipeline as follows with the provided wrapper script:

```
submitjob 72 -c 4 -m 20 -q premium -P acc_pintod02b \
   ~/opt/isoPropeller-collapse/run-isoPropeller-collapse \
   -i sample-input.csv
```



## Overview of pipeline outputs

The isoPropeller-collapse pipeline is organized as a series of tasks, each of which produces their own output folder. An overview of each task and the outputs it produces is provided below.



### 01_mapping

**Description:** This directory contains the initial processing of full-length non-chimeric (FLNC) long-read sequencing data. It handles the parallel concatenation of input files, performs high-accuracy splice-aware alignment to the reference genome, and prepares the data for the collapse stages.

**Contents:**

- **`flnc_merged.fastq.gz`**: The consolidated Full-Length Non-Chimeric (FLNC) reads for each sample. This ensures that even if the data was delivered in multiple parts, it is treated as a single cohesive unit for alignment.
- **`mapped.bam`**: The raw alignment file generated by **minimap2** using the `splice:hq` preset. This identifies the genomic coordinates of every read, accounting for large introns and splice junctions.
- **`{sample}_mapped_labeled.bam`**: A "TALON-ready" BAM file. This version of the alignment has been processed by `talon_label_reads` to identify and tag potential internal priming artifacts or specific read-level attributes.
- **`{sample}_mapped_fa_read_labels.tsv.gz`**: A compressed report detailing the labeling results for every read in the sample. This is used in subsequent steps to filter out low-confidence alignments or technical artifacts.



### 02_transcriptclean

**Description:** This directory contains the outputs of the transcript cleaning process performed using [TranscriptClean](https://github.com/mortazavilab/TranscriptClean). The pipeline splits large alignment files into manageable chunks, corrects non-canonical splice junctions and micro-indels against the reference genome, and performs a controlled downsampling of mitochondrial reads for samples containing a high fraction of mitochondrial sequences.

**Contents:**

- **`split/`**: A temporary subdirectory (managed via Snakemake checkpoints) containing the original BAM file partitioned into SAM chunks. Each chunk is housed in its own numbered subdirectory to ensure isolated processing.
- **`split_logs/`**: Compressed logs for each individual chunk, tracking the specific corrections made by `transcriptclean` during the scatter phase.
- **`{sample}_mapped_labeled_tclean_temp.bam`**: The intermediate merged result. This file represents the consolidated output of all corrected chunks before final processing.
- **`{sample}_mapped_labeled_tclean.bam`**: The final, high-quality alignment file. This version has been corrected for structural errors and includes a downsampled **chrM** (mitochondrial) fraction to prevent alignment depth bottlenecks.
- **`{sample}_final_clean.log.gz`**: A comprehensive, consolidated log file summarizing all corrections (mismatches, indels, and splice sites) made across the entire sample.
- **`{sample}_final_clean.TE.log.gz`**: A dedicated log focusing on **Transcript Error** statistics, providing a detailed breakdown of the error profile for the original library.



### 03_isoPropeller

**Description:** This directory contains the primary unfiltered transcript models for each sample. The process uses **isoPropeller** to cluster reads based on their splice junction chains and genomic coordinates. For each splice chain it also retains transcription start site (TSS) information to better define the 5' ends of transcripts before filtering out low-support "singleton" isoforms.

**Contents:**

- **`{sample}_all.gtf`**: The comprehensive, unfiltered set of all transcript models discovered in the sample. This includes potential rare isoforms and technical noise.
- **`{sample}_all_end_dist.txt`**: A statistical report detailing the distance of transcript ends relative to reference TSS/TTS sites, useful for assessing the completeness of the full-length reads.
- **`{sample}_all_stat.txt`**: A summary file providing high-level metrics on the number of genes and isoforms discovered per sample.
- **`{sample}_depth-gt1.gtf`**: The high-confidence per-sample transcriptome. This file is filtered to keep only those transcripts supported by **two or more reads** (depth > 1), significantly reducing the rate of false-positive isoforms.
- **`{sample}_all.log`**: The internal execution log from the isoPropeller engine, detailing the clustering parameters and any warnings encountered during processing.



### 04_isoPropeller-merge

**Description:** This directory contains the consolidated results for the entire sample set. The pipeline merges individual sample GTF files, aggregates expression counts across samples, and clusters transcription start sites (TSS) and transcription termination sites (TTS) to refine the structural boundaries of the final transcript models.

**Contents:**

- **`{prefix}_{suffix}.gtf`**: The initial merged transcriptome containing all unique isoforms identified across the cohort.
- **`{prefix}_{suffix}_exp.txt`**: A comprehensive expression matrix providing read counts (depth) for every transcript across every sample in the study.
- **`{prefix}_{suffix}_id.txt`**: A simple list of the unified transcript IDs used for downstream filtering and end-region analysis.
- **`{prefix}_{suffix}_tss.bed` / `_tts.bed`**: Genomic coordinates of clustered start and end regions, formatted for easy visualization in genome browsers.
- **`{prefix}_{suffix}_tss_count.txt` / `_tts_count.txt`**: Quantitative data for each identified TSS and TTS, used to determine which ends are the most biologically prevalent (modal).
- **`{prefix}_{suffix}_modal_ends.gtf`**: The refined cohort transcriptome. This file is an updated version of the merged GTF where transcript boundaries have been adjusted to the "modal" (most common) TSS and TTS positions found across the samples.



### 05_isoPropeller-filter

**Description:** This directory contains the results of the isoform filtering phase. The pipeline applies a series of modular filters to identify and remove transcripts based on structural flaws, potential artifacts (like template switching), and low expression support. Transcripts that fail any of the active filters are logged as "fail" IDs, while the remaining transcripts are aggregated into a final, high-quality "pass" dataset. The filters can be configured and turned on or off based on parameters defined in the `config.yaml` configuration file.

**Contents:**

- **`filt_\*/`**: Individual subdirectories for each active filter (e.g., `filt_monoexon_tss`, `filt_noncanonical_splice`). These contain:
  - **`isoqc_fail_\*.ids`**: A list of transcript IDs that failed that specific filter.
  - **`isoqc_fail_\*.bed`**: A BED file for visualizing the failed transcripts in a genome browser.
- **`{prefix}_{suffix}_isoqc_pass.gtf`**: The final, filtered transcriptome containing only high-confidence isoforms with boundaries standardized to the max-transcript ends.
- **`{prefix}_{suffix}_isoqc_pass_modal_ends.gtf`**: The filtered transcriptome with structural boundaries standardized to the most frequent (modal) start and end sites.
- **`{prefix}_{suffix}_isoqc_pass_exp.txt`**: A cleaned expression matrix containing only the "pass" isoforms, with headers standardized for the **isoPropeller-annotate** pipeline.
- **`{prefix}_{suffix}_isoqc_fail.ids`**: A master list of every transcript ID that was rejected by at least one filter.
- **`{prefix}_{suffix}_isoqc_pass.trackgroups`**: A default group assignment file (mapping all samples to a group labeled "ALL") to facilitate downstream visualization.



### 06_qc-reports

**Description:** This directory serves as the centralized repository for all Quality Control (QC) outputs. It spans the entire workflow from raw FASTQ statistics and alignment accuracy to cohort-level filtering summaries and proteogenomic validation prep. These reports provide information for identifying sample outliers, sequencing artifacts, and the overall efficiency of the isoform discovery process.

**Contents:**

#### flnc-fastqc / flnc-seqkit-stats

- **`flnc_merged_fastqc.html/zip`**: Standard quality metrics (per-base quality, GC content, overrepresented sequences) for the merged long reads.
- **`seqkit_flnc_wide.stats.tsv`**: A cohort-wide summary of read lengths, N50 values, and total base counts, facilitating quick comparisons of library preparation efficiency across samples.

#### mapped-rnaseqc / mapped-picard-RnaSeqMetrics

- **`rna_seqc_summary_wide.tsv`**: High-level alignment metrics, including gene body coverage, 3' bias, and the distribution of reads across genomic features (exonic, intronic, intergenic).
- **`{sample}.RnaSeqMetrics.pdf`**: Picard's visualization of read distribution and coverage, providing a classic "ribbon" plot of transcript coverage.
- **`collapsed_reference.gtf`**: A specialized, non-overlapping version of the reference annotation used to calculate unambiguous gene-level metrics.

#### mapped-bamqc

- **`flagstat.txt` / `idxstats.txt`**: Samtools summaries of alignment totals and per-chromosome read distributions.
- **`chrM_status.txt`**: A specific QC flag indicating if a sample passed the mitochondrial read depth threshold (`MAXCHRMREADS`).
- **`bam_qc_summary_wide.tsv`**: An aggregated cohort table summarizing mapping percentages and duplication rates.

#### flnc-longreadsum

- **`flnc-longreadsum-fastq-report/`**: Specialized long-read QC focusing on read length distributions and quality scores tailored for PacBio/Nanopore data.
- **`mapped-longreadsum-rnaseq-report/`**: Context-aware QC that evaluates how well long reads span known gene structures and identifies potential truncation or internal priming.

#### isoform-filtering

- **`isoform-filter-stats.totals.tsv`**: A "balance sheet" for the discovery process. It tracks how many unique isoforms were originally identified and how many were removed by each specific filter (e.g., non-canonical splice, template switching).
- **`isoform-filter-stats.per-sample-wide.tsv`**: A crucial table showing the conversion rate of raw reads into finalized "pass" isoforms for every sample, including the percentage of reads successfully assigned to the final transcriptome.

#### multiqc

- **`multiqc_report.html`**: The definitive project summary. This interactive report aggregates results from FastQC, RNA-SeQC, Picard, and LongReadSum into a single searchable dashboard.



### 07_isoPropeller-defrag

**Description:** This directory contains the "defragmented" or consolidated transcriptome. The pipeline identifies "contained" transcripts—those whose splice chains are exact subsets of longer "parent transcripts" and merges them. An additional count file is provided where expression counts from these fragments are not simply discarded but are proportionally redistributed to their parent transcripts.

**Contents:**

- **`{prefix}_{suffix}_isoqc_pass_containment_map.tsv`**: A detailed mapping file showing which truncated transcript IDs were collapsed into which full-length parent transcripts.
- **`{prefix}_{suffix}_isoqc_pass_defrag.gtf`**: The final, non-redundant transcriptome. All redundant fragments have been removed, leaving only the most complete representative for each unique splice chain structure.
- **`{prefix}_{suffix}_isoqc_pass_defrag_exp.txt`**: The updated expression matrix. In this file, the read counts originally assigned to fragments have been added to the counts of their respective parent transcripts.
- **`{prefix}_{suffix}_isoqc_pass_defrag_exp_redist.txt`**: A log of the redistribution process, detailing how many counts were moved and the logic used for the proportional allocation.
- **`{prefix}_{suffix}_isoqc_pass_defrag_modal_ends.gtf`**: The defragmented transcriptome utilizing standardized "modal" 5' and 3' boundaries for the final high-confidence models.



This stage of the **isoPropeller-collapse** pipeline focuses on noise reduction at the locus level. It implements a relative expression filter to remove rare isoforms that contribute minimally to the total expression of a gene or cluster, further refining the transcriptome to its most biologically relevant components.



### 08_isoPropeller-defrag-pruned

**Description:** This directory contains an additional "pruned" version of the collapsed transcriptome. After the defragmentation step (07) has consolidated fragments, this rule applies a relative abundance filter within each transcript cluster (locus). It identifies and removes "rare" isoforms that fall below a specific percentage of the total locus expression, ensuring that the final dataset isn't dominated by low-frequency splice variants that may represent biological noise or minor sequencing artifacts.

**Contents:**

- **`{prefix}_{suffix}_isoqc_pass_defrag_pruned.gtf`**: The high-confidence, non-redundant transcriptome. This is the primary structural output of the collapse pipeline, having passed through mapping, cleaning, clustering, structural filtering, defragmentation, and finally, abundance-based pruning.
- **`{prefix}_{suffix}_isoqc_pass_defrag_pruned_exp.txt`**: The final expression matrix for the project. Only transcripts that passed the pruning threshold are retained here.
- **`{prefix}_{suffix}_isoqc_pass_defrag_pruned_clusters.txt`**: A reference file mapping transcripts to their respective clusters (loci), used by the pruning script to calculate total locus expression.
- **`{prefix}_{suffix}_isoqc_pass_defrag_pruned_dropped.txt`**: A log of the specific transcripts that were pruned during this step, along with their relative expression metrics for audit purposes.
- **`{prefix}_{suffix}_isoqc_pass_defrag_pruned_tss.bed` / `_tts.bed`**: The coordinate-accurate start and end sites for the remaining pruned transcripts.



This stage of the **isoPropeller-collapse** pipeline is dedicated to the identification of gene fusion events. Using PacBio’s specialized **pbfusion** toolset, the pipeline discovers chimeric transcripts that represent the joining of two or more distinct genomic loci, a critical feature in cancer research and complex genomic rearrangements.



### 09_pbfusion

**Description:** This directory contains the discovery and characterization of fusion transcripts using the PacBio [pbfusion](https://github.com/PacificBiosciences/pbfusion) tool. The process begins with a specialized alignment using `pbmm2` (configured for Iso-Seq data), followed by the generation of a binary reference cache to accelerate fusion detection. The final step identifies high-confidence breakpoints and the specific transcripts supporting each fusion event.

**Contents:**

- **`{sample}_pbfusion.bam`**: A specialized alignment file generated by `pbmm2` using the `ISOSEQ` preset. This alignment is optimized to allow for split-read mapping, which is essential for identifying transcripts that originate from different chromosomes or distant genomic regions.
- **`reference_cache/gtf.bin`**: A binary representation of the reference GTF. Pre-caching the reference significantly reduces the memory footprint and runtime of the fusion discovery engine.
- **`{sample}.breakpoints.bed`**: A BED file containing the exact genomic coordinates of the discovered fusion breakpoints.
- **`{sample}.breakpoints.groups.bed`**: Groups together multiple breakpoints that appear to belong to the same complex rearrangement event.
- **`{sample}.transcripts`**: A detailed text report listing every transcript that provides evidence for a fusion, including quality scores and read support.
- **`{sample}.unannotated.bed` / `.clusters.bed`**: Reports on chimeric reads that do not overlap with known gene annotations, potentially representing novel fusion partners or intergenic rearrangements.
