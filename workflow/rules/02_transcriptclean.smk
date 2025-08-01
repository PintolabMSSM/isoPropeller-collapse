# ───────────────────────────────────────────────
# Rule: Convert BAM to SAM for transcriptclean
# ───────────────────────────────────────────────
rule bam_to_tc_sam:
    message: "Convert BAM to SAM for transcriptclean: {wildcards.sample}"
    input:
        bam = "01_mapping/{sample}/{sample}_mapped_labeled.bam"
    output:
        sam = temp("02_transcriptclean/{sample}/{sample}.sam")
    log:
        "logs/02_transcriptclean/{sample}_bam_to_sam.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_bam_to_sam.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        "mkdir -p 02_transcriptclean/{wildcards.sample} && "
        "samtools view -h {input.bam} > {output.sam} 2>> {log}"

# ───────────────────────────────────────────────
# Rule: Run transcriptclean
# ───────────────────────────────────────────────
rule run_transcriptclean:
    message: "Running transcriptclean on {wildcards.sample}"
    input:
        sam = "02_transcriptclean/{sample}/{sample}.sam",
        ref = GENOMEFASTA,
        fai = GENOMEFASTA + ".fai"
    output:
        cleaned_sam   = temp("02_transcriptclean/{sample}/{sample}_clean.sam"),
        clean_log     = "02_transcriptclean/{sample}/{sample}_clean.log.gz",
        clean_te_log  = "02_transcriptclean/{sample}/{sample}_clean.TE.log.gz"
    log:
        "logs/02_transcriptclean/{sample}_transcriptclean.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_transcriptclean.txt"
    params:
        junctions_file = SPLICEJUNCTIONS,
        outprefix      = "02_transcriptclean/{sample}/{sample}"
    threads: 24
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        """
        transcriptclean \
            -s {input.sam} \
            -g {input.ref} \
            -j {params.junctions_file} \
            -t {threads} \
            --primaryOnly \
            --canonOnly \
            --deleteTmp \
            -o {params.outprefix} 2>> {log}
        
        # post-run cleanup
        pigz -p {threads} {params.outprefix}_clean.log 2>> {log}
        pigz -p {threads} {params.outprefix}_clean.TE.log 2>> {log}
        rm -f {params.outprefix}_clean.fa
        """

# ───────────────────────────────────────────────
# Rule: Convert cleaned SAM to BAM + index
# ───────────────────────────────────────────────
rule transcriptclean_sam_to_bam:
    message: "Convert cleaned SAM to BAM for {wildcards.sample}"
    input:
        sam = "02_transcriptclean/{sample}/{sample}_clean.sam"
    output:
        bam = temp("02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam"),
        bai = temp("02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam.bai")
    log:
        "logs/02_transcriptclean/{sample}_sam_to_bam.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_sam_to_bam.txt"
    threads: 8
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """

# ───────────────────────────────────────────────
# Rule: Downsample reads on chrM to avoid slowdowns
# ───────────────────────────────────────────────
rule downsample_chrM:
    message: "Downsample mitochondrial reads for {wildcards.sample}"
    input:
        bam = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam",
        bai = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam.bai"
    output:
        bam = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam",
        bai = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam.bai"
    params:
        max_chrM_reads = MAXCHRMREADS
    log:
        "logs/02_transcriptclean/{sample}_downsample_chrM.log"
    benchmark:
         "benchmarks/02_transcriptclean/{sample}_downsample_chrM.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/seqtk.yaml"
    shell:
        r"""
        {{
            set -euo pipefail
    
            echo "Starting downsample_chrM for {input.bam}"
    
            mito_name=$(samtools idxstats {input.bam} | cut -f 1 | grep -E '^chrM$|^chrMT$|^MT$|^M$' || true)
    
            if [ -n "$mito_name" ]; then
                echo "Mitochondrial contig detected: $mito_name"
    
                # Extract and downsample mitochondrial reads
                samtools view -b {input.bam} "$mito_name" > chrM.bam
                samtools view chrM.bam | seqtk sample -s42 - {params.max_chrM_reads} | \
                    samtools view -Sb - > chrM_downsampled.bam
    
                # Extract all non-mitochondrial reads
                samtools view -h {input.bam} | awk -v mito="$mito_name" '$3 != mito || $1 ~ /^@/' | \
                    samtools view -Sb - > no_chrM.bam
    
                # Merge and sort
                samtools merge -@ {threads} -f merged.bam no_chrM.bam chrM_downsampled.bam
                samtools sort -@ {threads} -o {output.bam} merged.bam
                samtools index {output.bam}
    
                rm chrM.bam chrM_downsampled.bam no_chrM.bam merged.bam
    
            else
                echo "No mitochondrial contig found. Copying original BAM."
                cp {input.bam} {output.bam}
                samtools index {output.bam}
            fi
    
            echo "Finished downsample_chrM for {wildcards.sample}"
        }} >> {log} 2>&1
        """
