# ───────────────────────────────────────────────
# Rule: Convert FLNC BAMs to FASTQ (intermediates)
# ───────────────────────────────────────────────
rule bam_to_fastq:
    message: "Converting BAM to FASTQ: {wildcards.sample}, part {wildcards.part}"
    input:
        lambda wc: FLNC_BAM_PARTS[(wc.sample, int(wc.part))]
    output:
        fastq = temp("01_mapping/{sample}/flnc_parts/flnc_part{part}.fastq.gz")
    log:
        "logs/01_mapping/{sample}_bam_to_fastq_{part}.log"
    threads: 12
    benchmark:
        "benchmarks/01_mapping/{sample}_bam_to_fastq_{part}.txt"
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Converting bam to fastq"
        
        samtools fastq {input} | pigz -p {threads} > {output.fastq}
        
        echo "Finished converting bam to fastq"
        ) &> {log}
        """

# ───────────────────────────────────────────────
# Rule: Merge FASTQ parts
# ───────────────────────────────────────────────
rule merge_flnc_fastqs:
    message: "Merging FASTQ parts for sample {wildcards.sample}"
    input:
        lambda wc: [
            f"01_mapping/{wc.sample}/flnc_parts/flnc_part{i}.fastq.gz"
            for i in range(len(FLNC_BAMS[wc.sample]))
        ]
    output:
        fastq = temp("01_mapping/{sample}/flnc_merged.fastq.gz")
    log:
        "logs/01_mapping/{sample}_merge_fastq.log"
    threads: 1
    benchmark:
        "benchmarks/01_mapping/{sample}_fastq_merge.txt"
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Merging FASTQ parts"
        
        cat {input} > {output.fastq}        
        
        echo "Finished merging FASTQ parts"
        ) &> {log}
        """

# ───────────────────────────────────────────────
# Rule: Mapping with minimap2
# ───────────────────────────────────────────────
rule mapping:
    message: "Mapping FLNC merged FASTQ for sample {wildcards.sample}"
    input:
        fastq = "01_mapping/{sample}/flnc_merged.fastq.gz",
        ref   = GENOMEFASTA
    output:
        bam   = temp("01_mapping/{sample}/mapped.bam"),
        bai   = temp("01_mapping/{sample}/mapped.bam.bai")
    log:
        "logs/01_mapping/{sample}_mapping.log"
    threads: 12
    params:
        max_intron_length = MAXINTRONLEN,
        minimap_flags     = "-ax splice:hq --MD -uf --secondary=no"
    benchmark: 
        "benchmarks/01_mapping/{sample}_minimap2-mapping.txt"
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Mapping FLNC fastq"
        
        mkdir -p $(dirname {output.bam})
        minimap2 {params.minimap_flags} \
            -G {params.max_intron_length} \
            -t {threads} \
               {input.ref} \
               {input.fastq} \
            | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        
        echo "Finished mapping FLNC fastq"
        ) &> {log}
        """

# ───────────────────────────────────────────────
# Rule: Run talon_label_reads to prepare a talon-compatible mapping file
# ───────────────────────────────────────────────
rule talon_label_reads:
    message: "Running TALON label_reads on sample {wildcards.sample}"
    input:
        bam = "01_mapping/{sample}/mapped.bam",
        ref = GENOMEFASTA,
        fai = GENOMEFASTA + ".fai"
    output:
        sam     = temp("01_mapping/{sample}/{sample}_mapped_labeled.sam"),
        tsv_tmp = temp("01_mapping/{sample}/{sample}_mapped_fa_read_labels.tsv"),
        tsv     = "01_mapping/{sample}/{sample}_mapped_fa_read_labels.tsv.gz",
        bam     = "01_mapping/{sample}/{sample}_mapped_labeled.bam",
        bai     = "01_mapping/{sample}/{sample}_mapped_labeled.bam.bai",
    log:
        "logs/01_mapping/{sample}_talon_label_reads.log"
    threads: 12
    params:
        tmpdir    = "01_mapping/{sample}/talon_tmp",
        outprefix = "01_mapping/{sample}/{sample}_mapped_fa"
    benchmark: 
        "benchmarks/01_mapping/{sample}_talon_label_reads.txt"
    conda:
        SNAKEDIR + "envs/talon.yaml"
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Running TALON label_reads"
        
        mkdir -p {params.tmpdir}

        # Run talon_label_reads function, which produces sam output unfortunately...
        talon_label_reads \
            --f={input.bam} \
            --g={input.ref} \
            --t={threads} \
            --o={params.outprefix} \
            --tmpDir={params.tmpdir} \
            --deleteTmp

        # Convert SAM to BAM and remove the intermediate SAM
        samtools view -@ {threads} -b {output.sam} > {output.bam}
        samtools index {output.bam}
        
        # Compress the label read csv output
        pigz -p {threads} {output.tsv_tmp} > {output.tsv}
        
        echo "Finished running TALON label_reads"
        ) &> {log}
        """
