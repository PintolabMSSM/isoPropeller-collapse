# Rule: Convert BAM to SAM for transcriptclean
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


# Rule: Run transcriptclean
rule run_transcriptclean:
    message: "Running transcriptclean on {wildcards.sample}"
    input:
        sam = "02_transcriptclean/{sample}/{sample}.sam",
        ref = GENOMEFASTA
    output:
        cleaned_sam = temp("02_transcriptclean/{sample}/{sample}_clean.sam"),
        junctions   = temp("02_transcriptclean/{sample}/{sample}_clean_junctions.bed"),
        notes       = temp("02_transcriptclean/{sample}/{sample}_clean_notes.txt")
    log:
        "logs/02_transcriptclean/{sample}_transcriptclean.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_transcriptclean.txt"
    params:
        junctions_file = SPLICEJUNCTIONS,
        outprefix      = "02_transcriptclean/{sample}/{sample}_clean"
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
        """


# Rule: Convert cleaned SAM to BAM + index
rule transcriptclean_sam_to_bam:
    message: "Convert cleaned SAM to BAM for {wildcards.sample}"
    input:
        sam = "02_transcriptclean/{sample}/{sample}_clean.sam"
    output:
        bam = "02_transcriptclean/{sample}/{sample}_clean.bam",
        bai = "02_transcriptclean/{sample}/{sample}_clean.bam.bai"
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


# Rule: Cleanup and compress final outputs
rule compress_transcriptclean_outputs:
    message: "Compressing transcriptclean outputs for {wildcards.sample}"
    input:
        bed  = "02_transcriptclean/{sample}/{sample}_clean_junctions.bed",
        note = "02_transcriptclean/{sample}/{sample}_clean_notes.txt"
    output:
        bed_gz  = "02_transcriptclean/{sample}/{sample}_clean_junctions.bed.gz",
        note_gz = "02_transcriptclean/{sample}/{sample}_clean_notes.txt.gz"
    log:
        "logs/02_transcriptclean/{sample}_compress_outputs.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_compress_outputs.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        """
        gzip -c {input.bed} > {output.bed_gz}
        gzip -c {input.note} > {output.note_gz}
        """
