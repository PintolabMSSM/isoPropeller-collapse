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
        ref = GENOMEFASTA,
        fai = GENOMEFASTA + ".fai"
    output:
        cleaned_sam   = temp("02_transcriptclean/{sample}/{sample}_clean.sam")
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


# Rule: Convert cleaned SAM to BAM + index
rule transcriptclean_sam_to_bam:
    message: "Convert cleaned SAM to BAM for {wildcards.sample}"
    input:
        sam = "02_transcriptclean/{sample}/{sample}_clean.sam"
    output:
        bam = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam",
        bai = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam.bai"
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
