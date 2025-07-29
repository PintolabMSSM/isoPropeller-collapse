# Rule: Convert FLNC BAMs to FASTA (intermediates)
rule bam_to_fasta:
    message: "Converting BAM to FASTA: {wildcards.sample}, part {wildcards.part}"
    input:
        lambda wc: FLNC_BAM_PARTS[(wc.sample, int(wc.part))]
    output:
        fasta = temp("01_mapping/{sample}/flnc_parts/flnc_part{part}.fasta.gz")
    log:
        "logs/01_mapping/{sample}_bam_to_fasta_{part}.log"
    threads: 12
    benchmark:
        "benchmarks/01_mapping/{sample}_bam_to_fasta_{part}.txt"
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    shell:
        """samtools fasta {input} | pigz -p {threads} > {output.fasta} 2>> {log}"""


# Rule: Merge FASTA parts
rule merge_flnc_fastas:
    message: "Merging FASTA parts for sample {wildcards.sample}"
    input:
        lambda wc: [
            f"01_mapping/{wc.sample}/flnc_parts/flnc_part{i}.fasta.gz"
            for i in range(len(FLNC_BAMS[wc.sample]))
        ]
    output:
        fasta = temp("01_mapping/{sample}/flnc_merged.fasta.gz")
    log:
        "logs/01_mapping/{sample}_merge_fasta.log"
    threads: 1
    benchmark:
        "benchmarks/01_mapping/{sample}_fasta_merge.txt"
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    shell:
        """cat {input} > {output.fasta} 2>> {log}"""


# Rule: Mapping with minimap2
rule mapping:
    message: "Mapping FLNC merged FASTA for sample {wildcards.sample}"
    input:
        fasta = "01_mapping/{sample}/flnc_merged.fasta.gz",
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
        """
        mkdir -p $(dirname {output.bam})
        minimap2 {params.minimap_flags} -G {params.max_intron_length} -t {threads} {input.ref} {input.fasta} 2>> {log} | \
            samtools sort -@ {threads} -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """


# Rule: Create .fai index for genome fasta
rule index_genome_fasta:
    message: "Indexing genome FASTA: {input.fasta}"
    input:
        fasta = GENOMEFASTA
    output:
        fai   = GENOMEFASTA + ".fai"
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        "samtools faidx {input.fasta}"


# Rule: Run talon_label_reads to prepare a talon-compatible mapping file
rule talon_label_reads:
    message: "Running TALON label_reads on sample {wildcards.sample}"
    input:
        bam = "01_mapping/{sample}/mapped.bam",
        ref = GENOMEFASTA,
        fai = GENOMEFASTA + ".fai"
    output:
        bam = "01_mapping/{sample}/{sample}_mapped_labeled.bam",
        bai = "01_mapping/{sample}/{sample}_mapped_labeled.bam.bai"
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
        """
        mkdir -p {params.tmpdir}

        # Run talon_label_reads function, which produces sam output unfortunately...
        talon_label_reads \
            --f={input.bam} \
            --g={input.ref} \
            --t={threads} \
            --o={params.outprefix} \
            --tmpDir={params.tmpdir} \
            --deleteTmp 2>> {log}

        # Convert SAM to BAM and remove the intermediate SAM
        samtools view -@ {threads} -b {params.outprefix}_labeled.sam > {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        pigz -p {threads} {params.outprefix}_read_labels.tsv 2>> {log}
        rm {params.outprefix}_labeled.sam
        """
