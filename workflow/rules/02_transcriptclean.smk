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
        r"""
        (
        echo "Converting {input.bam} to sam for transcriptclean"
        
        mkdir -p "02_transcriptclean/{wildcards.sample}"
        samtools view -h "{input.bam}" > "{output.sam}"
        
        echo "Finished creating transcriptclean sam for {wildcards.sample}"
        ) &> "{log}"
        """

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
        cleaned_sam     = temp("02_transcriptclean/{sample}/{sample}_clean.sam"),
        clean_fa        = temp("02_transcriptclean/{sample}/{sample}_clean.fa"),
        clean_log       = temp("02_transcriptclean/{sample}/{sample}_clean.log"),
        clean_te_log    = temp("02_transcriptclean/{sample}/{sample}_clean.TE.log"),
        clean_log_gz    = "02_transcriptclean/{sample}/{sample}_clean.log.gz",
        clean_te_log_gz = "02_transcriptclean/{sample}/{sample}_clean.TE.log.gz"
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
        r"""
        (
        echo "Starting transcriptclean for {input.sam}"

        # Capture transcriptclean's stderr to a temporary file
        # The main output files will be handled separately
        transcriptclean \
            -s "{input.sam}" \
            -g "{input.ref}" \
            -j "{params.junctions_file}" \
            -t {threads} \
            --primaryOnly \
            --canonOnly \
            --deleteTmp \
            -o "{params.outprefix}"

        # Compress the log files and redirect output to the final paths
        pigz -f -c -p {threads} "{output.clean_log}"    > "{output.clean_log_gz}"
        pigz -f -c -p {threads} "{output.clean_te_log}" > "{output.clean_te_log_gz}"

        echo "Finished transcriptclean for {wildcards.sample}"
        ) &> "{log}"
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
        r"""
        (
        echo "Starting conversion of {input.sam} to BAM"
        
        samtools sort -@ {threads} -o "{output.bam}" "{input.sam}"
        samtools index "{output.bam}"
    
        echo "Finished conversion to BAM"
        ) &> "{log}"
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
        bam         = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam",
        bai         = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam.bai",
    params:
        max_chrM_reads = MAXCHRMREADS,
        seed           = 4232
    log:
        "logs/02_transcriptclean/{sample}_downsample_chrM.log"
    benchmark:
         "benchmarks/02_transcriptclean/{sample}_downsample_chrM.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/pysam.yaml"
    script:
        SNAKEDIR + "scripts/downsample_chrM_bam.py"
