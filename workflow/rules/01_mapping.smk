# ───────────────────────────────────────────────
# Rule: Merge FASTQ parts
# ───────────────────────────────────────────────
rule merge_flnc_fastqs:
    message: "Merging FASTQ parts for sample {wildcards.sample}"
    input:
        lambda wc: [
            f"01_mapping/{wc.sample}/flnc_parts/flnc_part{i}.fastq.gz"
            for i in range(len(PARTS[wc.sample]))
        ]
    output:
        fastq = "01_mapping/{sample}/flnc_merged.fastq.gz"
    log:
        "logs/01_mapping/{sample}_merge_fastq.log"
    threads: 1
    benchmark:
        "benchmarks/01_mapping/{sample}_fastq_merge.txt"
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    run:
        import shlex

        if len(input) > 1:
            # If there are multiple input files, concatenate them after quoting filenames for special chars
            num_files     = len(input)
            fastqs_quoted = " ".join([shlex.quote(str(f)) for f in input])

            shell(
                r"""
                (
                echo "Merging {num_files} FASTQ parts"
                
                cat {fastqs_quoted} > "{output.fastq}"
                
                echo "Finished merging FASTQ parts"
                ) &> "{log}"
                """
            )
        else:
            # If there is only one input file, simply move it to the output path.
            shell(
                r"""
                (
                echo "Only one FASTQ part found. Moving '{input[0]}' to the output path."
                
                mv "{input[0]}" "{output.fastq}"
                
                echo "Finished moving file"
                ) &> "{log}"
                """
            )

# ───────────────────────────────────────────────
# Rule: Mapping with minimap2
# ───────────────────────────────────────────────
rule mapping:
    message: "Mapping FLNC merged FASTQ for sample {wildcards.sample}"
    input:
        fastq = "01_mapping/{sample}/flnc_merged.fastq.gz",
        ref   = GENOMEFASTA
    output:
        bam   = "01_mapping/{sample}/{sample}_mapped.bam",
        bai   = "01_mapping/{sample}/{sample}_mapped.bam.bai"
    log:
        "logs/01_mapping/{sample}_mapping.log"
    threads: 24
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
        echo "Mapping FLNC fastq"
        
        mkdir -p $(dirname {output.bam})
        minimap2 {params.minimap_flags} \
            -G {params.max_intron_length} \
            -t {threads} \
               "{input.ref}" \
               "{input.fastq}" \
            | samtools sort -@ {threads} -o "{output.bam}"
        samtools index "{output.bam}"
        
        echo "Finished mapping FLNC fastq"
        ) &> "{log}"
        """

# ───────────────────────────────────────────────
# Rule: Helper function to determine mapping output file path
# ───────────────────────────────────────────────
def get_final_mapping_bam(wildcards):
    """
    Returns downsampled BAM path if MAXCHRMREADS > -1, otherwise returns mapped.bam.
    """
    if MAXCHRMREADS >= 0:
        return f"01_mapping/{wildcards.sample}/{wildcards.sample}_mapped_chrM_DS.bam"
    return f"01_mapping/{wildcards.sample}/{wildcards.sample}_mapped.bam"


# ───────────────────────────────────────────────
# Rule: Downsample reads on chrM to avoid slowdowns in subsequent correction steps
# ───────────────────────────────────────────────
if MAXCHRMREADS >= 0:
    rule downsample_chrM:
        message: "Downsample mitochondrial reads for {wildcards.sample}"
        input:
            bam     = "01_mapping/{sample}/{sample}_mapped.bam",
            bai     = "01_mapping/{sample}/{sample}_mapped.bam.bai"
        output:
            bam     = "01_mapping/{sample}/{sample}_mapped_chrM_DS.bam",
            bai     = "01_mapping/{sample}/{sample}_mapped_chrM_DS.bam.bai"
        params:
            max_chrM_reads = MAXCHRMREADS,
            seed           = 4232
        log:
            "logs/01_mapping/{sample}_downsample_chrM.log"
        benchmark:
            "benchmarks/01_mapping/{sample}_downsample_chrM.txt"
        threads: 24
        conda:
            SNAKEDIR + "envs/pysam.yaml"
        script:
            SNAKEDIR + "scripts/downsample_chrM_bam.py"
