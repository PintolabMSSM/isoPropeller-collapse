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
        fastq = temp("01_mapping/{sample}/flnc_merged.fastq.gz")
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
        bam   = temp("01_mapping/{sample}/mapped.bam"),
        bai   = temp("01_mapping/{sample}/mapped.bam.bai")
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
# Rule: Run talon_label_reads to prepare a talon-compatible mapping file
# ───────────────────────────────────────────────
rule talon_label_reads:
    message: "Running TALON label_reads on sample {wildcards.sample}"
    input:
        bam = "01_mapping/{sample}/mapped.bam",
        ref = GENOMEFASTA,
        fai = GENOMEFASTA + ".fai"
    output:
        tsv     = "01_mapping/{sample}/{sample}_mapped_fa_read_labels.tsv.gz",
        bam     = "01_mapping/{sample}/{sample}_mapped_labeled.bam",
        bai     = "01_mapping/{sample}/{sample}_mapped_labeled.bam.bai"
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
        echo "Running TALON label_reads"
        
        mkdir -p {params.tmpdir}

        # Run talon_label_reads function, which produces sam output unfortunately...
        talon_label_reads \
            --f="{input.bam}" \
            --g="{input.ref}" \
            --t={threads} \
            --o="{params.outprefix}" \
            --tmpDir="{params.tmpdir}" \
            --deleteTmp

        # Convert SAM to BAM and remove the intermediate SAM
        samtools view -@ {threads} -b "{params.outprefix}_labeled.sam" \
            | samtools sort -@ {threads} -o "{output.bam}"
        samtools index "{output.bam}"
        
        # Delete the intermediary sam file
        rm -f "{params.outprefix}_labeled.sam"

        # Compress the label read csv output
        pigz -f -c -p {threads} "{params.outprefix}_read_labels.tsv" > "{output.tsv}"
        
        # Delete the intermediary tsv file
        rm -f "{params.outprefix}_read_labels.tsv"
        
        echo "Finished running TALON label_reads"
        ) &> "{log}"
        """
