
# ───────────────────────────────────────────────
# Rule: Mapping with minimap2 for pbfusion analysis
# ───────────────────────────────────────────────
rule mapping_pbfusion:
    message: "Mapping FLNC merged FASTQ for sample {wildcards.sample}"
    input:
        fastq = "01_mapping/{sample}/flnc_merged.fastq.gz",
        ref   = GENOMEFASTA
    output:
        bam   = temp("09_pbfusion/{sample}_pbfusion.bam"),
        bai   = temp("09_pbfusion/{sample}_pbfusion.bam.bai")
    log:
        "logs/09_pbfusion/{sample}_mapping.log"
    threads: 24
    params:
        max_intron_length = MAXINTRONLEN,
        minimap_flags     = "--eqx -ax splice:hq --MD -uf --secondary=no"
    benchmark: 
        "benchmarks/09_pbfusion/{sample}_minimap2-mapping.txt"
    conda:
        SNAKEDIR + "envs/pbfusion.yaml"
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Mapping FLNC fastq with modified CIGAR string for pbfusion"
        
        mkdir -p "$(dirname {output.bam})"
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
# Rule: Generate a gff-cache for pbfusion
# ───────────────────────────────────────────────
rule pbfusion_cache:
    message: "Generating local gff-cache for pbfusion"
    input:
        gtf = REFGTF
    output:
        cache = "09_pbfusion/reference_cache/gtf.bin"
    log:
        "logs/09_pbfusion/gtf_cache.log"
    benchmark: 
        "benchmarks/09_pbfusion/prepare-reference-cache.txt"
    conda:
        SNAKEDIR + "envs/pbfusion.yaml"
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Preparing gtf cache file for pbfusion"

        pbfusion gff-cache \
          --gtf      "{input.gtf}" \
          --gtf-out  "{output.cache}" \
          --gtf-transcript-allow-lncRNA \
          --verbose
        
        echo "Finished preparing gtf cache file"
        ) &> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Run pbfusion with the gtf cache file
# ───────────────────────────────────────────────
rule run_pbfusion:
    message: "Running pbfusion for sample {wildcards.sample}"
    input:
        bam    = "09_pbfusion/{sample}_pbfusion.bam",
        bai    = "09_pbfusion/{sample}_pbfusion.bam.bai",
        cache  = "09_pbfusion/reference_cache/gtf.bin"
    output:
        out_breakp       = "09_pbfusion/{sample}/{sample}.breakpoints.bed",
        out_transcr      = "09_pbfusion/{sample}/{sample}.transcripts",
        out_breakpgrp    = "09_pbfusion/{sample}/{sample}.breakpoints.groups",
        out_unannot      = "09_pbfusion/{sample}/{sample}.unannotated.bed",
        out_unannot_clst = "09_pbfusion/{sample}/{sample}.unannotated.clusters.bed",
    log:
        "logs/09_pbfusion/{sample}_pbfusion.log"
    threads: 24
    benchmark: 
        "benchmarks/09_pbfusion/{sample}_pbfusion.txt"
    conda:
        SNAKEDIR + "envs/pbfusion.yaml"
    params:
        out_prefix  = "09_pbfusion/{sample}/{sample}"
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Running pbfusion"
        
        mkdir -p "$(dirname "{params.out_prefix}")"
        
        pbfusion \
           --threads        {threads} \
           --gtf            "{input.cache}" \
           --output-prefix  "{params.out_prefix}" \
           --verbose \
           "{input.bam}"
        
        echo "Finished running pbfusion"
        ) &> "{log}"
        """
