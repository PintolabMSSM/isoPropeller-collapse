
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
    benchmark: 
        "benchmarks/09_pbfusion/{sample}_pbmm2-mapping.txt"
    conda:
        SNAKEDIR + "envs/pbfusion.yaml"
    params:
        read_group = lambda wildcards: (
            f"@RG\\t"
            f"ID:{wildcards.sample}\\t"
            f"PL:PACBIO\\t"
            f"DS:READTYPE=CCS\\t"
            f"PM:SEQUELII\\t"
            f"SM:{wildcards.sample}"
        )
    shell:
        r"""
        (
          set -euo pipefail

          echo "Mapping Iso-Seq / polished transcripts with pbmm2 (ISOSEQ preset, --sort)"

          mkdir -p "$(dirname "{output.bam}")"

          pbmm2 align "{input.ref}" "{input.fastq}" "{output.bam}" \
            --preset ISOSEQ \
            --sort \
            --rg '{params.read_group}' \
            -j {threads}

          samtools index "{output.bam}"

          echo "Finished mapping with pbmm2"
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
        out_prefix                 = "09_pbfusion/{sample}/{sample}",
        min_fusion_quality         = PBFUSION_MIN_FUSION_QUALITY,
        min_coverage               = PBFUSION_MIN_COVERAGE,
        min_mean_identity          = PBFUSION_MIN_MEAN_IDENTITY,
        min_mean_mapq              = PBFUSION_MIN_MEAN_MAPQ,
        min_fusion_read_fraction   = PBFUSION_MIN_FUSION_READ_FRACTION,
        max_variability            = PBFUSION_MAX_VARIABILITY,
        max_readthrough            = PBFUSION_MAX_READTHROUGH,
        max_genes_in_event         = PBFUSION_MAX_GENES_IN_EVENT,
        min_fusion_fraction        = PBFUSION_MIN_FUSION_FRACTION,
        prom_filter                = PBFUSION_PROM_FILTER,
    shell:
        r"""
        (
        set -euo pipefail
        
        echo "Running pbfusion"
        
        mkdir -p "$(dirname "{params.out_prefix}")"
        
        pbfusion discover \
           --threads                  {threads} \
           --gtf                      "{input.cache}" \
           --output-prefix            "{params.out_prefix}" \
           --min-fusion-quality       "{params.min_fusion_quality}" \
           --min-coverage              {params.min_coverage} \
           --min-mean-identity         {params.min_mean_identity} \
           --min-mean-mapq             {params.min_mean_mapq} \
           --min-fusion-read-fraction  {params.min_fusion_read_fraction} \
           --max-variability           {params.max_variability} \
           --max-readthrough           {params.max_readthrough} \
           --max-genes-in-event        {params.max_genes_in_event} \
           --min-fusion-fraction       {params.min_fusion_fraction} \
           --prom-filter               {params.prom_filter} \
           --verbose \
           "{input.bam}"
        
        echo "Finished running pbfusion"
        ) &> "{log}"
        """
