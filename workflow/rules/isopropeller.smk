# Rule: Run isoPropeller
rule run_isopropeller:
    message: "Running isoPropeller on sample {wildcards.sample}"
    input:
        bam  = lambda wc: ( f"02_transcriptclean/{wc.sample}/{wc.sample}_mapped_labeled_tclean.bam"
                            if USE_TC
                            else f"01_mapping/{wc.sample}/{wc.sample}_mapped_labeled.bam"),
        ref  = GENOMEFASTA,
        cage = CAGEBED
    output:
        gtf   = "03_isoPropeller/{sample}/{sample}.gtf"
    log:
        "logs/03_isoPropeller/{sample}_run.log"
    benchmark:
        "benchmarks/03_isoPropeller/{sample}_run.txt"
    threads: 24
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    params:
        sample = "{sample}",
        outdir = "03_isoPropeller/{sample}"
    shell:
        """
        mkdir -p {params.outdir}
        isoPropeller \
            -i {input.bam} \
            -e \
            -p ISOP_{params.sample} \
            -o {params.outdir}/{params.sample} \
            -g {input.ref} \
            -f {input.cage} \
            -t {threads} 2>> {log}
        """
