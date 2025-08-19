# ───────────────────────────────────────────────
# Rule: Genotyping QC (samtools mpileup + VarScan)
# ───────────────────────────────────────────────
rule qc_genotyping_snps:
    message:
        "Genotyping QC (mpileup→VarScan) for sample {wildcards.sample}"
    input:
        bam = lambda wc: (
            f"01_mapping/{wc.sample}/{wc.sample}_mapped_labeled.bam"
            if BAM_CHOICE == "labeled"
            else f"01_mapping/{wc.sample}/mapped.bam"
        ),
        bai = lambda wc: (
            f"01_mapping/{wc.sample}/{wc.sample}_mapped_labeled.bam.bai"
            if BAM_CHOICE == "labeled"
            else f"01_mapping/{wc.sample}/mapped.bam.bai"
        ),
        ref = GENOMEFASTA
    output:
        snp = "06_qc-reports/genotyping/{sample}.varscan.snp.tsv"
    log:
        "logs/06_qc-reports/{sample}_genotyping_varscan.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_genotyping_varscan.txt"
    threads: 4
    params:
        min_cov       = VARSCAN_MIN_COV,
        min_reads2    = VARSCAN_MIN_READS2,
    conda:
        SNAKEDIR + "envs/varscan.yaml"
    shell:
        r"""
        (
        echo "Running samtools mpileup + VarScan for {wildcards.sample}"

        samtools mpileup -f "{input.ref}" "{input.bam}" \
            | awk '$4 != 0' \
            | varscan pileup2snp - \
                --min-coverage {params.min_cov} \
                --min-reads2 {params.min_reads2} \
            > "{output.snp}"

        echo "Finished genotyping QC for {wildcards.sample}"
        ) &> "{log}"
        """
