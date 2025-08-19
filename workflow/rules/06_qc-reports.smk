# ───────────────────────────────────────────────
# Rule: Genotyping QC (samtools mpileup + VarScan)
# ───────────────────────────────────────────────
rule qc_genotyping_snps:
    message:
        "Genotyping QC (mpileup→VarScan) for sample {wildcards.sample}"
    input:
        bam = "01_mapping/{sample}/{sample}_mapped_labeled.bam",
        bai = "01_mapping/{sample}/{sample}_mapped_labeled.bam.bai",
        ref = GENOMEFASTA
    output:
        snp = "06_qc-reports/genotyping/{sample}.varscan.snp.tsv"
    log:
        "logs/06_qc-reports/{sample}_genotyping_varscan.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_genotyping_varscan.txt"
    threads: 4
    params:
        min_cov     = VARSCAN_MIN_COV,
        min_reads2  = VARSCAN_MIN_READS2,
    conda:
        SNAKEDIR + "envs/varscan.yaml"
    shell:
        r"""
        set -o pipefail
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
