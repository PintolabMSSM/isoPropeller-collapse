# ───────────────────────────────────────────────
# Rule: Stage FLNC BAMs to FASTQ (intermediates)
# ───────────────────────────────────────────────
rule stage_flnc_part:
    message: "BAM→FASTQ: {wildcards.sample}, part {wildcards.part}"
    input:
        bam = lambda wc: PART_MAP[(wc.sample, int(wc.part))],
        pbi = lambda wc: PART_MAP[(wc.sample, int(wc.part))] + ".pbi"
    output:
        fastq = temp("01_mapping/{sample}/flnc_parts/flnc_part{part}.fastq.gz")
    benchmark:
        "benchmarks/01_mapping/{sample}_bam_to_fastq_{part}.txt"
    log:
        "logs/01_mapping/{sample}_bam_to_fastq_{part}.log"
    threads: 12
    conda:
        SNAKEDIR + "envs/mapping.yaml"
    shell:
        r'''
        (
          echo "Converting bam to fastq"
          samtools fastq -@ 2 "{input.bam}" | pigz -f -c -p {threads} > "{output.fastq}"
          echo "Done"
        ) &> "{log}"
        '''
