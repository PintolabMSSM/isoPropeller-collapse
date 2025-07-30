rule index_fasta:
    message: "Indexing FASTA file: {input}"
    input:
        fasta = "{base}.fasta"
    output:
        fai   = "{base}.fasta.fai"
    log:
        "logs/convert/{base}_faidx.log"
    benchmark:
        "benchmarks/convert/{base}_faidx.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fai})
        samtools faidx {input.fasta} 2>> {log}
        """

rule gtf_to_gff:
    message: "Converting GTF to GFF: {input}"
    input:
        gtf = "{base}.gtf"
    output:
        gff = "{base}.gff"
    params:
        snakedir = SNAKEDIR
    log:
        "logs/convert/{base}_gtf2gff.log"
    benchmark:
        "benchmarks/convert/{base}_gtf2gff.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gff})
        {params.snakedir}workflow/scripts/gtf2gff.pl {input.gtf} > {output.gff} 2>> {log}
        """

rule gff_to_bed:
    message: "Converting GFF to BED: {input}"
    input:
        gff = "{base}.gff"
    output:
        bed = "{base}.bed"
    params:
        snakedir = SNAKEDIR
    log:
        "logs/convert/{base}_gff2bed.log"
    benchmark:
        "benchmarks/convert/{base}_gff2bed.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        {params.snakedir}workflow/scripts/gff2bed.pl {input.gff} > {output.bed} 2>> {log}
        """
