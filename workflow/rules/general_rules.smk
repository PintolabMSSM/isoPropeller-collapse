# ───────────────────────────────────────────────
# Rule: Index FASTA → FAI
# ───────────────────────────────────────────────
rule index_fasta:
    message: "Indexing FASTA file: {input}"
    input:
        fasta = "{base}.fasta"
    output:
        fai   = "{base}.fasta.fai"
    threads: 1
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fai})
        samtools faidx {input.fasta} 2>> {log}
        """

# ───────────────────────────────────────────────
# Rule: GTF → GFF
# ───────────────────────────────────────────────
rule gtf_to_gff:
    message: "Converting GTF to GFF: {input}"
    input:
        gtf = "{base}.gtf"
    output:
        gff = "{base}.gff"
    params:
        snakedir = SNAKEDIR
    threads: 1
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gff})
        {params.snakedir}workflow/scripts/gtf2gff.pl {input.gtf} > {output.gff} 2>> {log}
        """

# ───────────────────────────────────────────────
# Rule: GFF → BED
# ───────────────────────────────────────────────
rule gff_to_bed:
    message: "Converting GFF to BED: {input}"
    input:
        gff = "{base}.gff"
    output:
        bed = "{base}.bed"
    params:
        snakedir = SNAKEDIR
    threads: 1
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        {params.snakedir}workflow/scripts/gff2bed.pl {input.gff} > {output.bed} 2>> {log}
        """
