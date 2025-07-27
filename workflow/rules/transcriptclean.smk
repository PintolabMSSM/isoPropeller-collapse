# Pre-processing: reformat FASTA headers or prepare auxiliary files
rule transcriptclean_prep:
    input:
        lambda wildcards: FLNC_PATHS[wildcards.sample]
    output:
        "01_transcriptclean/{sample}/prepared.fasta"
    shell:
        "prepare_fasta.py {input} > {output}"

# Main cleaning step
rule transcriptclean:
    input:
        "01_transcriptclean/{sample}/prepared.fasta"
    output:
        "01_transcriptclean/{sample}/cleaned.fasta"
    shell:
        "transcriptclean --input {input} --output {output}"

# Post-processing: generate final cleaned result
rule transcriptclean_post:
    input:
        "01_transcriptclean/{sample}/cleaned.fasta"
    output:
        "01_transcriptclean/{sample}/final.fasta"
    shell:
        "postprocess_cleaned.py {input} > {output}"
