# ───────────────────────────────────────────────
# Rule: Run isoPropeller
# ───────────────────────────────────────────────
rule run_isopropeller:
    message: "Running isoPropeller on sample {wildcards.sample}"
    input:
        bam  = lambda wc: ( f"02_transcriptclean/{wc.sample}/{wc.sample}_mapped_labeled_tclean.bam"
                            if USE_TC
                            else f"01_mapping/{wc.sample}/{wc.sample}_mapped_labeled.bam"),
        genfasta = GENOMEFASTA,
        reftss   = REFTSS
    output:
        gtf   = "03_isoPropeller/{sample}/{sample}_all.gtf",
        end   = "03_isoPropeller/{sample}/{sample}_all_end_dist.txt",
        stat  = "03_isoPropeller/{sample}/{sample}_all_stat.txt",
        log   = "03_isoPropeller/{sample}/{sample}_all.log"
    log:
        "logs/03_isoPropeller/{sample}_run.log"
    benchmark:
        "benchmarks/03_isoPropeller/{sample}_run.txt"
    resources:
        lsf_queue = "premium",
        mem_mb    = 90000,
        runtime   = 1440
    threads: 24
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    params:
        sample     = "{sample}",
        outdir     = "03_isoPropeller/{sample}",
        extra_args = ISOPROPEXTRAARGS
    shell:
        r"""
        (
        set -euo pipefail

        echo "Running isoPropeller"
        
        mkdir -p {params.outdir}
        isoPropeller \
            -i {input.bam} \
            -e \
            -p ISOP_{params.sample} \
            -o {params.outdir}/{params.sample}_all \
            -g {input.genfasta} \
            -f {input.reftss} \
            -t {threads} \
            {params.extra_args}
        
        echo "Finished running isoPropeller"
        ) &> {log}
        """


# ───────────────────────────────────────────────
# Rule: Filter per-sample GTF fiels to keep transcripts with at least 2 or more reads
# ───────────────────────────────────────────────
rule filter_isopropeller_gtf:
    message: "Filtering isoPropeller GTF for transcripts with depth > 1 for {wildcards.sample}"
    input:
        gtf          = "03_isoPropeller/{sample}/{sample}_all.gtf"
    output:
        filtered_gtf = "03_isoPropeller/{sample}/{sample}_depth-gt1.gtf",
        tmp_id_file  = temp("03_isoPropeller/{sample}/{sample}_temp_transcript_id.txt")
    log:
        "logs/03_isoPropeller/{sample}_filter_gtf.log"
    benchmark:
        "benchmarks/03_isoPropeller/{sample}_filter_gtf.txt"
    threads: 1
    params:
        attr        = "transcript_id",
        onlychr     = ONLYCHR
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    shell:
        r"""
        (
        set -euo pipefail

        echo "Filtering isoPropeller outputs"
        
        # Extract transcript IDs with depth > 1
        if [ "{params.onlychr}" = "True" ]; then
            awk 'BEGIN {{IGNORECASE=1}} $3 == "transcript" && $1 ~ /^chr([0-9]+|[XY])$/' {input.gtf} | grep -v 'depth "1"' | cut -d'"' -f4 > {output.tmp_id_file}
        else
            awk 'BEGIN {{IGNORECASE=1}} $3 == "transcript"' {input.gtf} | grep -v 'depth "1"' | cut -d'"' -f4 > {output.tmp_id_file}
        fi

        # Filter GTF with perl script
        select_gtf_by_attribute_list.pl {input.gtf} {output.filtered_gtf} {output.tmp_id_file} {params.attr}
        
        echo "Finished filtering isoPropeller outputs"
        ) &> {log}
        """
