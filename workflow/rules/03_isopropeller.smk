# ───────────────────────────────────────────────
# Rule: Run isoPropeller
# ───────────────────────────────────────────────
rule run_isopropeller:
    message: "Running isoPropeller on sample {wildcards.sample}"
    input:
        bam  = lambda wc: ( f"02_transcriptclean/{wc.sample}/{wc.sample}_mapped_labeled_tclean.bam"
                            if USE_TC
                            else f"01_mapping/{wc.sample}/{wc.sample}_mapped_labeled_chrM_DS.bam"),
        genfasta = GENOMEFASTA,
        reftss   = REFTSS
    output:
        gtf   = "03_isoPropeller/{sample}/{sample}_raw.gtf",
        end   = "03_isoPropeller/{sample}/{sample}_raw_end_dist.txt",
        stat  = "03_isoPropeller/{sample}/{sample}_raw_stat.txt",
        log   = "03_isoPropeller/{sample}/{sample}_raw.log"
    log:
        "logs/03_isoPropeller/{sample}_run.log"
    benchmark:
        "benchmarks/03_isoPropeller/{sample}_run.txt"
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
        echo "Running isoPropeller"
        
        mkdir -p "{params.outdir}"
        isoPropeller \
            -i "{input.bam}" \
            -e \
            -p "ISOP_{params.sample}" \
            -o "{params.outdir}/{params.sample}_raw" \
            -g "{input.genfasta}" \
            -f "{input.reftss}" \
            -t {threads} \
            {params.extra_args}
        
        echo "Finished running isoPropeller"
        ) &> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Post-process the raw per-sample isoPropeller GTF into the consumable set
#       selected by config 'isoform_depth_suffix'. Both products apply the same
#       chromosome filtering (when limit_output_to_chrNXY is true); they differ only
#       in whether depth-1 transcripts are removed:
#         - "all"       : keep every transcript (all read depths)
#         - "depth-gt1" : keep only transcripts with read depth > 1
# ───────────────────────────────────────────────
rule postprocess_isopropeller_gtf:
    message: "Post-processing isoPropeller GTF ({wildcards.suffix}) for {wildcards.sample}"
    input:
        gtf          = "03_isoPropeller/{sample}/{sample}_raw.gtf"
    output:
        filtered_gtf = "03_isoPropeller/{sample}/{sample}_{suffix}.gtf",
        tmp_id_file  = temp("03_isoPropeller/{sample}/{sample}_temp_transcript_id_{suffix}.txt")
    log:
        "logs/03_isoPropeller/{sample}_postprocess_gtf_{suffix}.log"
    benchmark:
        "benchmarks/03_isoPropeller/{sample}_postprocess_gtf_{suffix}.txt"
    threads: 1
    wildcard_constraints:
        # Only the two consumable products are built here; the raw GTF is never
        # matched by this rule (it is produced directly by run_isopropeller).
        suffix = "all|depth-gt1"
    params:
        attr        = "transcript_id",
        onlychr     = ONLYCHR
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    shell:
        r"""
        (
        echo "Post-processing isoPropeller outputs ({wildcards.suffix})"

        # For the "depth-gt1" product, drop transcripts with read depth == 1.
        # For the "all" product, keep every transcript (pass-through with cat).
        if [ "{wildcards.suffix}" = "depth-gt1" ]; then
            depth_filter() {{ grep -v 'depth "1"'; }}
        else
            depth_filter() {{ cat; }}
        fi

        # Extract the retained transcript IDs, optionally restricting to chrN/chrX/chrY
        if [ "{params.onlychr}" = "True" ]; then
            awk 'BEGIN {{IGNORECASE=1}} $3 == "transcript" && $1 ~ /^chr([0-9]+|[XY])$/' "{input.gtf}" | depth_filter | cut -d'"' -f4 > "{output.tmp_id_file}"
        else
            awk 'BEGIN {{IGNORECASE=1}} $3 == "transcript"' "{input.gtf}" | depth_filter | cut -d'"' -f4 > "{output.tmp_id_file}"
        fi

        # Filter GTF with perl script
        select_gtf_by_attribute_list.pl "{input.gtf}" "{output.filtered_gtf}" "{output.tmp_id_file}" "{params.attr}"

        echo "Finished post-processing isoPropeller outputs ({wildcards.suffix})"
        ) &> "{log}"
        """
