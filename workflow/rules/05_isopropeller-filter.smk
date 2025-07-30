# Define the filter tag used to mark subdirectories and group outputs
FILTERTAG = f"tpm{FILT_TPM_MIN_COUNT}f{FILT_TPM_MIN_FRACT}"

# ───────────────────────────────────────────────
# Define active filters (and their filter rules)
# ───────────────────────────────────────────────
FILTER_FAIL_ID_PATHS = []

if REMOVE_MONOEXONS:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_monoexon_tss/isoqc_fail_{wc.prefix}_{wc.suffix}_monoexon-no-reftss-overlap.ids"
    )

# ───────────────────────────────────────────────
# Rule: Filter Monoexon Transcripts by TSS 
# ───────────────────────────────────────────────
rule filter_monoexon_tss_overlap:
    message: "Filtering monoexons missing TSS overlap ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed  = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        isoform_tss  = "04_isoPropeller-merge/{prefix}_{suffix}_tss.bed",
        reftss_bed   = REFTSS,
        genome_index = GENOMEFASTA + ".fai"
    output:
        fail_bed = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_tss/isoqc_fail_{prefix}_{suffix}_monoexon-no-reftss-overlap.bed",
        fail_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_tss/isoqc_fail_{prefix}_{suffix}_monoexon-no-reftss-overlap.ids"
    params:
        maxdist   = FILT_MONOEXON_TSS_MAXDIST_BP,
        filtertag = FILTERTAG,
        snakedir  = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_tss/filter_monoexon_tss_overlap.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_tss/filter_monoexon_tss_overlap.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fail_bed})
        python {params.snakedir}scripts/filter_monoexon-tss-overlap.py \
            --isoform_bed {input.isoform_bed} \
            --isoform_tss_bed {input.isoform_tss} \
            --reftss_bed {input.reftss_bed} \
            --genome_index {input.genome_index} \
            --max_distance {params.maxdist} \
            --out_bed {output.fail_bed} \
            --out_ids {output.fail_ids} \
            2>> {log}
        """

# ───────────────────────────────────────────────
# Rule: Aggregate Filtered Final Outputs 
# ───────────────────────────────────────────────
rule filter_aggregate_final_outputs:
    message: "Aggregating final isoform outputs after filtering ({wildcards.prefix}_{wildcards.suffix})"
    input:
        gtf = "04_isoPropeller-merge/{prefix}_{suffix}.gtf",
        exp = "04_isoPropeller-merge/{prefix}_{suffix}_exp.txt",
        ids = "04_isoPropeller-merge/{prefix}_{suffix}_id.txt",
        tss = "04_isoPropeller-merge/{prefix}_{suffix}_tss.bed",
        tts = "04_isoPropeller-merge/{prefix}_{suffix}_tts.bed",
        fail_ids = lambda wc: [func(wc) for func in FILTER_FAIL_ID_PATHS]
    output:
        gtf = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass.gtf",
        exp = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_exp.txt",
        ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_id.txt",
        tss = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_tss.bed",
        tts = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_tts.bed",
        qcf = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_fail.ids"
    params:
        snakedir = SNAKEDIR,
        filtertag = FILTERTAG
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/aggregate_outputs.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/aggregate_outputs.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    run:
        from pathlib import Path
        Path(output.gtf).parent.mkdir(parents=True, exist_ok=True)

        if not input.fail_ids:
            shell("""
                cp {input.gtf} {output.gtf} 2>> {log}
                cp {input.exp} {output.exp} 2>> {log}
                cp {input.ids} {output.ids} 2>> {log}
                cp {input.tss} {output.tss} 2>> {log}
                cp {input.tts} {output.tts} 2>> {log}
            """)
        else:
            shell("""
                cat {input.fail_ids} | sort | uniq                       > {output.qcf} 2>> {log}
                gtf-filter-attributes.pl -m {output.qcf} -v {input.gtf}  > {output.gtf} 2>> {log}
                diff-by-ids -ff {input.exp} -if {output.qcf} -fc 1       > {output.exp} 2>> {log}
                diff-by-ids -ff {input.ids} -if {output.qcf} -fc 1       > {output.ids} 2>> {log}
                diff-by-ids -ff {input.tss} -if {output.qcf} -fc 4       > {output.tss} 2>> {log}
                diff-by-ids -ff {input.tts} -if {output.qcf} -fc 4       > {output.tts} 2>> {log}
            """)
