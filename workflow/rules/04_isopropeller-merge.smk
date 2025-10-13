# ───────────────────────────────────────────────
# Rule: Prepare a temporary list with all the gtf input files
# ───────────────────────────────────────────────
rule prepare_gtf_list:
    message: "Preparing GTF list for merge: {wildcards.suffix}"
    input:
        gtfs = lambda wc: expand("03_isoPropeller/{sample}/{sample}_{suffix}.gtf", sample=SAMPLES, suffix=[wc.suffix])
    output:
        gtf_list = temp("04_isoPropeller-merge/gtf_list_{suffix}.txt")
    log:
        "logs/04_isoPropeller-merge/prepare_gtf_list_{suffix}.log"
    benchmark:
        "benchmarks/04_isoPropeller-merge/prepare_gtf_list_{suffix}.txt"
    shell:
        r"""
        (
        echo "Preparing GTF list for merge"
        
        mkdir -p 04_isoPropeller-merge
        
        # Use a for loop to handle any characters in file names, including spaces.
        for f in {input.gtfs}; do
            echo "$f"
        done > "{output.gtf_list}"

        echo "Finished preparing GTF list for merge"
        ) &> "{log}"
        """

# ───────────────────────────────────────────────
# One rule: merge to a temp prefix, then either consolidate or move to final outputs
# ───────────────────────────────────────────────
rule merge_isopropeller_gtfs:
    message: "Merging isoPropeller GTFs"
    input:
        gtf_list = "04_isoPropeller-merge/gtf_list_{suffix}.txt"
    output:
        gtf = "04_isoPropeller-merge/{prefix}_{suffix}.gtf",
        exp = "04_isoPropeller-merge/{prefix}_{suffix}_exp.txt",
        ids = "04_isoPropeller-merge/{prefix}_{suffix}_id.txt"
    log:
        "logs/04_isoPropeller-merge/merge_{prefix}_{suffix}.log"
    benchmark:
        "benchmarks/04_isoPropeller-merge/merge_{prefix}_{suffix}.txt"
    threads: 24
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    params:
        prefix_val  = MERGEDISOPREFIX,
        consolidate = CONSOLIDATE_CONTAINED_SPLICECHAINS,
        redist_py   = SNAKEDIR + "scripts/consolidate_splice_chain_fragments.py",
        tss_bed     = REFTSS,
        protect_win = CONSOLIDATE_PROTECT_WINDOW,
        round_mode  = CONSOLIDATE_ROUND_MODE,
        tmp_prefix  = "04_isoPropeller-merge/{prefix}_{suffix}_before-consolidate-splice-fragments"
    shell:
        r"""
        (
        set -euo pipefail

        echo "[merge] Running isoPropeller_merge to temp prefix"

        isoPropeller_merge \
          -i "{input.gtf_list}" \
          -o "{params.tmp_prefix}" \
          -p "{params.prefix_val}" \
          -e depth \
          -t {threads}

        if [[ "{params.consolidate}" == "True" || "{params.consolidate}" == "true" ]]; then
          echo "[consolidate] Cascade redistribution + GTF filtering"

          # Run consolidation; output expression is **TSV**
          python "{params.redist_py}" \
            --gtf  "{params.tmp_prefix}.gtf"     \
            --expr "{params.tmp_prefix}_exp.txt" \
            --tx-col 'transcript_id' \
            --expr-col '*' \
            --terminal-only \
            --minimal-superset \
            --protect-tss-bed "{params.tss_bed}" \
            --protect-window {params.protect_win} \
            --proportional \
            --drop-contained \
            --counts --round-counts {params.round_mode} \
            --threads {threads} \
            --log-every 1000 \
            --out "{output.exp}.tmp.tsv" \
            --map-out "04_isoPropeller-merge/{params.prefix_val}_{wildcards.suffix}_before-consolidate-splice-fragments_containment_map.tsv"

          # Final expression: TSV -> overwrite declared output name
          mv "{output.exp}.tmp.tsv" "{output.exp}"

          # Final ID list: first column of TSV
          intersect-by-ids -ff "{params.tmp_prefix}_id.txt" -fc 4 -if  "{output.exp}" > "{output.ids}"

          # Filter GTF to surviving transcripts using the consolidated IDs
          echo "[consolidate] Filtering GTF by consolidated IDs"
          gtf-filter-attributes.pl -a transcript_id -m "{output.ids}" "{params.tmp_prefix}.gtf" > "{output.gtf}"

          echo "[consolidate] Done"

        else
          echo "[merge] Consolidation disabled — moving temp files to final outputs"
          mv "{params.tmp_prefix}.gtf"     "{output.gtf}"
          mv "{params.tmp_prefix}_exp.txt" "{output.exp}"
          mv "{params.tmp_prefix}_id.txt"  "{output.ids}"
        fi

        ) &> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Similar to what we did for the GTF list, we also prepare and end dist listing for all input files
# ───────────────────────────────────────────────
rule prepare_end_dist_list:
    message: "Preparing end distribution file list"
    input:
        end_dists = expand("03_isoPropeller/{sample}/{sample}_all_end_dist.txt", sample=SAMPLES)
    output:
        listfile  = temp("04_isoPropeller-merge/temp_end_dist_list.txt")
    log:
        "logs/04_isoPropeller-merge/prepare_end_dist_list.log"
    benchmark:
        "benchmarks/04_isoPropeller-merge/prepare_end_dist_list.txt"
    shell:
        r"""
        (
        echo "Preparing end distribution file list"
        
        mkdir -p 04_isoPropeller-merge
        
        # Use a for loop to handle any characters in file names, including spaces.
        for f in {input.end_dists}; do
            echo "$f"
        done > "{output.listfile}"

        echo "Finished preparing end distribution file list"
        ) &> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: And finally, we use this end dist listing together with the list of transcript IDs we want to retain
# to select the TSS and TTS regions to accompany the main isoform gtf file
# ───────────────────────────────────────────────
ruleorder: analyze_end_regions > gff_to_bed
rule analyze_end_regions:
    message: "Analyzing end regions for suffix {wildcards.suffix}"
    input:
        dist_list = "04_isoPropeller-merge/temp_end_dist_list.txt",
        id_list   = "04_isoPropeller-merge/{prefix}_{suffix}_id.txt"
    output:
        tss = "04_isoPropeller-merge/{prefix}_{suffix}_tss.bed",
        tts = "04_isoPropeller-merge/{prefix}_{suffix}_tts.bed"
    log:
        "logs/04_isoPropeller-merge/analyze_end_regions_{prefix}_{suffix}.log"
    benchmark:
        "benchmarks/04_isoPropeller-merge/analyze_end_regions_{prefix}_{suffix}.txt"
    threads: 24
    params:
        prefix_val = MERGEDISOPREFIX
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    shell:
        r"""
        (
        echo "Analyzing end regions"
        
        isoPropeller_end_region \
            -i "{input.dist_list}" \
            -o "04_isoPropeller-merge/{params.prefix_val}_{wildcards.suffix}" \
            -d "{input.id_list}" \
            -t {threads}

        echo "Finished analyzing end regions"
        ) &> "{log}"
        """
