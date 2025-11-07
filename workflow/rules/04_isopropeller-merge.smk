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
# Rule: Merge the GTF files together into a single output file
# ───────────────────────────────────────────────
rule merge_isopropeller_gtfs:
    message: "Merging isoPropeller GTFs ({wildcards.suffix})"
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
        prefix_val = MERGEDISOPREFIX
    shell:
        r"""
        (
        echo "Merging isoPropeller GTFs"
        
        isoPropeller_merge \
            -i "{input.gtf_list}" \
            -o "04_isoPropeller-merge/{params.prefix_val}_{wildcards.suffix}" \
            -p "{params.prefix_val}" \
            -e depth \
            -t {threads} 
        
        echo "Finished merging isoPropeller GTFs"
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
        id_list   = "04_isoPropeller-merge/{prefix}_{suffix}_id.txt",
        gtf = "04_isoPropeller-merge/{prefix}_{suffix}.gtf",
    output:
        tss       = "04_isoPropeller-merge/{prefix}_{suffix}_tss.bed",
        tts       = "04_isoPropeller-merge/{prefix}_{suffix}_tts.bed",
        tsscount  = "04_isoPropeller-merge/{prefix}_{suffix}_tss_count.txt",
        ttscount  = "04_isoPropeller-merge/{prefix}_{suffix}_tts_count.txt",
        gtf_modal = "04_isoPropeller-merge/{prefix}_{suffix}_modal_ends.gtf",
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
        
        # Generate tts and tss bed files
        isoPropeller_end_region \
            -i "{input.dist_list}" \
            -o "04_isoPropeller-merge/{params.prefix_val}_{wildcards.suffix}" \
            -d "{input.id_list}" \
            -t {threads}

        # TSS clustering and quantification
        isoPropeller_TSS_quantification.pl \
            -i "{input.dist_list}" \
            -o "{output.tsscount}" \
            -d "{input.id_list}" \
            -t {threads}

        # TTS clustering and qantification
        isoPropeller_TTS_quantification.pl \
            -i "{input.dist_list}" \
            -o "{output.ttscount}" \
            -d "{input.id_list}" \
            -t {threads}

        # Update GTF file with modal ends
        isoPropeller_end_update.pl \
            -i "{input.gtf}" \
            -o "{output.gtf_modal}" \
            -a "{output.tsscount}" \
            -b "{output.ttscount}" \
            -t {threads}

        echo "Finished analyzing end regions"
        ) &> "{log}"
        """
