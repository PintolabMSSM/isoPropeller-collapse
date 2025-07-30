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
        """
        mkdir -p 04_isoPropeller-merge
        printf "%s\n" {input.gtfs} > {output.gtf_list}
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
        """
        isoPropeller_merge \
            -i {input.gtf_list} \
            -o 04_isoPropeller-merge/{params.prefix_val}_{wildcards.suffix} \
            -p {params.prefix_val} \
            -e depth \
            -t {threads} 2>> {log}
        """

# ───────────────────────────────────────────────
# Rule: Similar to what we did for the GTF list, we also prepare and end dist listing for all input files
# ───────────────────────────────────────────────
rule prepare_end_dist_list:
    message: "Preparing end distribution file list"
    input:
        end_dists = expand("03_isoPropeller/{sample}/{sample}_all_end_dist.txt", sample=SAMPLES)
    output:
        listfile = temp("04_isoPropeller-merge/temp_end_dist_list.txt")
    log:
        "logs/04_isoPropeller-merge/prepare_end_dist_list.log"
    benchmark:
        "benchmarks/04_isoPropeller-merge/prepare_end_dist_list.txt"
    shell:
        """
        mkdir -p 04_isoPropeller-merge
        printf "%s\n" {input.end_dists} > {output.listfile}
        """

# ───────────────────────────────────────────────
# Rule: And finally, we use this end dist listing together with the list of transcript IDs we want to retain
# to select the TSS and TTS regions to accompany the main isoform gtf file
# ───────────────────────────────────────────────
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
        """
        isoPropeller_end_region \
            -i {input.dist_list} \
            -o 04_isoPropeller-merge/{params.prefix_val}_{wildcards.suffix} \
            -d {input.id_list} \
            -t {threads} 2>> {log}
        """
