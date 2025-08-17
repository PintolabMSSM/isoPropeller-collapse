# ───────────────────────────────────────────────
# Define active filters (and their filter rules)
# ───────────────────────────────────────────────
FILTER_FAIL_ID_PATHS = []

if REMOVE_MONOEXONS_NO_TSS:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_monoexon_tss/isoqc_fail_{wc.prefix}_{wc.suffix}_monoexon-no-reftss-overlap.ids"
    )
if REMOVE_MONOEXON_PRE_MRNAS:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_monoexon_premrna/isoqc_fail_{wc.prefix}_{wc.suffix}_monoexon-likely-premrnas.ids"
    )

if REMOVE_NONCANONICAL_SPLICE:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_noncanonical_splice/isoqc_fail_{wc.prefix}_{wc.suffix}_multiexonic-noncanonical-splices.ids"
    )

if REMOVE_TSWITCH_ARTIFACTS:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_template_switch/isoqc_fail_{wc.prefix}_{wc.suffix}_multiexonic-rt-switching.ids"
    )

if REMOVE_ANTISENSE_SPLICEMATCH:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_antisense_match/isoqc_fail_{wc.prefix}_{wc.suffix}_multiexonic-antisense-splicechain-match.ids"
    )

if REMOVE_CONTAINED_IN_REPEATS:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_repeat_overlap/isoqc_fail_{wc.prefix}_{wc.suffix}_repeatmasker-overlap.ids"
    )

if REMOVE_PAR_OVERLAP:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_par_overlap/isoqc_fail_{wc.prefix}_{wc.suffix}_PAR-overlap.ids"
    )

if REMOVE_BELOW_TPM:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_min_tpm/isoqc_fail_{wc.prefix}_{wc.suffix}_min-TPM.ids"
    )

if REMOVE_TERMINAL_EXONS_SEGDUP:
    FILTER_FAIL_ID_PATHS.append(
        lambda wc: f"05_isoPropeller-filter/{wc.prefix}_{wc.suffix}_{FILTERTAG}/filt_terminal_exon_segdup/isoqc_fail_{wc.prefix}_{wc.suffix}_mismapped-terminal-exon-in-segdup.ids"
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
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_monoexon-tss-overlap.py" \
            --isoform_bed     "{input.isoform_bed}" \
            --isoform_tss_bed "{input.isoform_tss}" \
            --reftss_bed      "{input.reftss_bed}" \
            --genome_index    "{input.genome_index}" \
            --max_distance    {params.maxdist} \
            --out_bed         "{output.fail_bed}" \
            --out_ids         "{output.fail_ids}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter Monoexon pre-mRNA fragments 
# ───────────────────────────────────────────────
rule filter_monoexon_premrna_fragments:
    message: "Filtering monoexons likely to be pre-mRNA fragments ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        reference_bed = REFGTF.replace(".gtf", ".bed")  # assumes REFGTF is defined
    output:
        fail_bed = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_premrna/isoqc_fail_{prefix}_{suffix}_monoexon-likely-premrnas.bed",
        fail_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_premrna/isoqc_fail_{prefix}_{suffix}_monoexon-likely-premrnas.ids"
    params:
        min_intron_ovlp = FILT_MONOEXON_MIN_INTRON_OVLP_BP,
        filtertag = FILTERTAG,
        snakedir = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_premrna/filter_monoexon_premrna_fragments.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_monoexon_premrna/filter_monoexon_premrna_fragments.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_monoexon-premrna-fragments.py" \
            --isoform_bed12       "{input.isoform_bed}" \
            --reference_bed12     "{input.reference_bed}" \
            --min_intron_overlap  {params.min_intron_ovlp} \
            --out_bed             "{output.fail_bed}" \
            --out_ids             "{output.fail_ids}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter multiexon isoforms with non-canonical splice junctions 
# ───────────────────────────────────────────────
rule filter_noncanonical_splice_junctions:
    message: "Filtering multiexonic isoforms with noncanonical splice sites ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        genome_fasta = GENOMEFASTA
    output:
        fail_bed = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_noncanonical_splice/isoqc_fail_{prefix}_{suffix}_multiexonic-noncanonical-splices.bed",
        fail_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_noncanonical_splice/isoqc_fail_{prefix}_{suffix}_multiexonic-noncanonical-splices.ids",
        motifs    = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_noncanonical_splice/isoqc_fail_{prefix}_{suffix}_multiexonic-noncanonical-splices.motifs.txt"
    params:
        filtertag = FILTERTAG,
        snakedir = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_noncanonical_splice/filter_noncanonical_splice_junctions.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_noncanonical_splice/filter_noncanonical_splice_junctions.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_multiexon-noncanonical-splices.py" \
            --isoform_bed12 "{input.isoform_bed}" \
            --genome_fasta  "{input.genome_fasta}" \
            --out_bed       "{output.fail_bed}" \
            --out_ids       "{output.fail_ids}" \
            --out_motifs    "{output.motifs}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter multiexon isoforms with template-switching artifacts
# ───────────────────────────────────────────────
rule filter_template_switching_artifacts:
    message: "Filtering isoforms with potential template switching artifacts ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        genome_fasta = GENOMEFASTA
    output:
        fail_rts  = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_template_switch/isoqc_fail_{prefix}_{suffix}_multiexonic-rt-switching.repeats.txt",
        fail_ids  = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_template_switch/isoqc_fail_{prefix}_{suffix}_multiexonic-rt-switching.ids",
        fail_bed  = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_template_switch/isoqc_fail_{prefix}_{suffix}_multiexonic-rt-switching.bed"
    params:
        filtertag = FILTERTAG,
        snakedir  = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_template_switch/filter_template_switching_artifacts.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_template_switch/filter_template_switching_artifacts.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_multiexon-rt-switching.py" \
            --isoform_bed12 "{input.isoform_bed}" \
            --genome_fasta  "{input.genome_fasta}" \
            --out_rts_tsv   "{output.fail_rts}" \
            --out_ids       "{output.fail_ids}" \
            --out_bed       "{output.fail_bed}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter multiexon isoforms that have antisense splice chain overlaps with sense transcripts
# ───────────────────────────────────────────────
rule filter_antisense_splicechain_match:
    message: "Filtering antisense isoforms with perfect splicechain match to known transcript ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        reference_bed = REFGTF.replace(".gtf", ".bed")
    output:
        fail_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_antisense_match/isoqc_fail_{prefix}_{suffix}_multiexonic-antisense-splicechain-match.ids",
        fail_bed = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_antisense_match/isoqc_fail_{prefix}_{suffix}_multiexonic-antisense-splicechain-match.bed"
    params:
        filtertag = FILTERTAG,
        snakedir  = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_antisense_match/filter_antisense_splicechain_match.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_antisense_match/filter_antisense_splicechain_match.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_multiexon-antisense-splicechain-match.py" \
            --isoform_bed12   "{input.isoform_bed}" \
            --reference_bed12 "{input.reference_bed}" \
            --out_ids         "{output.fail_ids}" \
            --out_bed         "{output.fail_bed}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter isoforms with near-complete overlap to repeatmasker regions
# ───────────────────────────────────────────────
rule filter_repeat_region_overlap:
    message: "Filtering isoforms overlapping repeatmasker regions ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        repeat_bed  = RMSKBED
    output:
        fail_ids  = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_repeat_overlap/isoqc_fail_{prefix}_{suffix}_repeatmasker-overlap.ids",
        fail_bed  = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_repeat_overlap/isoqc_fail_{prefix}_{suffix}_repeatmasker-overlap.bed",
        fail_stats = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_repeat_overlap/isoqc_fail_{prefix}_{suffix}_repeatmasker-overlap.stats.txt"
    params:
        min_frac = FILT_RMSK_MIN_OVLP_FRACT,
        filtertag = FILTERTAG,
        snakedir  = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_repeat_overlap/filter_repeat_region_overlap.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_repeat_overlap/filter_repeat_region_overlap.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_reference-overlap-on-exons.py" \
            --isoform_bed12         "{input.isoform_bed}" \
            --reference_bed12       "{input.repeat_bed}" \
            --min_overlap_fraction  {params.min_frac} \
            --out_ids               "{output.fail_ids}" \
            --out_bed               "{output.fail_bed}" \
            --out_stats             "{output.fail_stats}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter isoforms that overlap PAR regions
# ───────────────────────────────────────────────
rule filter_par_region_overlap:
    message: "Filtering isoforms overlapping PAR regions ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        par_bed     = PARBED
    output:
        fail_ids   = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_par_overlap/isoqc_fail_{prefix}_{suffix}_PAR-overlap.ids",
        fail_bed   = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_par_overlap/isoqc_fail_{prefix}_{suffix}_PAR-overlap.bed",
        fail_stats = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_par_overlap/isoqc_fail_{prefix}_{suffix}_PAR-overlap.stats.txt"
    params:
        min_frac  = FILT_PAR_MIN_OVLP_FRACT,
        filtertag = FILTERTAG,
        snakedir  = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_par_overlap/filter_par_region_overlap.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_par_overlap/filter_par_region_overlap.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_reference-overlap-on-exons.py" \
            --isoform_bed12         "{input.isoform_bed}" \
            --reference_bed12       "{input.par_bed}" \
            --min_overlap_fraction  {params.min_frac} \
            --out_ids               "{output.fail_ids}" \
            --out_bed               "{output.fail_bed}" \
            --out_stats             "{output.fail_stats}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter isoforms that do not meet defined threshold for minTPM in a certain fraction of samples
# ───────────────────────────────────────────────
rule filter_tpm_expression:
    message: "Filtering isoforms with low expression support (TPM) ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        count_matrix = "04_isoPropeller-merge/{prefix}_{suffix}_exp.txt"
    output:
        fail_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_min_tpm/isoqc_fail_{prefix}_{suffix}_min-TPM.ids",
        fail_bed = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_min_tpm/isoqc_fail_{prefix}_{suffix}_min-TPM.bed"
    params:
        min_tpm    = FILT_TPM_MIN_COUNT,
        min_frac   = FILT_TPM_MIN_FRACT,
        filtertag  = FILTERTAG,
        snakedir   = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_min_tpm/filter_tpm_expression.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_min_tpm/filter_tpm_expression.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/isoform-filter.yaml"
    shell:
        """
        mkdir -p "$(dirname "{output.fail_bed}")"
        python "{params.snakedir}scripts/filter_TPM-fraction.py" \
            --count_matrix         "{input.count_matrix}" \
            --isoform_bed12        "{input.isoform_bed}" \
            --min_tpm              {params.min_tpm} \
            --min_fraction_samples {params.min_frac} \
            --out_ids              "{output.fail_ids}" \
            --out_bed              "{output.fail_bed}" \
            2>> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Filter multi-exonic isoforms with terminal exons mismapped due to segdups 
# ───────────────────────────────────────────────
rule filter_terminal_exons_in_segdup:
    message: "Filtering isoforms with mismapped terminal exons in segmental duplications ({wildcards.prefix}_{wildcards.suffix})"
    input:
        isoform_bed = "04_isoPropeller-merge/{prefix}_{suffix}.bed",
        reference_gtf = REFGTF,
        segdup_bed = SEGDUP
    output:
        fail_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_terminal_exon_segdup/isoqc_fail_{prefix}_{suffix}_mismapped-terminal-exon-in-segdup.ids",
        fail_bed = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_terminal_exon_segdup/isoqc_fail_{prefix}_{suffix}_mismapped-terminal-exon-in-segdup.bed"
    params:
        filtertag = FILTERTAG,
        snakedir = SNAKEDIR
    log:
        "logs/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_terminal_exon_segdup/filter_terminal_exons_in_segdup.log"
    benchmark:
        "benchmarks/05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/filt_terminal_exon_segdup/filter_terminal_exons_in_segdup.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/omics-toolkit.yaml"
    shell:
        r"""
        (
        echo "Filtering isoforms with mismapped terminal exons"
        
        outdir="$(dirname "{output.fail_bed}")"
        mkdir -p "$outdir"

        echo "[INFO] Running bed2intronexongff.pl..."
        bed2intronexongff.pl -v 1 "{input.isoform_bed}" > "$outdir/isoqc_temp_corrected.intronexon.gff"

        echo "[INFO] Running gtf-get-gene-regions.pl..."
        gtf-get-gene-regions.pl "{input.reference_gtf}" > "$outdir/isoqc_temp_reference-gene-regions.gtf"

        for LEVEL in 1 2 3 4; do
            echo "[INFO] Level $LEVEL filtering..."
            "{params.snakedir}scripts/filter_segdup-mismapped-terminal-exons.pl" \
                -i "$outdir/isoqc_temp_corrected.intronexon.gff" \
                -g "$outdir/isoqc_temp_reference-gene-regions.gtf" \
                -s "{input.segdup_bed}" \
                -l $LEVEL \
                > "$outdir/isoqc_temp_terminal-exons-in-segdup_$LEVEL.txt"

            echo "[INFO] Parsing LEVEL $LEVEL results..."
            awk -F $'\t' -v OFS=$'\t' '($2=="no" && $3=="yes" && $6>0 && $7>100000) || ($2=="no" && $8=="yes" && $11>0 && $12>100000) {{print $1}}' \
                "$outdir/isoqc_temp_terminal-exons-in-segdup_$LEVEL.txt" \
                > "$outdir/isoqc_temp_mismapped-terminal-exon-in-segdup_$LEVEL.txt"
        done

        echo "[INFO] Merging .ids files..."
        cat "$outdir"/isoqc_temp_mismapped-terminal-exon-in-segdup_*.txt | sort | uniq \
            > "{output.fail_ids}"
        echo "[INFO] Final .ids output:"
        head "{output.fail_ids}"

        echo "[INFO] Extracting BED regions for failed isoforms..."
        intersect-by-ids -ff "{input.isoform_bed}" -fc 4 -if "{output.fail_ids}" > "{output.fail_bed}"
        echo "[INFO] Final .bed output:"
        head "{output.fail_bed}"

        echo "[INFO] Cleaning up intermediate files..."
        rm -f "$outdir"/isoqc_temp_*

        echo "Finished preparing GTF list for merge"
        ) &> "{log}"
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
        qcf = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_fail.ids",
        trk = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass.trackgroups",
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
        import shlex # Import the shlex module

        Path(output.gtf).parent.mkdir(parents=True, exist_ok=True)

        if not input.fail_ids:
            shell("""
                cp "{input.gtf}" "{output.gtf}" 2>> "{log}"
                cp "{input.ids}" "{output.ids}" 2>> "{log}"
                cp "{input.tss}" "{output.tss}" 2>> "{log}"
                cp "{input.tts}" "{output.tts}" 2>> "{log}"
                
                # Reprocess the header of the expression matrix to be compatible with the isoPropeller-annotate pipeline
                awk 'BEGIN{{OFS="\\t"}} NR==1 {{ $1="#TranscriptID"; for(i=2; i<=NF; i++) {{ split($i, parts, "/"); $i=parts[2] }} }} 1' \
                "{input.exp}" > "{output.exp}" 2>> "{log}"
                
                # Generate a default trackgroups file
                awk 'BEGIN{{OFS="\\t"}} NR==1{{for(i=2; i<=NF; i++) print $i, "ALL"; exit}}' "{input.exp}" > "{output.trk}" 2>> "{log}"
                
            """)
        else:
            # Join the list of files into a single, properly quoted string
            # shlex.quote() correctly handles spaces and other special characters
            fail_id_files_quoted = " ".join([shlex.quote(str(f)) for f in input.fail_ids])
            
            shell("""
                sort -u {fail_id_files_quoted}                              > "{output.qcf}" 2>> "{log}"
                gtf-filter-attributes.pl -m "{output.qcf}" -v "{input.gtf}" > "{output.gtf}" 2>> "{log}"
                diff-by-ids -ff "{input.ids}" -if "{output.qcf}" -fc 1      > "{output.ids}" 2>> "{log}"
                diff-by-ids -ff "{input.tss}" -if "{output.qcf}" -fc 4      > "{output.tss}" 2>> "{log}"
                diff-by-ids -ff "{input.tts}" -if "{output.qcf}" -fc 4      > "{output.tts}" 2>> "{log}"
                
                # Reprocess the header of the expression matrix to be compatible with the isoPropeller-annotate pipeline
                diff-by-ids -ff "{input.exp}" -if "{output.qcf}" -fc 1  \
                    | awk 'BEGIN{{OFS="\\t"}} NR==1 {{ $1="#TranscriptID"; for(i=2; i<=NF; i++) {{ split($i, parts, "/"); $i=parts[2] }} }} 1' \
                    > "{output.exp}" 2>> "{log}"

                # Generate a default trackgroups file
                awk 'BEGIN{{OFS="\\t"}} NR==1{{for(i=2; i<=NF; i++) print $i, "ALL"; exit}}' "{input.exp}" > "{output.trk}" 2>> "{log}"

            """)
