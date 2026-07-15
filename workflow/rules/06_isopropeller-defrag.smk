
# ───────────────────────────────────────────────
# Remove transcript fragments from isopropeller outputs
# ───────────────────────────────────────────────
ruleorder: isopropeller_defrag > gff_to_bed
rule isopropeller_defrag:
    message: "Defragmenting isoPropeller outputs"
    input:
        gtf = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_pass.gtf",
        mod = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_pass_modal_ends.gtf",
        exp = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_pass_exp.txt",
        ids = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_pass_id.txt",
        tss = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_pass_tss.bed",
        tts = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_pass_tts.bed",
        qcf = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_fail.ids",
        trk = "04_isoPropeller-filter/{prefix}_{suffix}/" + TPM_FILTERTAG + "/isoqc_pass.trackgroups",
    output:
        ctm = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_containment_map.tsv",
        red = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag_exp_redist.txt",
        gtf = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag.gtf",
        mod = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag_modal_ends.gtf",
        exp = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag_exp.txt",
        ids = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag_id.txt",
        tss = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag_tss.bed",
        tts = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag_tts.bed",
        trk = "06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/isoqc_pass_defrag.trackgroups",
    log:
        "logs/06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/defrag_{prefix}_{suffix}.txt"
    benchmark:
        "benchmarks/06_isoPropeller-defrag/{prefix}_{suffix}/" + TPM_FILTERTAG + "/" + DEFRAG_FILTERTAG + "/defrag_{prefix}_{suffix}.txt"
    threads: 12
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    params:
        prefix_val  = MERGEDISOPREFIX,
        consolidate = CONSOLIDATE_CONTAINED_SPLICECHAINS,
        redist_py   = SNAKEDIR + "scripts/consolidate_splice_chain_fragments.py",
        min_frac    = CONSOLIDATE_MERGE_MIN_FRAC,
        min_samples = CONSOLIDATE_MERGE_MIN_SAMPLES,
        round_mode  = CONSOLIDATE_ROUND_MODE,
        extra_args  = CONSOLIDATE_EXTRAARGS
    shell:
        r"""
        (
        set -euo pipefail

        echo "[consolidate] Cascade redistribution + GTF filtering"

        # Run consolidation; output expression is **TSV**
        python "{params.redist_py}" \
          --gtf  "{input.gtf}"     \
          --expr "{input.exp}" \
          --tx-col '#TranscriptID' \
          --expr-col '*' \
          --terminal-only \
          --minimal-superset \
          --merge-min-frac {params.min_frac} --merge-min-samples {params.min_samples} \
          --proportional \
          --drop-contained \
          --counts --round-counts {params.round_mode} \
          --threads {threads} \
          --log-every 1000 \
          --out "{output.red}" \
          --map-out "{output.ctm}" \
          {params.extra_args}

        # Filter the 
        gtf-filter-attributes.pl -m "{output.red}" "{input.gtf}"      > "{output.gtf}"
        gtf-filter-attributes.pl -m "{output.red}" "{input.mod}"      > "{output.mod}"
        intersect-by-ids -ff "{input.exp}" -fc 1 -if "{output.red}"   > "{output.exp}"
        intersect-by-ids -ff "{input.ids}" -fc 4 -if "{output.red}"   > "{output.ids}"
        intersect-by-ids -ff "{input.tss}" -fc 4 -if "{output.red}"   > "{output.tss}"
        intersect-by-ids -ff "{input.tts}" -fc 4 -if "{output.red}"   > "{output.tts}"
        cp "{input.trk}" "{output.trk}"

        echo "[consolidate] Done"

        ) &> "{log}"
        """

