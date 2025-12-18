
# ───────────────────────────────────────────────
# Prune low-expressed isoforms from isoform clusters
# ───────────────────────────────────────────────
ruleorder: isopropeller_defrag_prune > gff_to_bed
rule isopropeller_defrag_prune:
    message: "Defragmenting isoPropeller outputs"
    input:
        gtf = "07_isoPropeller-defrag/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag.gtf",
        mod = "07_isoPropeller-defrag/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_modal_ends.gtf",
        exp = "07_isoPropeller-defrag/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_exp.txt",
        ids = "07_isoPropeller-defrag/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_id.txt",
        tss = "07_isoPropeller-defrag/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_tss.bed",
        tts = "07_isoPropeller-defrag/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_tts.bed",
        trk = "07_isoPropeller-defrag/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag.trackgroups",
    output:
        gtf = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned.gtf",
        mod = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned_modal_ends.gtf",
        exp = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned_exp.txt",
        ids = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned_id.txt",
        tss = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned_tss.bed",
        tts = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned_tts.bed",
        trk = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned.trackgroups",
        cls = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned_clusters.txt",
        drp = "08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_defrag_pruned_dropped.txt",
    log:
        "logs/08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/prune_{prefix}_{suffix}.txt"
    benchmark:
        "benchmarks/08_isoPropeller-defrag-pruned/{prefix}_{suffix}_{filtertag}/prune_{prefix}_{suffix}.txt"
    threads: 2
    conda:
        SNAKEDIR + "envs/isopropeller.yaml"
    params:
        prefix_val  = MERGEDISOPREFIX,
        prune_py    = SNAKEDIR + "scripts/prune_rare_isoforms_from_clusters.py",
        retain_pct  = PRUNE_LOW_EXPRESSED_ISOFORMS_RETAIN_PCT,
        min_samples = PRUNE_LOW_EXPRESSED_ISFORMS_MIN_SAMPLES,
        match_mode  = PRUNE_LOW_EXPRESSED_ISFORMS_MATCH_MODE,
    shell:
        r"""
        (
        set -euo pipefail

        echo "[prune] Prune low-expressed isoforms in clusters + GTF filtering"

        # Run consolidation; output expression is **TSV**
        python "{params.prune_py}" \
          --gtf  "{input.gtf}" \
          --expr "{input.exp}" \
          --tx-col '#TranscriptID' \
          --expr-col '*' \
          --match-mode "{params.match_mode}" \
          --sample-support-filter \
          --retain-locus-expr-pct-per-sample {params.retain_pct} \
          --min-support-samples {params.min_samples} \
          --min-keep 1 \
          --out filtered.tsv \
          --clusters-out "{output.cls}" \
          --dropped-out "{output.drp}"
        
        # Filter the 
        gtf-filter-attributes.pl -m "{output.exp}" "{input.gtf}"      > "{output.gtf}"
        gtf-filter-attributes.pl -m "{output.exp}" "{input.mod}"      > "{output.mod}"
        intersect-by-ids -ff "{input.ids}" -fc 4 -if "{output.exp}"   > "{output.ids}"
        intersect-by-ids -ff "{input.tss}" -fc 4 -if "{output.exp}"   > "{output.tss}"
        intersect-by-ids -ff "{input.tts}" -fc 4 -if "{output.exp}"   > "{output.tts}"
        cp "{input.trk}" "{output.trk}"

        echo "[prune] Done"

        ) &> "{log}"
        """

