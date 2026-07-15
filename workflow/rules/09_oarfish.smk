_OARF_P   = MERGEDISOPREFIX
_OARF_SUF = DEPTHSUFFIX          # "depth-gt1" or "all", from config isoform_depth_suffix

# Sample-set folder shared by all stages, e.g. ISOP_depth-gt1
_OARF_SAMPLESET = f"{_OARF_P}_{_OARF_SUF}"

# Per-stage nested filter-tag path (relative to the sample-set folder). Depth
# grows filter (tpm) < defrag (tpm/defrag) < defrag_pruned (tpm/defrag/prune),
# matching each upstream stage's own output folder.
_OARF_STAGE_TAGPATH = {
    "filter":        TPM_FILTERTAG,
    "defrag":        f"{TPM_FILTERTAG}/{DEFRAG_FILTERTAG}",
    "defrag_pruned": f"{TPM_FILTERTAG}/{DEFRAG_FILTERTAG}/{PRUNE_FILTERTAG}",
}

# Scratch subdir under 09_oarfish/<stage>/ — uniform full-depth path so runs with
# different defrag/prune parameters never collide (the <stage> folder already
# separates content, so the extra depth for filter/defrag is harmless).
_OARF_SUBDIR = f"{_OARF_SAMPLESET}/{TPM_FILTERTAG}/{DEFRAG_FILTERTAG}/{PRUNE_FILTERTAG}"

# stage -> {parent, stem}.  `parent` is the stage output folder that holds the
# per-filtertag subdirectory; `stem` is the basename of that stage's primary GTF
# (without extension).  Both depend only on prefix/suffix, never on the filtertag.
OARFISH_STAGE_SPECS = {
    "filter": {
        "parent": "04_isoPropeller-filter",
        "stem":   "isoqc_pass",
    },
    "defrag": {
        "parent": "06_isoPropeller-defrag",
        "stem":   "isoqc_pass_defrag",
    },
    "defrag_pruned": {
        "parent": "07_isoPropeller-defrag-pruned",
        "stem":   "isoqc_pass_defrag_pruned",
    },
}


# ── Path helpers ──────────────────────────────────────────────────────────────
def _oarf_stage_dir(stage):
    """Upstream stage folder holding that stage's outputs, e.g.
    04_isoPropeller-filter/ISOP_depth-gt1/tpm10s3."""
    spec = OARFISH_STAGE_SPECS[stage]
    return f"{spec['parent']}/{_OARF_SAMPLESET}/{_OARF_STAGE_TAGPATH[stage]}"

def _oarf_gtf(stage):
    spec = OARFISH_STAGE_SPECS[stage]
    return f"{_oarf_stage_dir(stage)}/{spec['stem']}.gtf"

def _oarf_count_matrix(stage):
    spec = OARFISH_STAGE_SPECS[stage]
    return f"{_oarf_stage_dir(stage)}/{spec['stem']}_oarfish_counts.tsv"

def _oarf_quant_files(stage):
    return [
        f"09_oarfish/{stage}/{_OARF_SUBDIR}/quant/{s}/{s}.quant" for s in SAMPLES
    ]

# Final deliverables consumed by `rule all` (one count matrix per stage,
# written into each stage's own filtertag folder).
OARFISH_COUNT_MATRICES = [_oarf_count_matrix(st) for st in OARFISH_STAGE_SPECS]

wildcard_constraints:
    # longest-first so "defrag_pruned" is not shadowed by "defrag"
    stage  = "defrag_pruned|defrag|filter",
    sample = r"[^/]+",


# ───────────────────────────────────────────────
# Helper: merge per-sample oarfish .quant files into a count matrix
# ───────────────────────────────────────────────
def _merge_oarfish_counts(quant_paths, sample_names, out_path):
    """Combine per-sample oarfish `.quant` files into a transcript x sample
    matrix of `num_reads`. Transcript order follows the first sample's file.
    """
    import os

    per_sample = {}          # sample -> {tname: value(str)}
    tname_order = []         # preserve order from the first non-empty file
    seen = set()

    for s, path in zip(sample_names, quant_paths):
        counts = {}
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            per_sample[s] = counts
            continue
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            # transcript-name column ("tname" if present, else first column)
            name_idx = header.index("tname") if "tname" in header else 0
            if "num_reads" not in header:
                raise ValueError(f"'num_reads' column not found in {path}; header={header}")
            nr_idx = header.index("num_reads")
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(name_idx, nr_idx):
                    continue
                tname = parts[name_idx]
                counts[tname] = parts[nr_idx]
                if tname not in seen:
                    seen.add(tname)
                    tname_order.append(tname)
        per_sample[s] = counts

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as out:
        out.write("#TranscriptID\t" + "\t".join(sample_names) + "\n")
        for tname in tname_order:
            row = [tname] + [per_sample[s].get(tname, "0") for s in sample_names]
            out.write("\t".join(row) + "\n")


# ───────────────────────────────────────────────
# Rule: Build a transcriptome FASTA from a stage GTF
# ───────────────────────────────────────────────
rule oarfish_build_transcriptome:
    message: "Building transcriptome FASTA for oarfish ({wildcards.stage})"
    input:
        gtf        = lambda wc: _oarf_gtf(wc.stage),
        genome     = GENOMEFASTA,
        genome_fai = GENOMEFASTA + ".fai"
    output:
        fa = f"09_oarfish/{{stage}}/{_OARF_SUBDIR}/transcriptome.fa"
    log:
        f"logs/09_oarfish/{{stage}}/{_OARF_SUBDIR}/build_transcriptome.log"
    benchmark:
        f"benchmarks/09_oarfish/{{stage}}/{_OARF_SUBDIR}/build_transcriptome.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/oarfish.yaml"
    shell:
        r"""
        (
        set -euo pipefail

        echo "Extracting spliced transcript sequences with gffread ({wildcards.stage})"
        mkdir -p "$(dirname "{output.fa}")"

        # -w: write spliced exon (transcript) sequences, named by transcript_id
        # -g: genome FASTA (uses the provided .fai)
        gffread -w "{output.fa}" -g "{input.genome}" "{input.gtf}"

        echo "Finished building transcriptome FASTA"
        ) &> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Build the oarfish/minimap2 index once per stage
# ───────────────────────────────────────────────
rule oarfish_index:
    message: "Indexing transcriptome for oarfish ({wildcards.stage})"
    input:
        fa = f"09_oarfish/{{stage}}/{_OARF_SUBDIR}/transcriptome.fa"
    output:
        mmi = f"09_oarfish/{{stage}}/{_OARF_SUBDIR}/transcriptome.mmi"
    log:
        f"logs/09_oarfish/{{stage}}/{_OARF_SUBDIR}/index.log"
    benchmark:
        f"benchmarks/09_oarfish/{{stage}}/{_OARF_SUBDIR}/index.txt"
    threads: 8
    conda:
        SNAKEDIR + "envs/oarfish.yaml"
    params:
        seq_tech = OARFISH_SEQ_TECH,
        src_flag = OARFISH_TX_SRC_FLAG
    shell:
        r"""
        (
        set -euo pipefail

        echo "Building minimap2 index via oarfish --only-index ({wildcards.stage})"

        oarfish \
            --only-index \
            {params.src_flag} "{input.fa}" \
            --index-out      "{output.mmi}" \
            --seq-tech       "{params.seq_tech}" \
            -j {threads}

        echo "Finished indexing transcriptome"
        ) &> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Quantify one sample against a stage transcriptome
# ───────────────────────────────────────────────
rule oarfish_quant_sample:
    message: "oarfish quant: {wildcards.sample} ({wildcards.stage})"
    input:
        reads = "01_mapping/{sample}/flnc_merged.fastq.gz",
        index = f"09_oarfish/{{stage}}/{_OARF_SUBDIR}/transcriptome.mmi"
    output:
        quant = f"09_oarfish/{{stage}}/{_OARF_SUBDIR}/quant/{{sample}}/{{sample}}.quant",
        meta  = f"09_oarfish/{{stage}}/{_OARF_SUBDIR}/quant/{{sample}}/{{sample}}.meta_info.json"
    log:
        f"logs/09_oarfish/{{stage}}/{_OARF_SUBDIR}/{{sample}}_quant.log"
    benchmark:
        f"benchmarks/09_oarfish/{{stage}}/{_OARF_SUBDIR}/{{sample}}_quant.txt"
    threads: 16
    conda:
        SNAKEDIR + "envs/oarfish.yaml"
    params:
        out_prefix = f"09_oarfish/{{stage}}/{_OARF_SUBDIR}/quant/{{sample}}/{{sample}}",
        seq_tech   = OARFISH_SEQ_TECH,
        extra      = OARFISH_EXTRA_ARGS
    shell:
        r"""
        (
        set -euo pipefail

        echo "Running oarfish (read mode) for {wildcards.sample} against {wildcards.stage} transcriptome"
        mkdir -p "$(dirname "{params.out_prefix}")"

        oarfish \
            --reads    "{input.reads}" \
            --index    "{input.index}" \
            --seq-tech "{params.seq_tech}" \
            -j {threads} \
            -o "{params.out_prefix}" \
            {params.extra}

        echo "Finished oarfish quant"
        ) &> "{log}"
        """


# ───────────────────────────────────────────────
# Rule: Merge per-sample counts into a matrix (one rule per stage so the
# matrix is written into the stage's own filtertag folder)
# ───────────────────────────────────────────────
rule oarfish_counts_filter:
    message: "Merging oarfish counts (filter)"
    input:
        quants = _oarf_quant_files("filter")
    output:
        matrix = _oarf_count_matrix("filter")
    log:
        f"logs/09_oarfish/filter/{_OARF_SUBDIR}/merge_counts.log"
    threads: 1
    run:
        _merge_oarfish_counts(list(input.quants), list(SAMPLES), output.matrix)


rule oarfish_counts_defrag:
    message: "Merging oarfish counts (defrag)"
    input:
        quants = _oarf_quant_files("defrag")
    output:
        matrix = _oarf_count_matrix("defrag")
    log:
        f"logs/09_oarfish/defrag/{_OARF_SUBDIR}/merge_counts.log"
    threads: 1
    run:
        _merge_oarfish_counts(list(input.quants), list(SAMPLES), output.matrix)


rule oarfish_counts_defrag_pruned:
    message: "Merging oarfish counts (defrag_pruned)"
    input:
        quants = _oarf_quant_files("defrag_pruned")
    output:
        matrix = _oarf_count_matrix("defrag_pruned")
    log:
        f"logs/09_oarfish/defrag_pruned/{_OARF_SUBDIR}/merge_counts.log"
    threads: 1
    run:
        _merge_oarfish_counts(list(input.quants), list(SAMPLES), output.matrix)
