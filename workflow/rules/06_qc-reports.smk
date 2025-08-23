#────────────────────────────────────────────────
# Helper functions
#────────────────────────────────────────────────

def _parts_for(sample):
    return [f"01_mapping/{sample}/flnc_parts/flnc_part{i}.fastq.gz"
            for i, _ in enumerate(PARTS[sample])]

def _all_parts():
    return [f"01_mapping/{s}/flnc_parts/flnc_part{i}.fastq.gz"
            for s in SAMPLES for i, _ in enumerate(PARTS[s])]

def _bam_for(sample):
    return (f"01_mapping/{sample}/{sample}_mapped_labeled.bam"
            if BAM_CHOICE == "labeled"
            else f"01_mapping/{sample}/mapped.bam")

def _bai_for(sample):
    return _bam_for(sample) + ".bai"

def _as_int(x):
    try:
        return int(x if not isinstance(x, (list, tuple)) else x[0])
    except Exception:
        return int(str(x).rstrip(",").strip())

#────────────────────────────────────────────────
# FastQC rules
#────────────────────────────────────────────────

rule fastqc_part:
    message: "FastQC: {wildcards.sample}, part {wildcards.part}"
    input:
        fq = "01_mapping/{sample}/flnc_parts/flnc_part{part}.fastq.gz"
    output:
        zip  = temp("06_qc-reports/fastqc/{sample}/parts/fastqc_part{part}_fastqc.zip"),
        html = temp("06_qc-reports/fastqc/{sample}/parts/fastqc_part{part}_fastqc.html")
    log:
        "logs/06_qc-reports/{sample}_fastqc_part{part}.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_fastqc_part{part}.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r'''
        (
            echo "Running FastQC on {wildcards.sample} part {wildcards.part}"
            
            fastqc --quiet --threads {threads} --outdir "$(dirname "{output.zip}")" "{input.fq}"
            
            inbase=$(basename "{input.fq}" .fastq.gz)

            mv -f "$(dirname "{output.zip}")/"$inbase"_fastqc.zip" "{output.zip}"
            mv -f "$(dirname "{output.zip}")/"$inbase"_fastqc.html" "{output.html}"

        ) &> "{log}"
        '''


#────────────────────────────────────────────────
# Seqkit on FLNC bam
#────────────────────────────────────────────────

rule seqkit_stats_flnc_bam:
    message: "SeqKit stats (cohort) on all merged FASTQs"
    input:
        expand("01_mapping/{sample}/flnc_merged.fastq.gz", sample=SAMPLES)
    output:
        tsv = "06_qc-reports/seqkit/cohort_merged.stats.tsv"
    log:
        "logs/06_qc-reports/cohort_seqkit_merged.log"
    benchmark:
        "benchmarks/06_qc-reports/cohort_seqkit_merged.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r'''
        (
            echo "Running seqkit stats on merged FASTQs"

            tmp=$(mktemp)
            seqkit stats -a -T {input} > "$tmp"

            awk -F'\t' 'BEGIN{{OFS="\t"}}
                NR==1 {{$1="sample"; print; next}}
                {{n=split($1, a, "/"); $1=a[n-1]; print}}' "$tmp" > "{output.tsv}"
            rm -f "$tmp"

        ) &> "{log}"
        '''

rule seqkit_stats_flnc_bam_parts:
    message: "SeqKit stats (cohort) on all FASTQ parts"
    input:
        lambda wc: _all_parts()
    output:
        tsv = "06_qc-reports/seqkit/cohort_parts.stats.tsv"
    log:
        "logs/06_qc-reports/cohort_seqkit_parts.log"
    benchmark:
        "benchmarks/06_qc-reports/cohort_seqkit_parts.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r'''
        (
            echo "Running seqkit stats on all FASTQ parts"
            
            tmp=$(mktemp)
            seqkit stats -a -T {input} > "$tmp"

            awk -F'\t' 'BEGIN{{OFS="\t"}}
                NR==1 {{print "sample","part",$0; next}}
                {{
                    n=split($1, a, "/");
                    s=a[n-2];
                    fn=a[n];
                    match(fn,/flnc_part([0-9]+)\.fastq(\.gz)?/,m);
                    p=(m[1] != "" ? m[1] : "NA");
                    print s, p, $0
                }}' "$tmp" > "{output.tsv}"
            rm -f "$tmp"
        ) &> "{log}"
        '''


#────────────────────────────────────────────────
# RNA-SeQC v2 (per-sample; then cohort aggregate)
#────────────────────────────────────────────────

rule rnaseqc_sample:
    message: "RNA-SeQC: {wildcards.sample} (bam_choice={BAM_CHOICE})"
    input:
        bam = lambda wc: _bam_for(wc.sample),
        bai = lambda wc: _bai_for(wc.sample),
        gtf = REFGTF,
        ref = GENOMEFASTA
    output:
        metrics = "06_qc-reports/rnaseqc/{sample}/metrics.tsv"
    log:
        "logs/06_qc-reports/{sample}_rnaseqc.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_rnaseqc.txt"
    threads: 8
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    params:
        outdir = "06_qc-reports/rnaseqc/{sample}"
    shell:
        r'''
        (
            echo "Running RNA-SeQC for {wildcards.sample}"

            rnaseqc "{input.gtf}" "{input.bam}" "{params.outdir}" --reference "{input.ref}" --threads {threads}

        ) &> "{log}"
        '''

rule rnaseqc_metrics_cohort:
    message: "Aggregate RNA-SeQC metrics (cohort)"
    input:
        expand("06_qc-reports/rnaseqc/{sample}/metrics.tsv", sample=SAMPLES)
    output:
        tsv = "06_qc-reports/rnaseqc/cohort_metrics.tsv"
    log:
        "logs/06_qc-reports/cohort_rnaseqc_metrics.log"
    benchmark:
        "benchmarks/06_qc-reports/cohort_rnaseqc_metrics.txt"
    threads: 1
    shell:
        r'''
        (
            echo "Aggregating RNA-SeQC metrics across cohort"

            first=1
            for f in {input}
            do
                s=$(basename $(dirname "$f"))
                if [ $first -eq 1 ]; then
                    awk -v S="$s" 'NR==1{{print $0"\tsample"; next}} {{print $0"\t"S}}' "$f" > "{output.tsv}"
                    first=0
                else
                    awk -v S="$s" 'NR>1{{print $0"\t"S}}' "$f" >> "{output.tsv}"
                fi
            done

        ) &> "{log}"
        '''

#────────────────────────────────────────────────
# BAM QC (per-sample) + cohort summary
#────────────────────────────────────────────────

rule bam_flagstat:
    message: "samtools flagstat: {wildcards.sample} ({BAM_CHOICE})"
    input:
        bam = lambda wc: _bam_for(wc.sample),
        bai = lambda wc: _bai_for(wc.sample)
    output:
        txt = "06_qc-reports/bamqc/{sample}/flagstat.txt"
    log:
        "logs/06_qc-reports/{sample}_flagstat.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_flagstat.txt"
    threads: 2
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r'''
        (
            echo "Running samtools flagstat for {wildcards.sample}"

            samtools flagstat -@ {threads} "{input.bam}" > "{output.txt}"

        ) &> "{log}"
        '''

rule bam_idxstats:
    message: "samtools idxstats: {wildcards.sample} ({BAM_CHOICE})"
    input:
        bam = lambda wc: _bam_for(wc.sample),
        bai = lambda wc: _bai_for(wc.sample)
    output:
        txt = "06_qc-reports/bamqc/{sample}/idxstats.txt"
    log:
        "logs/06_qc-reports/{sample}_idxstats.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_idxstats.txt"
    threads: 2
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r'''
        (
            echo "Running samtools idxstats for {wildcards.sample}"

            samtools idxstats "{input.bam}" > "{output.txt}"

        ) &> "{log}"
        '''

rule bam_stats:
    message: "samtools stats: {wildcards.sample} ({BAM_CHOICE})"
    input:
        bam = lambda wc: _bam_for(wc.sample),
        bai = lambda wc: _bai_for(wc.sample)
    output:
        txt = "06_qc-reports/bamqc/{sample}/stats.txt"
    log:
        "logs/06_qc-reports/{sample}_stats.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_stats.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r'''
        (
            echo "Running samtools stats for {wildcards.sample}"

            samtools stats -@ {threads} "{input.bam}" > "{output.txt}"

        ) &> "{log}"
        '''

rule bam_chrM_count:
    message: "chrM count: {wildcards.sample} (threshold={params.chrm_threshold})"
    input:
        idx = "06_qc-reports/bamqc/{sample}/idxstats.txt"
    output:
        count_txt  = "06_qc-reports/bamqc/{sample}/chrM_count.txt",
        status = "06_qc-reports/bamqc/{sample}/chrM_status.txt"
    log:
        "logs/06_qc-reports/{sample}_chrM_count.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_chrM_count.txt"
    params:
        chrm_threshold = MAXCHRMREADS
    threads: 1
    shell:
        r'''
        (
            echo "Counting chrM reads for {wildcards.sample}"

            mt=$(awk '($1=="chrM"||$1=="MT"||$1=="M"||$1=="ChrM"){{s+=$3}} END{{print (s==""?0:s)}}' "{input.idx}")
            echo "$mt" > "{output.count_txt}"
            
            thr="{params.chrm_threshold}"
            
            if [ "$mt" -le "$thr" ]; then
                echo "PASS ($mt <= $thr)" > "{output.status}"; 
            else
                echo "FAIL ($mt > $thr)" > "{output.status}"; 
            fi
        ) &> "{log}"
        '''

rule bam_qc_summary:
    message: "Aggregate BAM QC summary: {wildcards.sample}"
    input:
        flagstat = "06_qc-reports/bamqc/{sample}/flagstat.txt",
        idxstats = "06_qc-reports/bamqc/{sample}/idxstats.txt",
        stats    = "06_qc-reports/bamqc/{sample}/stats.txt",
        chrM_ct  = "06_qc-reports/bamqc/{sample}/chrM_count.txt",
        chrM_ok  = "06_qc-reports/bamqc/{sample}/chrM_status.txt"
    output:
        tsv = "06_qc-reports/bamqc/{sample}/summary.tsv"
    log:
        "logs/06_qc-reports/{sample}_bam_qc_summary.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_bam_qc_summary.txt"
    params:
        bam_choice     = BAM_CHOICE,
        chrm_threshold = MAXCHRMREADS
    threads: 1
    shell:
        r'''
        (
            echo "Summarizing BAM QC for {wildcards.sample}"

            total=$(awk '/in total/ {{print $1; exit}}' "{input.flagstat}")
            mapped=$(awk '/mapped \(/ && !/supplementary/ {{print $1; exit}}' "{input.flagstat}")
            primary=$(awk -F'\t' '$1=="SN" && $2=="reads mapped:" {{print $3}}' "{input.stats}")
            dup=$(awk -F'\t' '$1=="SN" && $2=="reads duplicated:" {{print $3}}' "{input.stats}")
            pct_mapped=$(awk -v m="$mapped" -v t="$total" 'BEGIN{{if(t>0) printf "%.4f", (m/t)*100; else print "NA"}}')
            chrM=$(cat "{input.chrM_ct}")
            status=$(cat "{input.chrM_ok}")
            
            {{
                echo -e "metric\tvalue"
                echo -e "bam_choice\t{params.bam_choice}"
                echo -e "total_reads\t$total"
                echo -e "mapped_reads\t$mapped"
                echo -e "primary_mapped_reads\t$primary"
                echo -e "duplicates\t$dup"
                echo -e "pct_mapped\t$pct_mapped"
                echo -e "chrM_mapped\t$chrM"
                echo -e "chrM_threshold\t{params.chrm_threshold}"
                echo -e "chrM_status\t$status"
            }} > "{output.tsv}"
        ) &> "{log}"
        '''

rule bam_qc_summary_cohort:
    message: "Aggregate BAM QC summaries (cohort)"
    input:
        expand("06_qc-reports/bamqc/{sample}/summary.tsv", sample=SAMPLES)
    output:
        tsv = "06_qc-reports/bamqc/cohort_summary.tsv"
    log:
        "logs/06_qc-reports/cohort_bam_qc_summary.log"
    benchmark:
        "benchmarks/06_qc-reports/cohort_bam_qc_summary.txt"
    threads: 1
    shell:
        r'''
        (
            echo "Aggregating BAM QC summaries across cohort"

            echo -e "sample\tmetric\tvalue" > "{output.tsv}"
            for f in {input};
            do
                s=$(basename $(dirname "$f"))
                awk -v S="$s" 'NR>1{{print S"\t"$1"\t"$2}}' "$f" >> "{output.tsv}";
            done

        ) &> "{log}"
        '''

#────────────────────────────────────────────────
# Genotyping QC
#────────────────────────────────────────────────

rule qc_genotyping_snps:
    message:
        "Genotyping QC (mpileup→VarScan) for sample {wildcards.sample}"
    input:
        bam = "01_mapping/{sample}/{sample}_mapped_labeled.bam",
        bai = "01_mapping/{sample}/{sample}_mapped_labeled.bam.bai",
        ref = GENOMEFASTA
    output:
        snp = "06_qc-reports/genotyping/{sample}.varscan.snp.tsv"
    log:
        "logs/06_qc-reports/{sample}_genotyping_varscan.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_genotyping_varscan.txt"
    threads: 4
    params:
        min_cov     = _as_int(VARSCAN_MIN_COV),
        min_reads2  = _as_int(VARSCAN_MIN_READS2)
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r"""
        (
            echo "Running samtools mpileup + VarScan for {wildcards.sample}"

            samtools mpileup -f "{input.ref}" "{input.bam}" \
                | awk '$4 != 0' \
                | varscan pileup2snp - \
                    --min-coverage {params.min_cov} \
                    --min-reads2 {params.min_reads2} \
                > "{output.snp}"

            echo "Finished genotyping QC for {wildcards.sample}"

        ) &> "{log}"
        """

#────────────────────────────────────────────────
# MultiQC (single cohort report)
#────────────────────────────────────────────────

rule multiqc_cohort:
    message: "MultiQC (cohort across all samples)"
    input:
        # Depend on the actual FastQC output files for every sample and every part
        fastqc_zips     = lambda wc: [
                                        f"06_qc-reports/fastqc/{s}/parts/fastqc_part{i}_fastqc.zip"
                                        for s in SAMPLES for i in range(len(PARTS.get(s, [])))
                                     ],
        rnaseqc_metrics = expand("06_qc-reports/rnaseqc/{sample}/metrics.tsv", sample=SAMPLES)
    output:
        report_html = "06_qc-reports/multiqc/cohort/multiqc_report.html",
        data        = "06_qc-reports/multiqc/cohort/multiqc_data/multiqc_data.json"
    log:
        "logs/06_qc-reports/cohort_multiqc.log"
    benchmark:
        "benchmarks/06_qc-reports/cohort_multiqc.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    params:
        outdir = "06_qc-reports/multiqc/cohort"
    shell:
        r'''
        (
            echo "Running MultiQC on cohort"

            multiqc --quiet -f -o "{params.outdir}" "06_qc-reports/fastqc" "06_qc-reports/rnaseqc"

        ) &> "{log}"
        '''

#────────────────────────────────────────────────
# Filtlong (per-sample) + cohort aggregation
#────────────────────────────────────────────────

def _filtlong_flags():
    import shlex
    p = []
    kp = config.get("filtlong_keep_percent")
    tb = config.get("filtlong_target_bases")
    ml = config.get("filtlong_min_length")
    lw = config.get("filtlong_length_weight")
    mw = config.get("filtlong_mean_q_weight")
    ww = config.get("filtlong_window_q_weight")
    if kp: p += ["--keep_percent", str(kp)]
    if tb: p += ["--target_bases", str(tb)]
    if ml: p += ["--min_length", str(ml)]
    if lw: p += ["--length_weight", str(lw)]
    if mw: p += ["--mean_q_weight", str(mw)]
    if ww: p += ["--window_q_weight", str(ww)]
    extra = config.get("filtlong_extra", "")
    if extra:
        p += shlex.split(str(extra))
    return " ".join(p)


rule filtlong_filter:
    message: "Filtlong filter: {wildcards.sample}"
    input:
        fq = "01_mapping/{sample}/flnc_merged.fastq.gz"
    output:
        filt_fq = "06_qc-reports/filtlong/{sample}/filtered.fastq.gz",
        report_tsv = "06_qc-reports/filtlong/{sample}/report.tsv",
        logtxt  = "06_qc-reports/filtlong/{sample}/filtlong.log"
    log:
        "logs/06_qc-reports/{sample}_filtlong.log"
    benchmark:
        "benchmarks/06_qc-reports/{sample}_filtlong.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    params:
        flags = _filtlong_flags(),
        outdir = "06_qc-reports/filtlong/{sample}"
    shell:
        r'''
        (
            echo "Running Filtlong on {wildcards.sample}"

            filtlong {params.flags} "{input.fq}" 2> "{output.logtxt}" \
                | gzip -c > "{output.filt_fq}"

            pre=$(mktemp); post=$(mktemp)
            seqkit stats -a -T "{input.fq}"      > "$pre"
            seqkit stats -a -T "{output.filt_fq}" > "$post"
            {{
                echo -e "sample\tstage\tfile\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tN50"
                awk 'NR==1{{next}} {{print S"\tpre\t"F"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' S="{wildcards.sample}" F="{input.fq}"  "$pre"
                awk 'NR==1{{next}} {{print S"\tpost\t"F"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' S="{wildcards.sample}" F="{output.filt_fq}" "$post"
            }} > "{output.report_tsv}"
            rm -f "$pre" "$post"

        ) &> "{log}"
        '''


rule filtlong_cohort_report:
    message: "Aggregate Filtlong reports (cohort)"
    input:
        expand("06_qc-reports/filtlong/{sample}/report.tsv", sample=SAMPLES)
    output:
        tsv = "06_qc-reports/filtlong/cohort_report.tsv"
    log:
        "logs/06_qc-reports/cohort_filtlong_report.log"
    benchmark:
        "benchmarks/06_qc-reports/cohort_filtlong_report.txt"
    threads: 1
    shell:
        r'''
        (
            echo "Aggregating Filtlong reports across cohort"

            first=1
            for f in {input}
            do
                if [ $first -eq 1 ]; then
                    cat "$f" > "{output.tsv}"; first=0
                else
                    awk 'NR>1' "$f" >> "{output.tsv}"
                fi
            done

        ) &> "{log}"
        '''

#────────────────────────────────────────────────
# Isoform Filtering QC — single cohort report
#────────────────────────────────────────────────

# Produces a single TSV with rows in long format:
# table   sample  metric              value   filter
# ------ ------  ---------------- ------- -----------------------------
# filter_counts  -       removed             N       <filter_name>
# filter_totals  -       unique_removed    N       -
# filter_totals  -       original_total    N       -
# filter_totals  -       remaining_pass    N       -
# per_sample     <s>     raw_reads           N       -
# per_sample     <s>     assigned_reads    N       -
# per_sample     <s>     retained_pct      float   -

# Helper to enumerate active filter ID files for a given (prefix, suffix)
def _active_filter_id_paths(prefix, suffix):
    # Reuse exactly the same logic you used to build FILTER_FAIL_ID_PATHS,
    # but as literal strings so Snakemake can expand inputs deterministically.
    paths = []
    if REMOVE_MONOEXONS_NO_TSS:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_monoexon_tss/isoqc_fail_{prefix}_{suffix}_monoexon-no-reftss-overlap.ids")
    if REMOVE_MONOEXON_PRE_MRNAS:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_monoexon_premrna/isoqc_fail_{prefix}_{suffix}_monoexon-likely-premrnas.ids")
    if REMOVE_NONCANONICAL_SPLICE:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_noncanonical_splice/isoqc_fail_{prefix}_{suffix}_multiexonic-noncanonical-splices.ids")
    if REMOVE_TSWITCH_ARTIFACTS:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_template_switch/isoqc_fail_{prefix}_{suffix}_multiexonic-rt-switching.ids")
    if REMOVE_ANTISENSE_SPLICEMATCH:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_antisense_match/isoqc_fail_{prefix}_{suffix}_multiexonic-antisense-splicechain-match.ids")
    if REMOVE_CONTAINED_IN_REPEATS:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_repeat_overlap/isoqc_fail_{prefix}_{suffix}_repeatmasker-overlap.ids")
    if REMOVE_PAR_OVERLAP:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_par_overlap/isoqc_fail_{prefix}_{suffix}_PAR-overlap.ids")
    if REMOVE_BELOW_TPM:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_min_tpm/isoqc_fail_{prefix}_{suffix}_min-TPM.ids")
    if REMOVE_TERMINAL_EXONS_SEGDUP:
        paths.append(f"05_isoPropeller-filter/{prefix}_{suffix}_{FILTERTAG}/filt_terminal_exon_segdup/isoqc_fail_{prefix}_{suffix}_mismapped-terminal-exon-in-segdup.ids")
    return paths

# Wrapper so we can use it inside rule.input lambdas
def _active_filter_id_paths_wc(wc):
    return _active_filter_id_paths(wc.prefix, wc.suffix)

# Main cohort report builder
rule iso_qc_cohort_report:
    message: "Isoform filtering QC report: {wildcards.prefix}_{wildcards.suffix} ({FILTERTAG})"
    input:
        # Original & final isoform lists / expression
        orig_ids = "04_isoPropeller-merge/{prefix}_{suffix}_id.txt",
        pass_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_id.txt",
        pass_exp = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_exp.txt",
        # All active filter fail lists
        fail_ids = _active_filter_id_paths_wc,
        # Cohort raw read counts from SeqKit (built once across all samples)
        seqkit_cohort = "06_qc-reports/seqkit/cohort_merged.stats.tsv"
    output:
        tsv = "06_qc-reports/isofilter/{prefix}_{suffix}_{filtertag}/cohort_iso_qc.tsv"
    log:
        "logs/06_qc-reports/{prefix}_{suffix}_{filtertag}_isofilter_cohort_qc.log"
    benchmark:
        "benchmarks/06_qc-reports/{prefix}_{suffix}_{filtertag}_isofilter_cohort_qc.txt"
    threads: 1
    params:
        filtertag = FILTERTAG
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    run:
        import os, re

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        # ---------- helper: safe line count ----------
        def uniq_count(path):
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                return 0
            # unique count of non-empty ids
            with open(path, "r") as fh:
                ids = {ln.strip() for ln in fh if ln.strip()}
            return len(ids)

        # ---------- per-filter counts ----------
        rows = []
        active_fail_files = list(input.fail_ids)
        # Filter name from path: ".../filt_<FILTERNAME>/isoqc_fail_....ids"
        filt_name_re = re.compile(r"/(filt_[^/]+)/")
        for f in active_fail_files:
            if not os.path.exists(f):
                continue
            m = filt_name_re.search(f)
            filt = m.group(1) if m else "unknown_filter"
            n = uniq_count(f)
            rows.append(("filter_counts", "-", "removed", str(n), filt))

        # ---------- union removed / totals ----------
        original_total = uniq_count(input.orig_ids)
        remaining_pass = uniq_count(input.pass_ids)
        # compute union across all filters
        removed_union = set()
        for f in active_fail_files:
            if not os.path.exists(f):
                continue
            with open(f) as fh:
                for ln in fh:
                    ln = ln.strip()
                    if ln:
                        removed_union.add(ln)
        rows.append(("filter_totals", "-", "unique_removed", str(len(removed_union)), "-"))
        rows.append(("filter_totals", "-", "original_total", str(original_total), "-"))
        rows.append(("filter_totals", "-", "remaining_pass", str(remaining_pass), "-"))

        # ---------- per-sample assigned reads (sum of columns in pass_exp) ----------
        # pass_exp first column is "#TranscriptID", subsequent columns are sample names (per your earlier header rewrite)
        assigned_by_sample = {}
        with open(input.pass_exp) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            samples = header[1:]
            sums = [0]*len(samples)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                # tolerate ragged lines
                for i in range(1, min(len(parts), len(samples)+1)):
                    try:
                        sums[i-1] += float(parts[i])
                    except ValueError:
                        pass
            for i, s in enumerate(samples):
                assigned_by_sample[s] = sums[i]

        # ---------- raw reads per sample from SeqKit cohort table ----------
        # Expect columns: sample, format, type, num_seqs, ...
        raw_by_sample = {}
        with open(input.seqkit_cohort) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            # find column index for 'sample' and 'num_seqs'
            col_idx = {name:i for i,name in enumerate(header)}
            i_sample = col_idx.get("sample", 0)
            i_numseq = col_idx.get("num_seqs", 3)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(i_sample, i_numseq):
                    continue
                s = parts[i_sample]
                try:
                    n = float(parts[i_numseq])
                except ValueError:
                    continue
                raw_by_sample[s] = n

        # ---------- assemble per-sample rows ----------
        for s in sorted(assigned_by_sample.keys()):
            raw = raw_by_sample.get(s, 0.0)
            asg = assigned_by_sample.get(s, 0.0)
            pct = (asg / raw * 100.0) if raw > 0 else 0.0
            rows.append(("per_sample", s, "raw_reads",      f"{int(raw)}", "-"))
            rows.append(("per_sample", s, "assigned_reads", f"{int(asg)}", "-"))
            rows.append(("per_sample", s, "retained_pct",   f"{pct:.4f}", "-"))

        # ---------- write output ----------
        with open(output.tsv, "w") as out:
            out.write("table\tsample\tmetric\tvalue\tfilter\n")
            for r in rows:
                out.write("\t".join(r) + "\n")

