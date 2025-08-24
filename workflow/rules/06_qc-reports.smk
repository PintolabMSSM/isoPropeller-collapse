#────────────────────────────────────────────────
# Helper functions
#────────────────────────────────────────────────

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


# ────────────────────────────────────────────────
# FastQC rules
# ────────────────────────────────────────────────

rule fastqc_merged:
    message: "FastQC (merged): {wildcards.sample}"
    input:
        fq = "01_mapping/{sample}/flnc_merged.fastq.gz"
    output:
        fq   = "06_qc-reports/flnc-fastqc/{sample}/{sample}.fastq.gz",
        zip  = "06_qc-reports/flnc-fastqc/{sample}/{sample}_fastqc.zip",
        html = "06_qc-reports/flnc-fastqc/{sample}/{sample}_fastqc.html"
    log:
        "logs/06_qc-reports/flnc-fastqc/{sample}_flnc-fastqc.log"
    benchmark:
        "benchmarks/06_qc-reports/flnc-fastqc/{sample}_flnc-fastqc.txt"
    threads: 1
    params:
        outdir = "06_qc-reports/flnc-fastqc/{sample}"
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    shell:
        r'''
        (
            echo "Running FastQC on merged FASTQ for {wildcards.sample}"

            # Fastqc derives the sample name from the input file
            # We create a symlink here to include the sample name in the fastq.gz file
            ln -s "../../../{input.fq}" "{output.fq}"

            # Run fastqc on the fastq.gz symlink
            fastqc --quiet --threads {threads} --outdir "{params.outdir}" "{output.fq}"

        ) &> "{log}"
        '''


# ────────────────────────────────────────────────
# Seqkit on FLNC fastq
# ────────────────────────────────────────────────

rule seqkit_stats_flnc:
    message: "SeqKit stats (cohort) on all merged FASTQs"
    input:
        expand("01_mapping/{sample}/flnc_merged.fastq.gz", sample=SAMPLES)
    output:
        tsv = "06_qc-reports/flnc-seqkit-stats/seqkit_flnc_wide.stats.tsv"
    log:
        "logs/06_qc-reports/flnc-seqkit-stats/flnc-seqkit-stats.log"
    benchmark:
        "benchmarks/06_qc-reports/flnc-seqkit-stats/flnc-seqkit-stats.txt"
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


# ────────────────────────────────────────────────
# RNA-SeQC v2 (per-sample; then cohort aggregate)
# ────────────────────────────────────────────────

rule rnaseqc_collapse_gtf:
    message: "Collapse reference GTF into unique regions for RNA-SeQC"
    input:
        gtf = REFGTF
    output:
        collapsed = "06_qc-reports/mapped-rnaseqc/reference/collapsed_reference.gtf"
    log:
        "logs/06_qc-reports/mapped-rnaseqc/rnaseqc_collapse_gtf.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-rnaseqc/rnaseqc_collapse_gtf.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    params:
        script = SNAKEDIR + "scripts/gtex-pipeline_collapse_annotation.py",
        extra  = config.get("rnaseqc_collapse_extra", "")
    shell:
        r'''
        (
            echo "Collapsing {input.gtf} → {output.collapsed}"

            python "{params.script}" "{input.gtf}" "{output.collapsed}" {params.extra}

        ) &> "{log}"
        '''


rule rnaseqc_sample:
    message: "RNA-SeQC: {wildcards.sample} (bam_choice={BAM_CHOICE})"
    input:
        bam = lambda wc: _bam_for(wc.sample),
        bai = lambda wc: _bai_for(wc.sample),
        gtf = "06_qc-reports/mapped-rnaseqc/reference/collapsed_reference.gtf",
        ref = GENOMEFASTA
    output:
        metrics = "06_qc-reports/mapped-rnaseqc/{sample}/{sample}.metrics.tsv"
    log:
        "logs/06_qc-reports/mapped-rnaseqc/{sample}_rnaseqc.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-rnaseqc/{sample}_rnaseqc.txt"
    threads: 2
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    params:
        outdir = "06_qc-reports/mapped-rnaseqc/{sample}",
    shell:
        r'''
        (
            echo "Running RNA-SeQC for {wildcards.sample}"

            rnaseqc "{input.gtf}" "{input.bam}" "{params.outdir}" --fasta "{input.ref}" --sample "{wildcards.sample}" --unpaired --mapping-quality 0 --stranded RF --coverage

        ) &> "{log}"
        '''


rule rnaseqc_metrics_cohort:
    message: "Aggregate RNA-SeQC metrics (cohort)"
    input:
        expand("06_qc-reports/mapped-rnaseqc/{sample}/{sample}.metrics.tsv", sample=SAMPLES)
    output:
        tsv = "06_qc-reports/mapped-rnaseqc/rna_seqc_summary_long.tsv"
    log:
        "logs/06_qc-reports/mapped-rnaseqc/rna_seqc_summary.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-rnaseqc/rna_seqc_summary.txt"
    threads: 1
    shell:
        r'''
        (
            echo "Aggregating RNA-SeQC metrics across cohort"

            echo -e "sample\tmetric\tvalue" > "{output.tsv}"
            for f in {input}
            do
                s=$(basename $(dirname "$f"))
                awk -v S="$s" 'NR>1{{print S"\t"$0}}' "$f" >> "{output.tsv}"
            done

        ) &> "{log}"
        '''


# ────────────────────────────────────────────────
# Picard CollectRnaSeqMetrics
# ────────────────────────────────────────────────

# Picard CollectRnaSeqMetrics per sample (uses refFlat from config)
rule picard_collect_rnaseqmetrics:
    message: "Picard CollectRnaSeqMetrics: {wildcards.sample}"
    input:
        bam     = lambda wc: _bam_for(wc.sample),
        bai     = lambda wc: _bai_for(wc.sample),
        refflat = lambda wc: config["reference_annotations_refflat"]  # <- from config.yaml
    output:
        metrics = "06_qc-reports/mapped-picard-RnaSeqMetrics/{sample}/{sample}.RnaSeqMetrics.txt",
        chart   = "06_qc-reports/mapped-picard-RnaSeqMetrics/{sample}/{sample}.RnaSeqMetrics.pdf"
    log:
        "logs/06_qc-reports/mapped-picard-RnaSeqMetrics/{sample}_picard_rnaseqmetrics.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-picard-RnaSeqMetrics/{sample}_picard_rnaseqmetrics.txt"
    threads: 2
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    params:
        strand = config.get("picard_strand", "SECOND_READ_TRANSCRIPTION_STRAND"),   # NONE | FIRST_READ_TRANSCRIPTION_STRAND | SECOND_READ_TRANSCRIPTION_STRAND
        tmpdir = "06_qc-reports/picard/{sample}/tmp"
    shell:
        r'''
        (
            echo "Running Picard CollectRnaSeqMetrics for {wildcards.sample}"

            picard -Xmx4g CollectRnaSeqMetrics \
                I="{input.bam}" \
                O="{output.metrics}" \
                REF_FLAT="{input.refflat}" \
                STRAND_SPECIFICITY="{params.strand}" \
                CHART_OUTPUT="{output.chart}" \
                VALIDATION_STRINGENCY=LENIENT \
                TMP_DIR="{params.tmpdir}"

        ) &> "{log}"
        '''


# ────────────────────────────────────────────────
# BAM QC (per-sample) + cohort summary
# ────────────────────────────────────────────────

rule bam_flagstat:
    message: "samtools flagstat: {wildcards.sample} ({BAM_CHOICE})"
    input:
        bam = lambda wc: _bam_for(wc.sample),
        bai = lambda wc: _bai_for(wc.sample)
    output:
        txt = "06_qc-reports/mapped-bamqc/{sample}/flagstat.txt"
    log:
        "logs/06_qc-reports/mapped-bamqc/{sample}_flagstat.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-bamqc/{sample}_flagstat.txt"
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
        txt = "06_qc-reports/mapped-bamqc/{sample}/idxstats.txt"
    log:
        "logs/06_qc-reports/mapped-bamqc/{sample}_idxstats.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-bamqc/{sample}_idxstats.txt"
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
        txt = "06_qc-reports/mapped-bamqc/{sample}/stats.txt"
    log:
        "logs/06_qc-reports/mapped-bamqc/{sample}_stats.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-bamqc/{sample}_stats.txt"
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
        idx = "06_qc-reports/mapped-bamqc/{sample}/idxstats.txt"
    output:
        count_txt  = "06_qc-reports/mapped-bamqc/{sample}/chrM_count.txt",
        status     = "06_qc-reports/mapped-bamqc/{sample}/chrM_status.txt"
    log:
        "logs/06_qc-reports/mapped-bamqc/{sample}_chrM_count.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-bamqc/{sample}_chrM_count.txt"
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
        flagstat = "06_qc-reports/mapped-bamqc/{sample}/flagstat.txt",
        idxstats = "06_qc-reports/mapped-bamqc/{sample}/idxstats.txt",
        stats    = "06_qc-reports/mapped-bamqc/{sample}/stats.txt",
        chrM_ct  = "06_qc-reports/mapped-bamqc/{sample}/chrM_count.txt",
        chrM_ok  = "06_qc-reports/mapped-bamqc/{sample}/chrM_status.txt"
    output:
        tsv = "06_qc-reports/mapped-bamqc/{sample}/summary.tsv"
    log:
        "logs/06_qc-reports/mapped-bamqc/{sample}_bam_qc_summary.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-bamqc/{sample}_bam_qc_summary.txt"
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
        expand("06_qc-reports/mapped-bamqc/{sample}/summary.tsv", sample=SAMPLES)
    output:
        tsv  = "06_qc-reports/mapped-bamqc/bam_qc_summary_long.tsv"
    log:
        "logs/06_qc-reports/mapped-bamqc/cohort_bam_qc_summary.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-bamqc/cohort_bam_qc_summary.txt"
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

# ────────────────────────────────────────────────
# Genotyping QC
# ────────────────────────────────────────────────

rule qc_genotyping_snps:
    message:
        "Genotyping QC (mpileup→VarScan) for sample {wildcards.sample}"
    input:
        bam = "01_mapping/{sample}/{sample}_mapped_labeled.bam",
        bai = "01_mapping/{sample}/{sample}_mapped_labeled.bam.bai",
        ref = GENOMEFASTA
    output:
        snp = "06_qc-reports/mapped-snp-genotypes/{sample}.varscan.snp.tsv"
    log:
        "logs/06_qc-reports/mapped-snp-genotypes/{sample}_genotyping_varscan.log"
    benchmark:
        "benchmarks/06_qc-reports/mapped-snp-genotypes/{sample}_genotyping_varscan.txt"
    threads: 4
    params:
        min_cov      = _as_int(VARSCAN_MIN_COV),
        min_reads2   = _as_int(VARSCAN_MIN_READS2)
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

# ───────────────────────────────────────────────
# MultiQC
# ───────────────────────────────────────────────
rule multiqc_cohort:
    message: "MultiQC (cohort across all samples)"
    input:
        fastqc_zips     = expand("06_qc-reports/flnc-fastqc/{sample}/{sample}_fastqc.zip", sample=SAMPLES),
        rnaseqc_metrics = expand("06_qc-reports/mapped-rnaseqc/{sample}/{sample}.metrics.tsv", sample=SAMPLES),
        picard_metrics  = expand("06_qc-reports/mapped-picard-RnaSeqMetrics/{sample}/{sample}.RnaSeqMetrics.txt", sample=SAMPLES)
    output:
        report_html = "06_qc-reports/multiqc/multiqc_report.html",
        data        = "06_qc-reports/multiqc/multiqc_data/multiqc_data.json"
    log:
        "logs/06_qc-reports/multiqc/multiqc.log"
    benchmark:
        "benchmarks/06_qc-reports/multiqc/multiqc.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    params:
        outdir = "06_qc-reports/multiqc"
    shell:
        r'''
        (
            echo "Running MultiQC on cohort"

            multiqc --quiet -f -o "{params.outdir}" "06_qc-reports/flnc-fastqc" "06_qc-reports/mapped-rnaseqc"

        ) &> "{log}"
        '''


# ────────────────────────────────────────────────
# Isoform Filtering QC — single cohort report
# ────────────────────────────────────────────────

def _active_filter_id_paths(prefix, suffix):
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

def _active_filter_id_paths_wc(wc):
    return _active_filter_id_paths(wc.prefix, wc.suffix)

rule iso_qc_cohort_report:
    message: "Isoform filtering QC report: {wildcards.prefix}_{wildcards.suffix} ({FILTERTAG})"
    input:
        orig_ids = "04_isoPropeller-merge/{prefix}_{suffix}_id.txt",
        pass_ids = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_id.txt",
        pass_exp = "05_isoPropeller-filter/{prefix}_{suffix}_{filtertag}/{prefix}_{suffix}_isoqc_pass_exp.txt",
        fail_ids = _active_filter_id_paths_wc,
        seqkit_cohort = "06_qc-reports/flnc-seqkit-stats/seqkit_flnc_wide.stats.tsv"
    output:
        tsv = "06_qc-reports/isoform-filtering/{prefix}_{suffix}_{filtertag}/isoform-filter-stats.tsv"
    log:
        "logs/06_qc-reports/isoform-filtering/{prefix}_{suffix}_{filtertag}_isoform-filter-stats.log"
    benchmark:
        "benchmarks/06_qc-reports/isoform-filtering/{prefix}_{suffix}_{filtertag}_isoform-filter-stats.txt"
    threads: 1
    params:
        filtertag = FILTERTAG
    conda:
        SNAKEDIR + "envs/qc-env.yaml"
    run:
        import os, re

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        def uniq_count(path):
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                return 0
            with open(path, "r") as fh:
                ids = {ln.strip() for ln in fh if ln.strip()}
            return len(ids)

        rows = []
        active_fail_files = list(input.fail_ids)
        filt_name_re = re.compile(r"/(filt_[^/]+)/")
        for f in active_fail_files:
            if not os.path.exists(f):
                continue
            m = filt_name_re.search(f)
            filt = m.group(1) if m else "unknown_filter"
            n = uniq_count(f)
            rows.append(("filter_counts", "-", "removed", str(n), filt))

        original_total = uniq_count(input.orig_ids)
        remaining_pass = uniq_count(input.pass_ids)
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

        assigned_by_sample = {}
        with open(input.pass_exp) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            samples = header[1:]
            sums = [0]*len(samples)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                for i in range(1, min(len(parts), len(samples)+1)):
                    try:
                        sums[i-1] += float(parts[i])
                    except ValueError:
                        pass
            for i, s in enumerate(samples):
                assigned_by_sample[s] = sums[i]

        raw_by_sample = {}
        with open(input.seqkit_cohort) as fh:
            header = fh.readline().rstrip("\n").split("\t")
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

        for s in sorted(assigned_by_sample.keys()):
            raw = raw_by_sample.get(s, 0.0)
            asg = assigned_by_sample.get(s, 0.0)
            pct = (asg / raw * 100.0) if raw > 0 else 0.0
            rows.append(("per_sample", s, "raw_reads",      f"{int(raw)}", "-"))
            rows.append(("per_sample", s, "assigned_reads", f"{int(asg)}", "-"))
            rows.append(("per_sample", s, "retained_pct",   f"{pct:.4f}", "-"))

        with open(output.tsv, "w") as out:
            out.write("table\tsample\tmetric\tvalue\tfilter\n")
            for r in rows:
                out.write("\t".join(r) + "\n")



# ────────────────────────────────────────────────
# Pivot any ..._long.tsv → ..._wide.tsv  (samples = rows, metrics = columns)
# ────────────────────────────────────────────────

rule pivot_long_to_wide:
    message: "Pivot to wide: {wildcards.base}_long.tsv → {wildcards.base}_wide.tsv"
    input:
        long = "{base}_long.tsv"
    output:
        wide = "{base}_wide.tsv"
    log:
        "logs/pivot/{base}.wide.log"
    benchmark:
        "benchmarks/pivot/{base}.wide.txt"
    threads: 1
    conda:
        SNAKEDIR + "envs/qc-env.yaml"   # must include python + pandas
    run:
        import os
        import pandas as pd

        # Read long table
        df = pd.read_csv(input.long, sep="\t")

        # Validate required columns
        required = {"sample", "metric", "value"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns in {input.long}: {missing}")

        # Preserve first-seen order
        sample_order = pd.unique(df["sample"])
        metric_order = pd.unique(df["metric"])

        # Coerce values to numeric where possible
        df["value"] = pd.to_numeric(df["value"], errors="ignore")

        # Aggregation strategy for duplicate (sample, metric) pairs
        agg_name = str(config.get("pivot_long_to_wide_agg", "first")).lower()
        aggfunc = {"first": "first", "sum": "sum", "mean": "mean"}.get(agg_name, "first")

        # Pivot
        wide = df.pivot_table(index="sample", columns="metric", values="value", aggfunc=aggfunc)
        wide = wide.reindex(index=sample_order, columns=metric_order)

        # Write
        wide.index.name = "sample"
        wide.reset_index(inplace=True)
        os.makedirs(os.path.dirname(output.wide), exist_ok=True)
        wide.to_csv(output.wide, sep="\t", index=False, na_rep="NA")

