import os
from snakemake.io import glob_wildcards

# ───────────────────────────────────────────────
# Rule 1: Convert BAM, split into chunks (Scatter)
# This is a checkpoint because the number of output
# chunks is not known before this rule is run.
# ───────────────────────────────────────────────
checkpoint bam_to_sam_and_split:
    message: "Convert BAM to SAM and split into chunks for {wildcards.sample}"
    input:
        bam="01_mapping/{sample}/{sample}_mapped_labeled.bam"
    output:
        # The output is a directory that will contain the subdirectories for each chunk.
        # This directory also acts as a flag for completion for the checkpoint.
        directory("02_transcriptclean/{sample}/split/")
    log:
        "logs/02_transcriptclean/{sample}_bam_to_sam_split.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_bam_to_sam_split.txt"
    params:
        reads_per_chunk = TRANSCRIPTCLEAN_CHUNK_MAXREADS,
        # Define the path for the manifest file listing all created chunks.
        chunk_list_file=lambda wildcards, output: f"{output}/chunk_list.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        r"""
        (
        echo "Splitting {input.bam} for {wildcards.sample} into chunks of {params.reads_per_chunk} reads."

        # Define the main output directory for the final SAM chunks
        SPLIT_DIR="{output}"
        mkdir -p "$SPLIT_DIR"

        # MODIFIED: Create a temporary directory inside the workflow structure.
        ALIGN_DIR=$(mktemp -d -p "$SPLIT_DIR")
        # Define the location for the temporary header file
        HEADER_FILE="$ALIGN_DIR/header.sam"

        # 1. Save the SAM header to the temp file
        echo "Extracting SAM header."
        samtools view -H "{input.bam}" > "$HEADER_FILE"

        # 2. Split only the alignment records into temporary chunks
        echo "Splitting alignment records."
        samtools view -@ {threads} "{input.bam}" | \
            split -l {params.reads_per_chunk} --numeric-suffixes=1 --additional-suffix=.body - "$ALIGN_DIR/align.chunk_"

        # 3. Rebuild full SAM chunks by prepending the header to each alignment chunk
        echo "Rebuilding full SAM chunks with headers into dedicated subdirectories."
        for align_chunk in "$ALIGN_DIR"/align.chunk_*; do
            # Extract the numeric suffix (e.g., 01) from the temp file name (e.g., align.chunk_01.body)
            chunk_suffix=$(basename "$align_chunk" | cut -d_ -f2 | cut -d. -f1)

            # Create a dedicated subdirectory for each chunk
            CHUNK_SUBDIR="$SPLIT_DIR/$chunk_suffix"
            mkdir -p "$CHUNK_SUBDIR"

            # Construct the final output file name inside the new subdirectory
            final_chunk_file="$CHUNK_SUBDIR/{wildcards.sample}.chunk_$chunk_suffix.sam"

            # Concatenate the header and the alignment chunk into the final output file
            cat "$HEADER_FILE" "$align_chunk" > "$final_chunk_file"
        done

        # 4. Create a manifest of the chunks for downstream rules.
        # MODIFIED: This command now specifically finds only numeric directory names,
        # ignoring the temporary directory created by mktemp.
        echo "Creating chunk manifest file."
        find "$SPLIT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {{}} \; | grep '^[0-9][0-9]*$' | sort -n > "{params.chunk_list_file}"

        # 5. Clean up the temporary directory containing the header and alignment chunks
        echo "Cleaning up temporary split files."
        rm -rf "$ALIGN_DIR"

        echo "Finished splitting SAM for {wildcards.sample}"
        ) &> "{log}"
        """

# ───────────────────────────────────────────────
# Rule 2: Run transcriptclean on each chunk
# ───────────────────────────────────────────────
rule run_transcriptclean_on_chunk:
    message: "Running transcriptclean on {wildcards.sample}, chunk {wildcards.chunk}"
    input:
        # The input SAM path reflects the chunk subdirectory structure
        sam = "02_transcriptclean/{sample}/split/{chunk}/{sample}.chunk_{chunk}.sam",
        ref = GENOMEFASTA,
        fai = GENOMEFASTA + ".fai"
    output:
        # Outputs are kept in the unique subdirectory for each chunk to prevent conflicts.
        cleaned_bam     = temp("02_transcriptclean/{sample}/split/{chunk}/{sample}.chunk_{chunk}_clean.bam"),
        clean_log_gz    = temp("02_transcriptclean/{sample}/split_logs/{sample}.chunk_{chunk}_clean.log.gz"),
        clean_te_log_gz = temp("02_transcriptclean/{sample}/split_logs/{sample}.chunk_{chunk}_clean.TE.log.gz")
    log:
        "logs/02_transcriptclean/{sample}_transcriptclean_chunk_{chunk}.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_transcriptclean_chunk_{chunk}.txt"
    params:
        # The output prefix points to the unique chunk subdirectory.
        outprefix      = lambda wildcards: f"02_transcriptclean/{wildcards.sample}/split/{wildcards.chunk}/{wildcards.sample}.chunk_{wildcards.chunk}",
        junctions_file = SPLICEJUNCTIONS
    threads: TRANCRIPTCLEAN_CHUNK_THREADS
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        r"""
        (
        # The CHUNK_DIR is already created by the checkpoint rule, but mkdir -p is safe.
        CHUNK_DIR="02_transcriptclean/{wildcards.sample}/split/{wildcards.chunk}"
        mkdir -p "$CHUNK_DIR"

        echo "Starting transcriptclean for {input.sam}"

        # Define temporary paths for logs and the intermediate clean SAM, now inside the chunk dir
        TMP_CLEAN_SAM="{params.outprefix}_clean.sam"
        TMP_CLEAN_LOG="{params.outprefix}_clean.log"
        TMP_CLEAN_TE_LOG="{params.outprefix}_clean.TE.log"

        transcriptclean \
            -s "{input.sam}" \
            -g "{input.ref}" \
            -j "{params.junctions_file}" \
            -t {threads} \
            --primaryOnly \
            --canonOnly \
            --deleteTmp \
            -o "{params.outprefix}"

        # Convert the cleaned SAM to BAM for proper merging later.
        echo "Converting cleaned SAM to BAM for chunk {wildcards.chunk}"
        samtools view -@ {threads} -bS "$TMP_CLEAN_SAM" > "{output.cleaned_bam}"
        rm "$TMP_CLEAN_SAM" # Remove the intermediate SAM file

        # Ensure the directory for compressed logs exists
        mkdir -p "02_transcriptclean/{wildcards.sample}/split_logs/"

        # Compress the log files to their final destination
        pigz -f -c -p {threads} "$TMP_CLEAN_LOG" > "{output.clean_log_gz}"
        pigz -f -c -p {threads} "$TMP_CLEAN_TE_LOG" > "{output.clean_te_log_gz}"

        # Clean up the original, uncompressed log files
        rm "$TMP_CLEAN_LOG" "$TMP_CLEAN_TE_LOG"

        # Clean up the input SAM chunk as it is no longer needed.
        echo "Cleaning up input SAM chunk: {input.sam}"
        rm "{input.sam}"

        echo "Finished transcriptclean for chunk {wildcards.chunk}"
        ) &> "{log}"
        """

# ───────────────────────────────────────────────
# Helper functions to gather all chunked files 
# after the checkpoint is complete.
# ───────────────────────────────────────────────

def get_chunk_wildcards(wildcards):
    """Helper to get chunk numbers for a sample by reading the manifest file."""
    checkpoint_dir = checkpoints.bam_to_sam_and_split.get(sample=wildcards.sample).output[0]
    # Read the chunk numbers from the manifest file instead of globbing.
    chunk_list_file = os.path.join(checkpoint_dir, "chunk_list.txt")
    with open(chunk_list_file) as f:
        chunk_numbers = [line.strip() for line in f if line.strip()]
    return chunk_numbers

def get_cleaned_chunks(wildcards):
    """Gathers all cleaned BAM chunks for a sample."""
    # Path now includes the {chunk} subdirectory wildcard.
    return expand("02_transcriptclean/{sample}/split/{chunk}/{sample}.chunk_{chunk}_clean.bam",
                  sample=wildcards.sample,
                  chunk=get_chunk_wildcards(wildcards))

def get_clean_logs(wildcards):
    """Gathers all standard log chunks for a sample."""
    return expand("02_transcriptclean/{sample}/split_logs/{sample}.chunk_{chunk}_clean.log.gz",
                  sample=wildcards.sample,
                  chunk=get_chunk_wildcards(wildcards))

def get_clean_te_logs(wildcards):
    """Gathers all TE log chunks for a sample."""
    return expand("02_transcriptclean/{sample}/split_logs/{sample}.chunk_{chunk}_clean.TE.log.gz",
                  sample=wildcards.sample,
                  chunk=get_chunk_wildcards(wildcards))


# ───────────────────────────────────────────────
# Rule 3: Gather results, merge, and index (Gather)
# ───────────────────────────────────────────────
rule gather_results:
    message: "Gathering all results for {wildcards.sample}: merging, sorting, and indexing BAMs and logs."
    input:
        bams       = get_cleaned_chunks,
        clean_logs = get_clean_logs,
        te_logs    = get_clean_te_logs
    output:
        merged_bam    = temp("02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam"),
        merged_bai    = temp("02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam.bai"),
        merged_log    = "02_transcriptclean/{sample}/{sample}_final_clean.log.gz",
        merged_te_log = "02_transcriptclean/{sample}/{sample}_final_clean.TE.log.gz"
    log:
        "logs/02_transcriptclean/{sample}_gather_results.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_gather_results.txt"
    params:
        bams_count       = lambda wildcards, input: len(input.bams),
        clean_logs_count = lambda wildcards, input: len(input.clean_logs),
        te_logs_count    = lambda wildcards, input: len(input.te_logs),
        unsorted_bam     = lambda wildcards, output: f"{output.merged_bam}.unsorted.bam"
    threads: 4
    conda:
        SNAKEDIR + "envs/transcriptclean.yaml"
    shell:
        r"""
        (
        echo "Merging {params.bams_count} cleaned BAM chunks for {wildcards.sample}"
        # Merge to a temporary file first to make the process more robust.
        samtools merge -f -@ {threads} "{params.unsorted_bam}" {input.bams}
        echo "Finished merging BAMs."

        echo "Sorting merged BAM file."
        # Now sort the temporary merged file into the final output file.
        samtools sort -@ {threads} -o "{output.merged_bam}" "{params.unsorted_bam}"
        echo "Finished sorting BAM."

        # Clean up the unsorted temporary file
        rm "{params.unsorted_bam}"

        echo "Indexing temporary merged BAM file."
        samtools index -@ {threads} "{output.merged_bam}"
        echo "Finished indexing."

        echo "Consolidating {params.clean_logs_count} standard log files."
        zcat -f {input.clean_logs} | pigz -c -p {threads} > "{output.merged_log}"
        echo "Finished consolidating standard logs."

        echo "Consolidating {params.te_logs_count} TE log files."
        zcat -f {input.te_logs} | pigz -c -p {threads} > "{output.merged_te_log}"
        echo "Finished consolidating TE logs."

        # Clean up the now-empty intermediate chunk directory to save disk space.
        echo "Cleaning up intermediate chunk directory."
        rm -rf "02_transcriptclean/{wildcards.sample}/split/"
        echo "Cleanup complete."

        echo "All intermediate results gathered for {wildcards.sample}"
        ) &> "{log}"
        """

# ───────────────────────────────────────────────
# Rule 4: Downsample reads on chrM to avoid slowdowns
# ───────────────────────────────────────────────
rule downsample_chrM:
    message: "Downsample mitochondrial reads for {wildcards.sample}"
    input:
        bam = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam",
        bai = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean_temp.bam.bai"
    output:
        bam = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam",
        bai = "02_transcriptclean/{sample}/{sample}_mapped_labeled_tclean.bam.bai"
    params:
        max_chrM_reads = MAXCHRMREADS,
        seed           = 4232
    log:
        "logs/02_transcriptclean/{sample}_downsample_chrM.log"
    benchmark:
        "benchmarks/02_transcriptclean/{sample}_downsample_chrM.txt"
    threads: 4
    conda:
        SNAKEDIR + "envs/pysam.yaml"
    script:
        SNAKEDIR + "scripts/downsample_chrM_bam.py"
