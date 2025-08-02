import pysam
import random
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

# Snakemake variables (assuming they are defined in the rule)
bam_in_path = Path(snakemake.input.bam)
bam_out_path = Path(snakemake.output.bam)
max_reads = snakemake.params.max_chrM_reads
seed = snakemake.params.seed
threads = snakemake.threads
logfile = Path(snakemake.log[0])

# Set a consistent seed for reproducibility
random.seed(seed)

with open(logfile, "w") as log:
    log.write(f"Starting downsample_chrM for {bam_in_path}\n")

    # Open input BAM and inspect reference names
    try:
        bam_in = pysam.AlignmentFile(bam_in_path, "rb")
    except Exception as e:
        log.write(f"Error opening input BAM: {e}\n")
        raise
    
    mito_candidates = ["chrM", "MT", "chrMT", "M"]
    contigs = set(bam_in.references)
    mito_name = next((name for name in mito_candidates if name in contigs), None)

    if mito_name:
        log.write(f"Mitochondrial contig detected: {mito_name}\n")
        
        # Use a temporary directory for robust cleanup
        with TemporaryDirectory() as tempdir:
            temp_dir_path = Path(tempdir)
            temp_unsorted_path = temp_dir_path / "temp_unsorted.bam"

            bam_unsorted = pysam.AlignmentFile(str(temp_unsorted_path), "wb", header=bam_in.header)
            
            # Initialize reservoir and counters for a single-pass approach
            chrM_reads = []
            mito_read_count_total = 0
            other_read_count = 0

            log.write("Performing reservoir sampling in a single pass.\n")

            for read in bam_in.fetch(until_eof=True):
                # Skip unmapped, secondary, or supplementary reads
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                
                if read.reference_name == mito_name:
                    mito_read_count_total += 1
                    if len(chrM_reads) < max_reads:
                        # Reservoir not yet full, add the read
                        chrM_reads.append(read)
                    else:
                        # Reservoir is full, randomly decide whether to replace an existing read
                        j = random.randrange(mito_read_count_total)
                        if j < max_reads:
                            chrM_reads[j] = read
                else:
                    # Non-mito read, write directly to the temporary file
                    bam_unsorted.write(read)
                    other_read_count += 1
            
            # Close the input BAM file
            bam_in.close()

            keep_n = min(max_reads, mito_read_count_total)

            # Log for clarity on whether downsampling was performed
            if mito_read_count_total <= max_reads:
                log.write(f"Number of mitochondrial reads ({mito_read_count_total}) is less than or equal to max_reads ({max_reads}). No downsampling performed.\n")
            else:
                log.write(f"Number of mitochondrial reads ({mito_read_count_total}) exceeds max_reads ({max_reads}). Downsampling to {keep_n} reads.\n")

            # Write the downsampled mito reads from the reservoir
            for read in chrM_reads:
                bam_unsorted.write(read)

            bam_unsorted.close()

            log.write(f"Wrote {other_read_count} non-mito reads.\n")
            log.write(f"Wrote {keep_n} downsampled mito reads.\n")

            # Sort and index the combined output
            log.write(f"Sorting and indexing the combined output BAM to {bam_out_path}.\n")
            pysam.sort("-@", str(threads), "-o", str(bam_out_path), str(temp_unsorted_path))
            pysam.index(str(bam_out_path))

            log.write(f"Output successfully written.\n")
        
        # The temporary directory and its contents are automatically removed here
        # (upon exiting the 'with' block).

    else:
        log.write("No mitochondrial contig found. Copying original BAM.\n")
        shutil.copy(bam_in_path, bam_out_path)
        pysam.index(str(bam_out_path))
