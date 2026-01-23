#!/bin/sh

# 19.01.2026 12:22:52 EST

#----------------------------------------------------------------------
# Utility functions
#----------------------------------------------------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

die() {
    echo "ERROR: $*" >&2
    [[ $- == *i* ]] && return 1
    exit 1
}


#----------------------------------------------------------------------
# Specify minimum input files
#----------------------------------------------------------------------

usage() {
    echo -e "\nUsage: $0 -a <ANNOTATE_GTF> -b <COLLAPSE_BASE>\n"
    echo -e "  -a: Path to the annotated GTF file from isoPropeller-annotate"
    echo -e "      (e.g. 06_tracks/ISOP_depth-gt1_isoqc_pass_defrag_patched_extra_stopfix.gtf)\n"
    echo -e "  -b: Base path + prefix of the isoPropeller-collapse source files"
    echo -e "      (e.g. 07_isoPropeller-defrag/ISOP_depth-gt1_tpm0.5s3/ISOP_depth-gt1_isoqc_pass_defrag)\n"
    exit 1
}

# Parse options
while getopts "a:b:" opt; do
    case $opt in
        a) ISOP_ANNOTATE_GTF="$OPTARG" ;;
        b) ISOP_COLLAPSE_BASE="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$ISOP_ANNOTATE_GTF" || -z "$ISOP_COLLAPSE_BASE" ]]; then
    log "Error: Missing required arguments."
    usage
fi


#----------------------------------------------------------------------
# Derive other required inputs (or specify)
#----------------------------------------------------------------------

# Derive the other required paths from the input files
OUTPUT_FOLDER="$(dirname ${ISOP_COLLAPSE_BASE}})_novelgene"
FILT_PREFIX="$(basename $(dirname ${ISOP_COLLAPSE_BASE}}))"
FILT_FOLDER="$(dirname ${ISOP_COLLAPSE_BASE}} | sed 's/\(\/07_\|\/08_\).*//' )/05_isoPropeller-filter/${FILT_PREFIX}"
UNFILT_PREFIX="$(dirname ${ISOP_COLLAPSE_BASE}} | sed 's/\(\/07_\|\/08_\).*//' )/04_isoPropeller-merge/ISOP_depth-gt1"

# Start logging
log "Validating inputs"

# Validate inputs
for VALIDATE in "${ISOP_ANNOTATE_GTF}" \
                "${ISOP_COLLAPSE_BASE}_exp.txt" \
                "${ISOP_COLLAPSE_BASE}.gtf" \
                "${ISOP_COLLAPSE_BASE}_id.txt" \
                "${ISOP_COLLAPSE_BASE}_modal_ends.gtf" \
                "${ISOP_COLLAPSE_BASE}.trackgroups" \
                "${ISOP_COLLAPSE_BASE}_tss.bed" \
                "${ISOP_COLLAPSE_BASE}_tts.bed" \
                "${FILT_FOLDER}" \
                "${UNFILT_PREFIX}.bed" \
                "${UNFILT_PREFIX}_exp.txt" \
                "${UNFILT_PREFIX}.gff" \
                "${UNFILT_PREFIX}.gtf" \
                "${UNFILT_PREFIX}_id.txt" \
                "${UNFILT_PREFIX}_modal_ends.gtf" \
                "${UNFILT_PREFIX}_tss.bed" \
                "${UNFILT_PREFIX}_tss_count.txt" \
                "${UNFILT_PREFIX}_tts.bed" \
                "${UNFILT_PREFIX}_tts_count.txt"
do
    if [[ ! -e "$VALIDATE" ]]; then
        die "${VALIDATE} not found"
    fi
done

log "Completed input validation"

# If we passed validations, create the output folder
mkdir -p "${OUTPUT_FOLDER}/tmp"

#----------------------------------------------------------------------
# Get isoforms corresponding to novel genes from the annotate pipeline
#----------------------------------------------------------------------

log "Getting novelgene isoform IDs from the annotate pipeline outputs"

# Map isoforms to isoPropeller collapse status
gtf-count-attributes.pl -a transcript_id,status "${ISOP_ANNOTATE_GTF}" > "${OUTPUT_FOLDER}/tmp/mapping-tid-status.txt"

# Extract all isoforms with antisense, divergent, intergenic, intronic, ncRNA_host_gene or overlapping status
awk 'BEGIN{FS=OFS="\t"}  $2=="antisense" || $2=="divergent" || $2=="intergenic" || $2=="intronic" || $2=="ncRNA_host_gene" || $2=="overlapping" {print $1}' \
    "${OUTPUT_FOLDER}/tmp/mapping-tid-status.txt" \
    > "${OUTPUT_FOLDER}/tmp/novelgene-tid.txt"


#----------------------------------------------------------------------
# From the unfiltered file, get other isoforms overlapping novel genes
#----------------------------------------------------------------------

log "Gather all other isoforms overlapping novelgenes by splice junction or exon boundary overlap"

# From the query file we grab all exon coordinates including strand for overlap matching
gtf-filter-attributes.pl -a transcript_id -m "${OUTPUT_FOLDER}/tmp/novelgene-tid.txt" "${ISOP_ANNOTATE_GTF}" \
    | awk 'BEGIN{FS=OFS="\t"} $3=="exon"' \
    > "${OUTPUT_FOLDER}/tmp/novelgene-tid.gtf"

gtf2gff.pl "${OUTPUT_FOLDER}/tmp/novelgene-tid.gtf" \
    | awk 'BEGIN{FS=OFS="\t"} {print $1":"$4"-"$5":"$7, $9}' \
    | sort -u \
    > "${OUTPUT_FOLDER}/tmp/novelgene-coordinates.txt"

# Now we do the same for the unfiltered file
awk 'BEGIN{FS=OFS="\t"} {print $1":"$4"-"$5":"$7, $9}' "${UNFILT_PREFIX}.gff" \
    | sort -u \
    > "${OUTPUT_FOLDER}/tmp/unfilt-coordinates.txt"

# Find all isoform IDs in the unfiltered file that overlap the query file
intersect-by-ids -ff "${OUTPUT_FOLDER}/tmp/unfilt-coordinates.txt" -if "${OUTPUT_FOLDER}/tmp/novelgene-coordinates.txt" \
    | cut -f 2 \
    | sort -u \
    > "${OUTPUT_FOLDER}/tmp/new-isoforms_all-tid.txt"

# Now remove from this list all isoforms that are already in the reference
ids-uniq2left \
    "${OUTPUT_FOLDER}/tmp/new-isoforms_all-tid.txt" \
    "${OUTPUT_FOLDER}/tmp/novelgene-tid.txt" \
    > "${OUTPUT_FOLDER}/tmp/new-isoforms_extra-tid.txt"


#----------------------------------------------------------------------
# Re-apply all filters except the TPM filter
#----------------------------------------------------------------------

log "Re-apply all isoform filters except the TPM filter"

# Aggregate all filter outputs
rm -f "${OUTPUT_FOLDER}/tmp/isoforms-to-filter_tid.tmp"
for FILT in "${FILT_FOLDER}/"filt_*/
do
    if [[ "$(basename ${FILT})" != "filt_min_tpm" ]]
    then
        for FILE in "${FILT}/"*.ids
        do
            cat "${FILE}" >> "${OUTPUT_FOLDER}/tmp/isoforms-to-filter_tid.tmp"
        done
    else
        echo "Skipping filt_min_tpm"
    fi
done

# Compact to unique IDs only
sort -u "${OUTPUT_FOLDER}/tmp/isoforms-to-filter_tid.tmp" > "${OUTPUT_FOLDER}/tmp/isoforms-to-filter_tid.txt"
rm -f "${OUTPUT_FOLDER}/tmp/isoforms-to-filter_tid.tmp"

# Now subtract those IDs from the extra novel isoform files
ids-uniq2left \
    "${OUTPUT_FOLDER}/tmp/new-isoforms_extra-tid.txt" \
    "${OUTPUT_FOLDER}/tmp/isoforms-to-filter_tid.txt" \
    > "${OUTPUT_FOLDER}/tmp/new-isoforms_extra-tokeep-tid.txt"

#----------------------------------------------------------------------
# Combine the novel isoform IDs with the existing isoform IDs
#----------------------------------------------------------------------

log "Adding extra isoforms for novel genes to the existing set"

# Grab all the isoform IDs from the collapse folder
cut -f1 "${ISOP_COLLAPSE_BASE}_exp.txt" \
    | grep -v '^#)' \
    > "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.tmp"

# Combine them with the extra novel gene isoforms
cat \
    "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.tmp" \
    "${OUTPUT_FOLDER}/tmp/new-isoforms_extra-tokeep-tid.txt" \
    | sort -u \
    > "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt"

# Make sure we have transcript_id in there for matching
echo "transcript_id" >> "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt"

# Cleanup
rm -f "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.tmp" 

#----------------------------------------------------------------------
# And finally, generate a new set of supplemented output files
#----------------------------------------------------------------------

log "Generating isopropeller-collapse outputs, augmented with extra novelgene isoforms"

OUT_PREFIX="$(basename ${ISOP_COLLAPSE_BASE})"

# Generate the _exp.txt file
intersect-by-ids \
    -ff "${UNFILT_PREFIX}_exp.txt" \
    -if "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt" \
    | perl -pe 's/^transcript_id/#TranscriptID/' \
    > "${OUTPUT_FOLDER}/${OUT_PREFIX}_exp.tmp"

# Reprocess the header of the expression matrix to be compatible with the isoPropeller-annotate pipeline
awk 'BEGIN{OFS="\t"} NR==1 { $1="#TranscriptID"; for(i=2; i<=NF; i++) { split($i, parts, "/"); $i=parts[2] } } 1' \
    "${OUTPUT_FOLDER}/${OUT_PREFIX}_exp.tmp" > "${OUTPUT_FOLDER}/${OUT_PREFIX}_exp.txt"
rm -f "${OUTPUT_FOLDER}/${OUT_PREFIX}_exp.tmp"

# Generate the max ends GTF file
gtf-filter-attributes.pl \
    -a transcript_id \
    -m "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt"  \
    "${UNFILT_PREFIX}.gtf" \
    > "${OUTPUT_FOLDER}/${OUT_PREFIX}.gtf"

# Generate the modal ends GTF file
gtf-filter-attributes.pl \
    -a transcript_id \
    -m "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt"  \
    "${UNFILT_PREFIX}_modal_ends.gtf" \
    > "${OUTPUT_FOLDER}/${OUT_PREFIX}_modal_ends.gtf"

# Generate the ID file
intersect-by-ids \
    -ff "${UNFILT_PREFIX}_id.txt" \
    -fc 4 \
    -if "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt" \
    | perl -pe 's/^transcript_id/#TranscriptID/' \
    > "${OUTPUT_FOLDER}/${OUT_PREFIX}_id.txt"
    
# Generate the TSS file
intersect-by-ids \
    -ff "${UNFILT_PREFIX}_tss.bed" \
    -fc 4 \
    -if "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt" \
    > "${OUTPUT_FOLDER}/${OUT_PREFIX}_tss.bed"
    
# Generate the TTS file
intersect-by-ids \
    -ff "${UNFILT_PREFIX}_tts.bed" \
    -fc 4 \
    -if "${OUTPUT_FOLDER}/tmp/FINAL-ISOFORMS.txt" \
    > "${OUTPUT_FOLDER}/${OUT_PREFIX}_tts.bed"

# Copy the trackgroups file
cp "${ISOP_COLLAPSE_BASE}.trackgroups" "${OUTPUT_FOLDER}/${OUT_PREFIX}.trackgroups"


