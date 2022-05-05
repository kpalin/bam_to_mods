#!/bin/bash

#  Created on Thursday, 28 April  2022

# For debugging
#set -o verbose

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail
# Make temporary directory which will be cleaned at script exit
TEMPDIR=$(mktemp --directory)
function _cleanup {
    rm -r "$TEMPDIR"
}
trap _cleanup EXIT

MIN_C=125
MAX_C=0

usage() {
    echo -e "usage:
$0 [-c $MIN_C] [-C $MAX_C] [-l LOD] 

-o OUTFILE
-h          Show this message and exit." >&2
    exit 1

}

while getopts "l:c:C:h" flag; do
    case "$flag" in
    l)
        LOD="$OPTARG"

        MAX_C=$(echo "256*e(-$LOD)/(e(-$LOD)+1)" | bc -l | cut -d. -f1)
        MIN_C=$(echo "256*e($LOD)/(e($LOD)+1)" | bc -l | cut -d. -f1)

        ;;
    c)
        MIN_C="$OPTARG"
        ;;
    C)
        MAX_C="$OPTARG"
        ;;
    h | *)
        usage
        ;;
    esac
done
shift $((OPTIND - 1))
OPTIND=1

LIKE_LIM="c${MIN_C}_C${MAX_C}"
echo "$LIKE_LIM"

/usr/bin/time -v ./release/bam_to_mods --no_phase -C $MAX_C -c $MIN_C \
    -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta \
    --mod CG+m.0 -@ 5 \
    -i subset.cram | bgzip >"threshold_tst/subset_tsv/subset.mC.${LIKE_LIM}.tsv.gz"

zcat "threshold_tst/subset_tsv/subset.mC.${LIKE_LIM}.tsv.gz" | awk '0.2*$6>$8' | cut -f-3,11 >"${TEMPDIR}/subset.mC.${LIKE_LIM}.bedgraph"
python /home/kpalin/src/myscripts/nanopore/methylation_CBS_wave_plot.py \
    -i "${TEMPDIR}/subset.mC.${LIKE_LIM}.bedgraph" \
    -t 1 -o threshold_tst/wave/methylation_frequency.Fam_c461_1_19_0711NK_meg.subset.CpG_CTCF_wave.${LIKE_LIM}.png

python /home/kpalin/src/myscripts/nanopore/methylation_CBS_wave_signal.py \
    -i threshold_tst/wave/methylation_frequency.Fam_c461_1_19_0711NK_meg.subset.CpG_CTCF_wave.${LIKE_LIM}.png.tsv \
    -o threshold_tst/signal/methylation_frequency.Fam_c461_1_19_0711NK_meg.subset.CpG_CTCF_wave.${LIKE_LIM}.signal.png
