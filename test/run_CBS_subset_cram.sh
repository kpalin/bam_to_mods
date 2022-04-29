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

LIKE_LIM=${1:-127}

/usr/bin/time -v ./bam_to_mods --no_phase -c "${LIKE_LIM}" \
    -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta \
    --mod CG+m.0 \
    -i subset.cram | bgzip >"subset.mC.${LIKE_LIM}.tsv.gz"

zcat "subset.mC.${LIKE_LIM}.tsv.gz" | awk '0.2*$6>$8' | cut -f-3,11 >"subset.mC.${LIKE_LIM}.bedgraph"
python /home/kpalin/src/myscripts/nanopore/methylation_CBS_wave_plot.py \
    -i "subset.mC.${LIKE_LIM}.bedgraph" \
    -t 1 -o methylation_frequency.Fam_c461_1_19_0711NK_meg.subset.CpG_CTCF_wave.c${LIKE_LIM}.png
