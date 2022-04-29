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

bedtools slop -b 3000 -g /home/kpalin/src/myscripts/nanopore/data/GRCh38_no_alt.fasta.fai \
    -i /home/kpalin/src/myscripts/nanopore/data/CBSs_39bp_clean.GRCh38.bed |
    bedtools merge -i - >"${TEMPDIR}/regions.bed"

samtools view -F 1796 -q 20 -L "${TEMPDIR}/regions.bed" -o "subset.cram" /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.sorted.cram
samtools index subset.cram
