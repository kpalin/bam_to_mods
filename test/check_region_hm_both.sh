#!/bin/bash

#  Created on Tuesday, 05 April  2022

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

./bam_to_mods --split_strand -R chr9:170000-180000 -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta \
    --mod CG+m.0 \
    -i /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.cram \
    >$TEMPDIR/run_region_cur_m.out

awk '$9!="h"' test/run_region.out >$TEMPDIR/run_region_orig_m.out

echo "Methylation check:"

diff $TEMPDIR/run_region_cur_m.out $TEMPDIR/run_region_orig_m.out

./bam_to_mods --split_strand -R chr9:170000-180000 -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta \
    --mod CG+h.0 \
    -i /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.cram \
    >$TEMPDIR/run_region_cur_h.out

awk '$9!="m"' test/run_region.out >$TEMPDIR/run_region_orig_h.out

echo "Hydroxymethylation check:"

diff $TEMPDIR/run_region_cur_h.out $TEMPDIR/run_region_orig_h.out
