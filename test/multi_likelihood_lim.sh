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

for LIKE_LIM in $(seq 2 4 48); do
    LD_PRELOAD=/scratch/kpalin/miniconda3/envs/hts_dev/lib/libmimalloc.so.2.0 MIMALLOC_SHOW_ERRORS=1 \
        ./bam_to_mods --no_phase -c $LIKE_LIM -R chr20 -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta \
        --mod CG+m.0 \
        -i /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.cram | bgzip >longshot.Fam_c461_1_19_0711NK_meg.chr20.c${LIKE_LIM}.tsv.gz &
done

wait
