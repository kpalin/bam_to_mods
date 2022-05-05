#!/bin/bash

#  Created on Friday, 22 April  2022

# For debugging
set -o verbose

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail

LOD=${1:-2}
MAX_C=$(echo "256*e(-$LOD)/(e(-$LOD)+1)" | bc -l | cut -d. -f1)
MIN_C=$(echo "256*e($LOD)/(e($LOD)+1)" | bc -l | cut -d. -f1)

LD_PRELOAD=/scratch/kpalin/miniconda3/envs/hts_dev/lib/libmimalloc.so.2.0 MIMALLOC_SHOW_ERRORS=1 \
    ./bam_to_mods -C $MAX_C -c $MIN_C -R chr9:170000-180000 \
    -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta \
    --mod CG+m.0 --mod CG+h.0 \
    -i /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.cram
