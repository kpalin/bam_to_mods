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
LD_PRELOAD=/scratch/kpalin/miniconda3/envs/hts_dev/lib/libmimalloc.so.2.0 MIMALLOC_SHOW_ERRORS=1 MIMALLOC_VERBOSE=1 \
    ./bam_to_mods -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta -i /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.cram
