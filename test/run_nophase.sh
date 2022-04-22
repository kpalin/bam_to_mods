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
./bam_to_mods --mod CG+m.0 --mod CG+h.0 --no_phase -r /mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta -i /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.cram
