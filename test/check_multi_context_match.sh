#!/bin/bash

#  Created on Tuesday, 15 April  2025

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

BAMFILE="data/fibreseq_demo_ont.bam"

./debug/bam_to_mods -m CG+m.0 -i ${BAMFILE} -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/a.out
./debug/bam_to_mods --no_header -m GC+m.1 -i ${BAMFILE} -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/b.out
./debug/bam_to_mods --no_header -m CG+h.0 -i ${BAMFILE} -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/c.out
./debug/bam_to_mods --no_header -m A+a.0 -i ${BAMFILE} -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/d.out

./debug/bam_to_mods -m GC+m.1 -m CG+m.0 -m A+a.0 -m CG+h.0 -i ${BAMFILE} -r data/ref.fa.gz -R chr20:46138245-46150899 | sort >${TEMPDIR}/ab.out

sort ${TEMPDIR}/?.out >${TEMPDIR}/a_b_sorted.out

cmp ${TEMPDIR}/a_b_sorted.out ${TEMPDIR}/ab.out

[[ $? == 0 ]]
