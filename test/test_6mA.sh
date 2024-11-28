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

./debug/bam_to_mods -r /mnt/cg8/reference-genomes/chm13v2.0_maskedY_rCRS/chm13v2.0_maskedY_rCRS.fa -i test/test_6aA_5mC_5hC.bam -m A+a.0 -m CG+m.0 \
    >$TEMPDIR/run_Aa.out ||true

diff test/test_6aA_5mC_5hC.Aa.out $TEMPDIR/run_Aa.out
