#!/bin/bash

#  Created on Friday, 22 April  2022

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

WRKDIR=$(dirname "$(readlink -f $0)")
bash $WRKDIR/run_region.sh >$TEMPDIR/run_region_cur.out

diff $TEMPDIR/run_region_cur.out $WRKDIR/run_region.out
