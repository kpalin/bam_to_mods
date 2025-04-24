#!/bin/bash

#  Created on Monday, 14 April  2025

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

./debug/bam_to_mods --split_strand -m CG+m.0 -m GC+m.1 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/tmp_strand.out

BADLINES=$(
    duckdb <<EOF

create table withsplit as
select
    *
from
    read_csv('${TEMPDIR}/tmp_strand.out');

.mode tabs
.header off

select
    count(*)
from
    withsplit a,
    withsplit b
where
    a.position = b.position
    and a.strand = b.strand
    and a.haplotype = b.haplotype
    and a.call_code < b.call_code
    and (
        a.called_reads != b.called_reads
        or a.uncalled_reads != b.uncalled_reads
        or a.modified_reads != b.modified_reads
        or a.mismatch_reads != b.mismatch_reads
    ) ;

EOF
)
#echo BAD LINES ${BADLINES}
[[ ${BADLINES:-1} == 0 ]]

#echo DONE
