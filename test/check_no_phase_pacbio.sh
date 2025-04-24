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

./debug/bam_to_mods -m CG+m.0 -m A+a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/tmp_strand.out
./debug/bam_to_mods --no_phase -m CG+m.0 -m A+a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/tmp.out

TXT=$(
    duckdb <<EOF
create table nosplit as
select
    *
from
    read_csv('${TEMPDIR}/tmp.out');

create table withsplit as
select
    *
from
    read_csv('${TEMPDIR}/tmp_strand.out');

.mode csv
.header off

select
    *
from
    nosplit a
    LEFT JOIN withsplit b ON row(a."start", a."haplotype") = row(b."start", b."haplotype")
where
    a.called_reads != b.called_reads
    or a.modified_reads != b.modified_reads
    or a.uncalled_reads != b.uncalled_reads
    or a.mismatch_reads != b.mismatch_reads;

EOF
)
test -z "$TXT" || {
    echo "$TXT"
    exit 1
}
#echo BAD LINES ${BADLINES:-x}
#[[ ${BADLINES:-1} == 0 ]]

#echo DONE
