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

./debug/bam_to_mods --split_strand -m CG+m.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/tmp_strand.out
./debug/bam_to_mods -m CG+m.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${TEMPDIR}/tmp.out

BADLINES=$(
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

create view sumsplit as
select
    (position - IF(strand = '-', 1, 0)) as "start",
    haplotype,
    sum("called_reads") as sum_called_reads,
    sum("uncalled_reads") as sum_uncalled_reads,
    sum("mismatch_reads") as sum_mismatched_reads,
    sum("modified_reads") as sum_modified_reads
from
    withsplit
group by
    "start",
    "haplotype";

.mode tabs
.header off
select
    count(*)
from
    (
        SELECT
            *
        FROM
            nosplit
        WHERE
            haplotype != '*'
    ) AS nosplit FULL
    OUTER JOIN sumsplit ON nosplit."start" = sumsplit."start"
    AND sumsplit."haplotype" = nosplit.haplotype
WHERE
    sum_called_reads != called_reads
    OR sum_uncalled_reads != uncalled_reads
    OR mismatch_reads != sum_mismatched_reads
    OR modified_reads != sum_modified_reads;
EOF
)
echo BAD LINES ${BADLINES:-x}
[[ ${BADLINES:-1} == 0 ]]

#echo DONE
