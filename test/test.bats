setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    # ... the remaining setup is unchanged

    # get the containing directory of this file
    # use $BATS_TEST_FILENAME instead of ${BASH_SOURCE[0]} or $0,
    # as those will point to the bats executable's location or the preprocessed file respectively
    DIR="$(cd "$(dirname "$BATS_TEST_FILENAME")" >/dev/null 2>&1 && pwd)"
    # make executables in src/ visible to PATH
    #PATH="$DIR/../src:$PATH"
}

@test "T-a and A+a match" {

    ./debug/bam_to_mods -m A+a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${DIR}/tmp/A+a.0.out
    ./debug/bam_to_mods -m T-a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${DIR}/tmp/T-a.0.out
    diff <(cut -f-8,10 ${DIR}/tmp/A+a.0.out) <(cut -f-8,10 ${DIR}/tmp/T-a.0.out)
    RET=$(duckdb <${DIR}/check_read_depths.sql)
    duckdb <${DIR}/check_read_depths.sql >${DIR}/tmp/T_depth.out
    diff ${DIR}/check_read_depths.expected ${DIR}/tmp/T_depth.out
}
