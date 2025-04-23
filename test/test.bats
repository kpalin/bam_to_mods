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

@test "Running -m T-a " {
    ./debug/bam_to_mods -m T-a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${DIR}/tmp/T-a.0.out

}

@test "Running -m A+a " {
    ./debug/bam_to_mods -m A+a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 >${DIR}/tmp/A+a.0.out

}

@test "T-a and A+a output match" {
    diff <(cut -f-8,10 ${DIR}/tmp/A+a.0.out) <(cut -f-8,10 ${DIR}/tmp/T-a.0.out)
}

@test "T-a and A+a read depth about match" {
    duckdb <${DIR}/check_read_depths.sql >${DIR}/tmp/T_depth.out
    diff ${DIR}/check_read_depths.expected ${DIR}/tmp/T_depth.out

}

@test "Running -m A+a -m T-a and they match individually " {

    ./debug/bam_to_mods -m A+a.0 -m T-a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 | sort -u >${DIR}/tmp/AT+-a.0.out

    sort -u ${DIR}/tmp/T-a.0.out ${DIR}/tmp/A+a.0.out >${DIR}/tmp/merged.AT+-a.0.out
    diff ${DIR}/tmp/merged.AT+-a.0.out ${DIR}/tmp/AT+-a.0.out

}
