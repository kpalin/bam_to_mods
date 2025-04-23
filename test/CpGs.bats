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

@test "Running ONT default region " {
    ./debug/bam_to_mods -i data/fibreseq_demo_ont.bam -r data/ref.fa.gz \
        -R chr20:46138245-46150899 >${DIR}/tmp/fibreseq_demo_ont.out

}

@test "ONT CG read depth about match" {
    duckdb <${DIR}/check_CG_read_depths.sql >${DIR}/tmp/CG_depth.out
    diff ${DIR}/check_CG_read_depths.expected ${DIR}/tmp/CG_depth.out
}

@test "ONT CG output match 0.4.2 version" {
    diff ${DIR}/fibreseq_demo_ont_master.expected ${DIR}/tmp/fibreseq_demo_ont.out
}
