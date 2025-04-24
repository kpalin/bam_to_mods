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

@test "can run default pacbio" {
    ./debug/bam_to_mods -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz
}

@test "can run CG+m.0 GC+m.1 split strand pacbio" {
    run ./test/check_CG_GC_count_match.sh
    assert_success
}

@test "can run CG+m.0 split strand pacbio" {
    run ./test/check_strand_split_scores.sh
    assert_success
}
@test "Summary unphased on pacbio working" {
    run ./test/check_no_phase_pacbio.sh
    assert_success
}

@test "can consistently run CG+m.0, GC+m.1, CG+h.0 and A+a.0 on ont" {
    run ./test/check_multi_context_match.sh
    assert_success
}
