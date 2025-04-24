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

@test "can launch the program" {
    ./debug/bam_to_mods -h
}

@test "can run region on ONT" {
    ./debug/bam_to_mods -i data/fibreseq_demo_ont.bam -r data/ref.fa.gz -R chr20:46138245-46150899
}

@test "can run unphased region on PacBio" {
    ./debug/bam_to_mods --no_phase -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899 | awk '$0!~/^#/ && $4!~/^[N*]$/ {print "Line",NR ":",$0 > "/dev/stderr" ; exit 1;}'
}
# @test "can run CG+m.0  pacbio" {
#     ./debug/bam_to_mods -m CG+m.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899
#     # Covered 963 sites.
# }

@test "can run A+a.0 split_strand fibreseq pacbio" {
    ./debug/bam_to_mods --split_strand -m A+a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899
}

@test "can run T-a.0 split_strand fibreseq pacbio" {
    ./debug/bam_to_mods --split_strand -m T-a.0 -i data/fibreseq_demo_pacbio.bam -r data/ref.fa.gz -R chr20:46138245-46150899
}
