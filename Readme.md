# bam to mods

`bam_to_mods` outputs methylation levels from ML/MM tag formatted bam files to tab separated file,
similar to [modbam2bed](https://github.com/epi2me-labs/modbam2bed) from Oxford Nanopore. The main
dependence for `bam_to_mods` is fairly recent [htslib](https://github.com/samtools/htslib) (at least version 1.16).

```text
usage: debug/bam_to_mods [-c 209] [-C 46]  [-b 7]  [-q 20] [-E 1796] [-m CG+m.0] [-R chr1:1-100] [--no_header] [--no_phase] [--split_strand] -r ref.fasta -i input.cram

 -c|--min_mod_prob Minimum probability of modification called modified (range 0-255)
 -C|--max_mod_prob Maximum probability of modification called unmodified (range 0-255)
 -b|--min_baseq    Minimum base quality considered.
 -q|--min_mapq     Minimum mapping quality considered.
 -E|--exclude      Exclude all reads matching any of these SAM flags.
 -R|--region       Genomic region to consider.
 -r|--reference_fasta   FAI indexed fasta file of the reference genome.
 -i|--input        Input SAM/BAM/CRAM file. Indexed if used with -R.
 -m|--mod          DNA modification to consider. Format: 'CG+m.0' for methylation 'm' of C:s on position '0' 
                   of context 'CG' in forward strand. Can be given multiple times. If none given, 'CG+m.0' is used.
 --split_strand    Report modifications on either strand separately.
 --no_header       Don't output header text.
 --no_phase        Don't use HP and PS tags to separate phased reads.
 -@                Extra threads for reading input.
 -h|--help         Print this help text.

Using htslib version 1.16.
```

## Installing htslib

```
wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -jxvf htslib-1.21.tar.bz2
cd htslib-1.21
./configure
make
make install
```