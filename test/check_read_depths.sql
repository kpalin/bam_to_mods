create table bam_mods as
select
    *
from
    read_csv('test/tmp/A+a.0.out');

create table depths as (
    SELECT
        *,
        '1' as haplotype,
        pos as 'end'
    from
        read_csv(
            "data/fibreseq_demo_pacbio_hp1.depth.gz",
            names = ["#chromosome","pos","depth"]
        )
    UNION
    SELECT
        *,
        '2' as haplotype,
        pos as 'end'
    from
        read_csv(
            "data/fibreseq_demo_pacbio_hp2.depth.gz",
            names = ["#chromosome","pos","depth"]
        )
);
.mode csv
select
    count(*),
(called_reads + uncalled_reads + mismatch_reads - depth) as gap
from
    bam_mods a,
    depths b
where
    b."end" = a."end"
    and a.haplotype = b.haplotype
    
group by 
    gap;


