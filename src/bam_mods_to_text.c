/*

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

//#include <config.h>

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "../klib/khash.h"
#include "../klib/kvec.h"
typedef struct
{
    samFile *fp;
    sam_hdr_t *fp_hdr;
} plp_dat;

static int readaln(void *data, bam1_t *b)
{
    plp_dat *g = (plp_dat *)data;
    int ret;

    while (1)
    {
        ret = sam_read1(g->fp, g->fp_hdr, b);
        if (ret < 0)
            break;
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))
            continue;
        break;
    }

    return ret;
}

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

// Initialise and destroy the base modifier state data. This is called
// as each new read is added or removed from the pileups.
int pileup_cd_create(void *data, const bam1_t *b, bam_pileup_cd *cd)
{
    hts_base_mod_state *m = hts_base_mod_state_alloc();
    bam_parse_basemod(b, m);
    cd->p = m;
    return 0;
}

int pileup_cd_destroy(void *data, const bam1_t *b, bam_pileup_cd *cd)
{
    hts_base_mod_state_free(cd->p);
    return 0;
}

#define MIN_BQ 10
#define MIN_MODQ 127

typedef union
{
    struct
    {
        uint8_t HPp1;
        uint64_t PS : 56;
    } hap;
    uint64_t data;
} hap_t;

struct _mod_count_t
{
    int modified_base;
    int canonical_base;
    int strand;
    int modified_count;
    int canonical_count;
    int other_count;
};

typedef struct _mod_count_t mod_count_t;

KHASH_MAP_INIT_INT64(mod_count_h, mod_count_t);

// Report a line of pileup, including base modifications inline with
// the sequence (including qualities), as [<strand><dir><qual>...]
void process_mod_pileup0(sam_hdr_t *h, const bam_pileup1_t *p,
                         int tid, int pos, int n, char ref_base, int modified_base,
                         int canonical_base, int read_rev,
                         char const HP, int const PS)
{

    // printf("%s\t%d\t%c\n", sam_hdr_tid2name(h, tid), pos, ref_base);
    int i;

    mod_count_t mod_counts = {
        .modified_base = modified_base,
        .canonical_base = canonical_base,
        .strand = '+',
        .modified_count = 0,
        .canonical_count = 0,
        .other_count = 0,
    };

    for (i = 0; i < n; i++, p++)
    {

        uint8_t *seq = bam_get_seq(p->b);
        uint8_t *qual = bam_get_qual(p->b);
        unsigned char c = seq_nt16_str[bam_seqi(seq, p->qpos)];
        uint8_t base_qual = qual[p->qpos];

        // printf("%d%c", base_qual, c);
        if (base_qual < MIN_BQ)
        { // Ignore bases with bad quality.
            continue;
        }

        if (bam_is_rev(p->b) != read_rev)
        {
            continue;
        }
        if (p->is_del)
        {
            mod_counts.other_count++;
            continue;
        }

        // // Get phase information, HPp1==0 is unphased, otherwise phased to haplotype HP where HPp1=HP+1

        // hap_t hap = {.data = 0};
        // uint8_t *_hp_tag = bam_aux_get(p->b, "HP");
        // if (_hp_tag != NULL)
        // {
        //     hap.hap.HPp1 = 1 + bam_aux2i(_hp_tag);
        //     hap.hap.PS = bam_aux2i(bam_aux_get(p->b, "PS"));
        // }
        // mod_count_t const *cnts;
        // mod_count_t new_cnts;
        // // khint_t k = kh_get_mod_count_h(cnts, hap.data);

        if (c != ref_base)
        {
            mod_counts.other_count++;
            continue;
        }
        // putchar('\n');
        // continue;
        //  Simple mod detection; assumes at most 5 mods
        hts_base_mod_state *m = p->cd.p;
        hts_base_mod mod[5];
        int nm;
        // putchar("+-"[bam_is_rev(p->b)]);
        if ((nm = bam_mods_at_qpos(p->b, p->qpos, m, mod, 5)) > 0)
        {
            int j;

            // putchar('[');
            for (j = 0; j < nm && j < 5; j++)
            {
                if (mod[j].modified_base == mod_counts.modified_base)
                {
                    if (mod[j].qual >= MIN_MODQ)
                    {
                        mod_counts.modified_count++;
                    }
                    else
                    {
                        mod_counts.canonical_count++;
                    }
                }
                // printf("%c%c%d", "+-"[mod[j].strand],
                //        mod[j].modified_base, mod[j].qual);
            }
            // putchar(']');
        }
        else
        {
            mod_counts.canonical_count++;
            // putchar(c);
        }
    }
    if ((mod_counts.modified_count + mod_counts.canonical_count) > 0)
    {
        printf("%s\t%d\t%c\t%d\t%d\t%c\t%c\t",
               sam_hdr_tid2name(h, tid), pos, ref_base,
               mod_counts.modified_count + mod_counts.canonical_count + mod_counts.other_count,
               PS, HP, "+-"[read_rev]);

        if (mod_counts.modified_base < 0)
            // ChEBI
            printf("%d\t", -mod_counts.modified_base);
        else
            printf("%c\t", mod_counts.modified_base);

        double mod_prop = mod_counts.modified_count * 1.0 / (mod_counts.modified_count + mod_counts.canonical_count);
        printf("%d\t%d\t%d\t%.3g\n", mod_counts.canonical_count, mod_counts.modified_count, mod_counts.other_count, mod_prop);
        // putchar('\n');
    }
}

KHASH_SET_INIT_STR(str)

typedef struct
{
    int ps;
    int hp;
} _hap_info;
typedef struct
{
    int modified_base;
    int canonical_base;
    int strand;
} _mod_info;
int find_phase_info(sam_hdr_t *h, const bam_pileup1_t *p, void *_haps, void *_mods,
                    int n)
{
    int i;
    kvec_t(_hap_info) *haps = _haps;
    kvec_t(_mod_info) *mods = _mods;
    _hap_info ch;
    char s[4096]; // max string length: 4095 characters
    khash_t(str) * hset;
    khint_t k;
    h = kh_init(str);

    ch.hp = 'N';
    ch.ps = 0;
    kv_push(_hap_info, *haps, ch);
    for (i = 0; i < n; i++, p++)
    {

        uint8_t *_hp_tag = bam_aux_get(p->b, "HP");
        if (_hp_tag != NULL)
        {
            int j;
            int do_add = 1;
            ch.hp = bam_aux2i(_hp_tag);
            ch.ps = bam_aux2i(bam_aux_get(p->b, "PS"));
            for (j = 0; j < kv_size(*haps); j++)
            {
                do_add &= !memcmp(&kv_A(*haps, j), &ch, sizeof(ch));
            }
            if (do_add)
            {
                kv_push(_hap_info, *haps, ch);
            }
        }

        hts_base_mod_state *m = p->cd.p;
        hts_base_mod mod[5];
        int nm;
        _mod_info mi;

        if ((nm = bam_mods_at_qpos(p->b, p->qpos, m, mod, 5)) > 0)
        {
            int j;

            // putchar('[');
            for (j = 0; j < nm && j < 5; j++)
            {
                int k;
                int do_add = 1;
                mi.canonical_base = mod[j].canonical_base;
                mi.modified_base = mod[j].modified_base;
                mi.strand = mod[j].strand;
                for (k = 0; j < kv_size(*mods); k++)
                {
                    do_add &= !memcmp(&kv_A(*mods, k), &ch, sizeof(mi));
                }
                if (do_add)
                {
                    kv_push(_mod_info, *mods, mi);
                }
            }
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        exit(1);
    }
    // First argument: reference genome
    faidx_t *ref_fai = fai_load(argv[1]);
    if (!ref_fai)
    {
        fprintf(stderr, "Can't open reference genome fasta '%s'.", argv[1]);
        exit(1);
    }
    argc--;
    argv++;

    samFile *in = sam_open(argc > 1 ? argv[1] : "-", "r");
    bam1_t *b = bam_init1();
    sam_hdr_t *h = sam_hdr_read(in);

    // Pileup iterator with constructor/destructor to parse base mod tags
    plp_dat dat = {
        .fp = in,
        .fp_hdr = h,
    };
    bam_plp_t iter = bam_plp_init(readaln, &dat);
    bam_plp_constructor(iter, pileup_cd_create);
    bam_plp_destructor(iter, pileup_cd_destroy);

    khash_t(mod_count_h) *hap_h = kh_init(mod_count_h);

    const bam_pileup1_t *p;
    int tid, pos, n;
    printf("Mod0\n");
    int prev_tid = -1;
    hts_pos_t ref_len = 0;
    char *ref_seq = NULL;
    kvec_t(_hap_info) haps;
    kvec_t(_mod_info) mods;
    kv_init(mods);
    kv_init(haps);

    while ((p = bam_plp_auto(iter, &tid, &pos, &n)) != 0)
    {
        if (tid != prev_tid)
        {
            if (ref_seq)
            {
                free(ref_seq);
            }
            ref_seq = fai_fetch64(ref_fai, sam_hdr_tid2name(dat.fp_hdr, tid), &ref_len);
            if (!ref_seq)
            {
                fprintf(stderr, "Couldn't fetch reference!");
                exit(1);
            }
            prev_tid = tid;
        }
        if (pos >= ref_len)
        {
            fprintf(stderr, "Trying to access reference beyond end!");
            exit(1);
        }
        // Only output CpG sites.

        if ((pos < ref_len && ref_seq[pos] == 'C' && ref_seq[pos + 1] == 'G') || (pos > 0 && ref_seq[pos - 1] == 'C' && ref_seq[pos] == 'G'))
        {
            find_phase_info(h, p, &haps, &mods, n);
            printf("#haps: %d #mods: %d\n", kv_size(haps), kv_size(mods));
            process_mod_pileup0(h, p, tid, pos, n, ref_seq[pos], 'm', 'C', ref_seq[pos] == 'G', '*', 0);
            process_mod_pileup0(h, p, tid, pos, n, ref_seq[pos], 'h', 'C', ref_seq[pos] == 'G', '*', 0);
            // kh_clear(mod_count_h, hap_h);
            while (kv_size(haps) > 0)
                kv_pop(haps);
        }
    }
    bam_plp_destroy(iter);

    sam_close(in);
    bam_destroy1(b);
    sam_hdr_destroy(h);
    free(ref_fai);
    return 0;
}
