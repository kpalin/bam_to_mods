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

typedef struct
{
    int modified_count;
    int canonical_count;
    int other_count;
} _mod_count_t;

typedef struct
{
    int pos;
    int modified_base;
    int canonical_base;
    int mod_isrev;
    int read_isrev;
    int ps;
    int hp;
} mod_key;
#include "../klib/kstring.h"

char *mod_to_key(mod_key *k)
{
    kstring_t r;
    ksprintf(&r, "%d %d %d %d %d %d", k->modified_base, k->canonical_base,
             k->mod_isrev, k->read_isrev, k->ps, k->hp);
    return ks_release(&r);
}

KHASH_MAP_INIT_STR(mod_count_h, _mod_count_t);

void put_mod(khash_t(mod_count_h) * mod_ret, mod_key *k, _mod_count_t *c)
{
    char *ks = mod_to_key(k);
    // kh_put(mod_count_h, mod_ret,ks,c);

    khint_t k_idx;

    int absent;
    k_idx = kh_put(mod_count_h, mod_ret, ks, &absent);
    if (absent)
        kh_key(mod_ret, k_idx) = strdup(ks);
    // else, the key is not touched; we do nothing

    free(ks);
}

void add_mod(khash_t(mod_count_h) * mod_ret, mod_key *k, _mod_count_t *c)
{
    char *ks = mod_to_key(k);
    // kh_put(mod_count_h, mod_ret,ks,c);

    khint_t k_idx = kh_get(mod_count_h, mod_ret, ks);

    if ((k_idx == kh_end(mod_ret))) // Missing
    {
        int absent;
        k_idx = kh_put(mod_count_h, mod_ret, ks, &absent);
        kh_value(mod_ret, k_idx) = *c;
        if (absent)
            kh_key(mod_ret, k_idx) = strdup(ks);
    }
    else
    {
        kh_value(mod_ret, k_idx).modified_count += c->modified_count;
        kh_value(mod_ret, k_idx).canonical_count += c->canonical_count;
        kh_value(mod_ret, k_idx).other_count += c->other_count;
    }
    // else, the key is not touched; we do nothing

    free(ks);
}

void destroy_mod(khash_t(mod_count_h) * mod_ret)
{
    // printf("# of distinct words: %d\n", kh_size(h));
    //  IMPORTANT: free memory allocated by strdup() above
    khint_t k;
    for (k = 0; k < kh_end(mod_ret); ++k)
        if (kh_exist(mod_ret, k))
            free((char *)kh_key(mod_ret, k));
    kh_destroy(mod_count_h, mod_ret);
}

// Report a line of pileup, including base modifications inline with
// the sequence (including qualities), as [<strand><dir><qual>...]
void process_mod_pileup0(sam_hdr_t *h, const bam_pileup1_t *p,
                         int tid, int pos, int n, char ref_base,
                         khash_t(mod_count_h) * mod_ret)
{

    // printf("%s\t%d\t%c\n", sam_hdr_tid2name(h, tid), pos, ref_base);
    int i;
    _mod_count_t mod_counts = {
        .modified_count = 0,
        .canonical_count = 0,
        .other_count = 0};

    // TODO: Fix this assumption
    // int read_rev = ref_base != canonical_base;

    for (i = 0; i < n; i++, p++)
    {
        mod_key mod_id = {.pos = pos, 0, 0, 0, 0, -1, -1};
        uint8_t *seq = bam_get_seq(p->b);
        uint8_t *qual = bam_get_qual(p->b);
        unsigned char q_base = seq_nt16_str[bam_seqi(seq, p->qpos)];
        uint8_t base_qual = qual[p->qpos];

        if (base_qual < MIN_BQ)
        { // Ignore bases with bad quality.
            continue;
        }

        // if (bam_is_rev(p->b) != read_rev)
        // {
        //     continue;
        // }

        uint8_t *_hp_tag = bam_aux_get(p->b, "HP");

        if (_hp_tag != NULL)
        {
            mod_id.hp = bam_aux2i(_hp_tag);
            mod_id.ps = bam_aux2i(bam_aux_get(p->b, "PS"));
        }

        if (p->is_del || q_base != ref_base)
        {
            // Only count the number of deletions and mismatches.
            // mod_id.canonical_base = q_base;
            // mod_id.modified_base = 0;
            // mod_id.mod_isrev = 0;
            // mod_id.read_isrev = bam_is_rev(p->b);

            mod_counts.other_count++;
            continue;
        }
        else
        {
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
                    if (mod[j].modified_base == 'm' /*mod_counts.modified_base*/)
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
            { // There is no modifications called at this locus
                mod_counts.canonical_count++;
                // putchar(c);
            }
        }
    }
    if ((mod_counts.modified_count + mod_counts.canonical_count) > 0)
    {
        // printf("%s\t%d\t%c\t%d\t%d\t%c\t%c\t",
        //        sam_hdr_tid2name(h, tid), pos, ref_base,
        //        mod_counts.modified_count + mod_counts.canonical_count + mod_counts.other_count,
        //        PS, HP, "+-"[read_rev]);

        // if (mod_counts.modified_base < 0)
        //     // ChEBI
        //     printf("%d\t", -mod_counts.modified_base);
        // else
        //     printf("%c\t", mod_counts.modified_base);

        // double mod_prop = mod_counts.modified_count * 1.0 / (mod_counts.modified_count + mod_counts.canonical_count);
        // printf("%d\t%d\t%d\t%.3g\n", mod_counts.canonical_count, mod_counts.modified_count, mod_counts.other_count, mod_prop);
        // putchar('\n');
    }
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s ref.fasta input.cram\n", argv[0]);
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

    // khash_t(mod_count_h) *hap_h = kh_init(mod_count_h);

    const bam_pileup1_t *p;
    int tid, pos, n;
    printf("Mod0\n");
    int prev_tid = -1;
    hts_pos_t ref_len = 0;
    char *ref_seq = NULL;

    khash_t(mod_count_h) * mod_ch;
    mod_ch = kh_init(mod_count_h);

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
            // fprintf(stderr, "Pos: %d\n", pos);
            process_mod_pileup0(h, p, tid, pos, n, ref_seq[pos], mod_ch);
            // kh_clear(mod_count_h, mod_ch);
        }
    }
    bam_plp_destroy(iter);

    sam_close(in);
    bam_destroy1(b);
    sam_hdr_destroy(h);
    free(ref_fai);
    return 0;
}
