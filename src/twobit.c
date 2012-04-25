
#include "twobit.h"
#include "kmer.h"
#include "misc.h"
#include <string.h>

struct twobit_t_
{
    size_t len; /* length of stored sequence */
    size_t n;   /* space (number of kmers) allocated in seq */
    kmer_t* seq;
};


/* the number of kmer_t elements needed to represent a sequence of the given
 * length */
static size_t kmers_needed(size_t len)
{
    if (len == 0) return 0;
    else          return (len - 1) / (4 * sizeof(kmer_t)) + 1;
}


void twobit_reserve(twobit_t* s, size_t seqlen)
{
    if (s->n < kmers_needed(s->len + seqlen)) {
        size_t new_n = s->n;
        while (new_n < kmers_needed(s->len + seqlen)) new_n *= 2;
        s->seq = realloc_or_die(s->seq, new_n * sizeof(kmer_t));
        memset(s->seq + s->n, 0, (new_n - s->n) * sizeof(kmer_t));
        s->n = new_n;
    }
}


void twobit_free_reserve(twobit_t* s)
{
    s->n = kmers_needed(s->len);
    s->seq = realloc_or_die(s->seq, s->n * sizeof(kmer_t));
}


twobit_t* twobit_alloc_n(size_t len)
{
    twobit_t* s = malloc_or_die(sizeof(twobit_t));
    s->len = 0;
    s->n   = kmers_needed(len);
    s->seq = malloc_or_die(s->n * sizeof(kmer_t));
    memset(s->seq, 0, s->n * sizeof(kmer_t));

    return s;
}


twobit_t* twobit_alloc()
{
    return twobit_alloc_n(512);
}


void twobit_free(twobit_t* s)
{
    if (s != NULL) {
        free(s->seq);
        free(s);
    }
}


twobit_t* twobit_dup(const twobit_t* s)
{
    twobit_t* t = twobit_alloc_n(s->len);
    t->len = s->len;
    memcpy(t->seq, s->seq, kmers_needed(t->len) * sizeof(kmer_t));
    return t;
}


void twobit_clear(twobit_t* s)
{
    s->len = 0;
    memset(s->seq, 0, s->n * sizeof(kmer_t));
}


size_t twobit_len(const twobit_t* s)
{
    return s->len;
}

void twobit_copy(twobit_t* s, const char* seqstr)
{
    twobit_clear(s);
    twobit_append(s, seqstr);
}


void twobit_copy_n(twobit_t* s, const char* seqstr, size_t seqlen)
{
    twobit_clear(s);
    twobit_append_n(s, seqstr, seqlen);
}


void twobit_append(twobit_t* s, const char* seqstr)
{
    twobit_append_n(s, seqstr, strlen(seqstr));
}


void twobit_append_n(twobit_t* s, const char* seqstr, size_t seqlen)
{
    twobit_reserve(s, seqlen);

    kmer_t c;
    size_t idx, off;
    size_t i;
    for (i = 0; i < seqlen; ++i) {
        c = chartokmer[(uint8_t) seqstr[i]];
        if (c > 3) continue;

        idx = (s->len + i) / (4 * sizeof(kmer_t));
        off = (s->len + i) % (4 * sizeof(kmer_t));

        s->seq[idx] = (s->seq[idx] & ~((kmer_t) 0x3 << (2 * off))) | (c << (2 * off));
    }

    s->len += seqlen;
}


void twobit_append_char(twobit_t* s, char a)
{
    twobit_reserve(s, 1);
    kmer_t c = chartokmer[(uint8_t) a];

    size_t idx = s->len / (4 * sizeof(kmer_t));
    size_t off = s->len % (4 * sizeof(kmer_t));

    s->seq[idx] = (s->seq[idx] & ~((kmer_t) 0x3 << (2 * off))) | (c << (2 * off));
    s->len++;
}


void twobit_append_kmer(twobit_t* s, kmer_t x, size_t k)
{
    /* expand, if more space is needed */
    if (s->n < kmers_needed(s->len + k)) {
        while (s->n < kmers_needed(s->len + k)) s->n *= 2;
        s->seq = realloc_or_die(s->seq, s->n * sizeof(kmer_t));
    }

    size_t idx, off;
    size_t i;
    for (i = 0; i < k; ++i) {

        idx = (s->len + i) / (4 * sizeof(kmer_t));
        off = (s->len + i) % (4 * sizeof(kmer_t));

        s->seq[idx] = (s->seq[idx] & ~((kmer_t) 0x3 << (2 * off))) | ((x & 0x3) << (2 * off));
        x >>= 2;
    }

    s->len += k;
}



void twobit_append_twobit(twobit_t* s, const twobit_t* t)
{
    /* expand, if more space is needed */
    if (s->n < kmers_needed(s->len + t->len)) {
        while (s->n < kmers_needed(s->len + t->len)) s->n *= 2;
        s->seq = realloc_or_die(s->seq, s->n * sizeof(kmer_t));
    }


    size_t idx, off;
    size_t i;
    for (i = 0; i < t->len; ++i) {

        idx = (s->len + i) / (4 * sizeof(kmer_t));
        off = (s->len + i) % (4 * sizeof(kmer_t));

        s->seq[idx] = (s->seq[idx] & ~((kmer_t) 0x3 << (2 * off))) | (twobit_get(t, i) << (2 * off));
    }

    s->len += t->len;
}




void twobit_reverse(twobit_t* s)
{
    kmer_t x;
    int i, j;
    for (i = 0, j = s->len - 1; i < j; ++i, --j) {
        x = twobit_get(s, i);
        twobit_set(s, i, twobit_get(s, j));
        twobit_set(s, j, x);
    }
}



void twobit_setc(twobit_t* s, size_t i, char seqc)
{
    size_t idx = i / (4 * sizeof(kmer_t));
    size_t off = i % (4 * sizeof(kmer_t));

    kmer_t c = chartokmer[(uint8_t) seqc]; 
    if (c > 3) return;

    s->seq[idx] = (s->seq[idx] & ~((kmer_t) 0x3 << (2 * off))) | (c << (2 * off));

}


void twobit_set(twobit_t* s, size_t i, kmer_t x)
{
    size_t idx = i / (4 * sizeof(kmer_t));
    size_t off = i % (4 * sizeof(kmer_t));

    s->seq[idx] = (s->seq[idx] & ~((kmer_t) 0x3 << (2 * off))) | (x << (2 * off));

}


kmer_t twobit_get(const twobit_t* s, size_t i)
{
    size_t idx = i / (4 * sizeof(kmer_t));
    size_t off = i % (4 * sizeof(kmer_t));

    return (s->seq[idx] >> (2 * off)) & 0x3;
}

kmer_t twobit_get_kmer_rev(const twobit_t* s, size_t i, size_t k)
{
    kmer_t x = 0;
    size_t j;
    for (j = i; j < i + k; ++j) {
        x = (x << 2) | twobit_get(s, j);
    }

    return x;
}

kmer_t twobit_get_kmer(const twobit_t* s, size_t i, size_t k)
{
    kmer_t x = 0;
    size_t j;
    for (j = 0; j < k; ++j) {
        x |= twobit_get(s, i + j) << (2 * j);
    }

    return x;
}


void twobit_print(const twobit_t* s, FILE* fout)
{
    size_t i;
    for (i = 0; i < s->len; ++i) {
        fputc(kmertochar[twobit_get(s, i)], fout);
    }
}

void twobit_print_stdout(const twobit_t* s)
{
    twobit_print(s, stdout);
    fputc('\n', stdout);
    fflush(stdout);
}



int twobit_cmp(const twobit_t* a, const twobit_t* b)
{
    if      (a->len < b->len) return -1;
    else if (a->len > b->len) return 1;

    return memcmp(a->seq, b->seq, kmers_needed(a->len) * sizeof(kmer_t));
}



void twobit_revcomp(twobit_t* dest, const twobit_t* src)
{
    twobit_reserve(dest, src->len);

    size_t i;
    for (i = 0; i < src->len; ++i) {
        twobit_set(dest, src->len - i - 1, kmer_comp(twobit_get(src, i), 1));
    }

    dest->len = src->len;
}



/*
 * Paul Hsieh's SuperFastHash
 * http://www.azillionmonkeys.com/qed/hash.html
 *
 * TODO: copyright
 */

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
    || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
        +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

uint32_t twobit_hash(const twobit_t* s)
{
    if (s->len == 0) return 0;
    size_t len = kmers_needed(s->n) * sizeof(kmer_t);
    uint8_t* data = (uint8_t*) s->seq;

    uint32_t hash = len, tmp;
    int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}



uint32_t twobit_mismatch_count(const twobit_t* subject,
                               const twobit_t* query,
                               size_t spos, uint32_t max_miss)
{
    /* This is made a bit tricky because the query's
     * limbs are probably not aligned with the subjects.
     * We have to compare each query limb to the end of
     * one subject limb and the beginning of the following.
     */

    size_t m = query->len;
    uint32_t mismatches = 0;

    const size_t k = 4 * sizeof(kmer_t);
    const size_t l = spos % k;

    size_t subj_idx  = spos / k;
    size_t query_idx = 0;

    kmer_t x, y;

    x = subject->seq[subj_idx] >> (2 * l);
    y = query->seq[query_idx];

    size_t i = 0; /* nucleotide position within query sequence */
    size_t j;     /* position within the 'section' */
    while (i < m) {
        for (j = 0; j < (k - l) && i < m; ++j, ++i) {
            if ((x & 0x3) != (y & 0x3)) ++mismatches;
            x >>= 2;
            y >>= 2;
        }
        if (i >= m || mismatches >= max_miss) break;

        x = subject->seq[++subj_idx];
        for (j = 0; j < l && i < m; ++j, ++i) {
            if ((x & 0x3) != (y & 0x3)) ++mismatches;
            x >>= 2;
            y >>= 2;
        }
        if (i >= m || mismatches >= max_miss) break;

        y = query->seq[++query_idx];
    }

    return mismatches;
}




