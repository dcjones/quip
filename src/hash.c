
#include "hash.h"
#include "misc.h"
#include <string.h>

static const size_t kmer_hash_initial_size = 1024;
static const size_t kmer_hash_max_load     = 5.0;


struct kmer_hash_t_
{
    size_t n;     /* number of slots */
    size_t m;     /* number of pairs stored */
    size_t m_max; /* number of stored pairs before we expand */

    size_t* slot_sizes;
    kmer_count_pair_t** slots;
};


kmer_hash_t* kmer_hash_alloc()
{
    kmer_hash_t* H = malloc_or_die(sizeof(kmer_hash_t));
    H->n = kmer_hash_initial_size;
    H->m = 0;
    H->m_max = (size_t) (kmer_hash_max_load * (double) H->n);

    H->slot_sizes = malloc_or_die(H->n * sizeof(size_t));
    memset(H->slot_sizes, 0, H->n * sizeof(size_t));

    H->slots = malloc_or_die(H->n * sizeof(kmer_count_pair_t*));
    memset(H->slots, 0, H->n * sizeof(kmer_count_pair_t*));

    return H;
}


void kmer_hash_free(kmer_hash_t* H)
{
    size_t i;
    for (i = 0; i < H->n; ++i) {
        free(H->slots[i]);
    }
    free(H->slots);

    free(H->slot_sizes);
    free(H);
}


size_t kmer_hash_size(kmer_hash_t* H)
{
    return H->m;
}



static void kmer_hash_expand(kmer_hash_t* H)
{
    size_t new_n = 2 * H->n;
    size_t* slot_sizes = malloc_or_die(new_n * sizeof(size_t));
    memset(slot_sizes, 0, new_n * sizeof(size_t));

    /* figure out slot sizes in the new table */
    uint64_t h;
    size_t i;
    size_t slot_size;
    kmer_count_pair_t *slot, *slot0;
    for (i = 0; i < H->n; ++i) {
        slot_size = H->slot_sizes[i];
        slot0 = H->slots[i];
        for (slot = slot0; (size_t) (slot - slot0) < slot_size; ++slot) {
            h = kmer_hash(slot->x) % new_n;
            ++slot_sizes[h];
        }
    }


    /* allocate slots */
    kmer_count_pair_t** slots = malloc_or_die(new_n * sizeof(kmer_count_pair_t*));
    for (i = 0; i < new_n; ++i) {
        if (slot_sizes[i] > 0) {
            slots[i] = malloc_or_die(slot_sizes[i] * sizeof(kmer_count_pair_t));
        }
        else slots[i] = NULL;
    }


    /* rehash values */
    kmer_count_pair_t** slots_next = malloc_or_die(new_n * sizeof(kmer_count_pair_t*));
    memcpy(slots_next, slots, new_n * sizeof(kmer_count_pair_t*));

    for (i = 0; i < H->n; ++i) {
        slot_size = H->slot_sizes[i];
        slot0 = H->slots[i];
        for (slot = slot0; (size_t) (slot - slot0) < slot_size; ++slot) {
            h = kmer_hash(slot->x) % new_n;

            slots_next[h]->x     = slot->x;
            slots_next[h]->count = slot->count;

            ++slots_next[h];
        }
    }

    free(slots_next);

    for (i = 0; i < H->n; ++i) free(H->slots[i]);
    free(H->slots);
    H->slots = slots;

    free(H->slot_sizes);
    H->slot_sizes = slot_sizes;

    H->n = new_n;
    H->m_max = (size_t) (kmer_hash_max_load * (double) H->n);
}


void kmer_hash_put(kmer_hash_t* H, kmer_t x, unsigned int count)
{
    /* if we are at capacity, preemptively expand */
    if (H->m >= H->m_max) kmer_hash_expand(H);

    uint64_t h = kmer_hash(x) % H->n;
    kmer_count_pair_t *slot, *slot0;
    slot0 = H->slots[h];
    size_t slot_size = H->slot_sizes[h];

    for (slot = slot0; (size_t) (slot - slot0) < slot_size; ++slot) {
        if (slot->x == x) {
            slot->count = count;
            return;
        }
    }


    /* kmer was not present, we must insert */
    H->slot_sizes[h] += 1;
    slot0 = realloc_or_die(slot0, H->slot_sizes[h] * sizeof(kmer_count_pair_t));
    slot0[H->slot_sizes[h] - 1].x     = x;
    slot0[H->slot_sizes[h] - 1].count = count;
    H->slots[h] = slot0;
    ++H->m;
}



unsigned int kmer_hash_get(kmer_hash_t* H, kmer_t x)
{
    uint64_t h = kmer_hash(x) % H->n;
    kmer_count_pair_t *slot, *slot0;
    slot0 = H->slots[h];
    size_t slot_size = H->slot_sizes[h];

    for (slot = slot0; (size_t) (slot - slot0) < slot_size; ++slot) {
        if (slot->x == x) {
            return slot->count;
        }
    }

    return 0;
}


int kmer_count_pair_cmp(const void* a, const void* b)
{
    unsigned int cnta = ((kmer_count_pair_t*) a)->count;
    unsigned int cntb = ((kmer_count_pair_t*) b)->count;

    if (cnta <= cntb) return cntb - cnta;
    else              return -(int) (cnta - cntb);
}


kmer_count_pair_t* kmer_hash_dump_sorted(kmer_hash_t* H)
{
    kmer_count_pair_t* dump = malloc_or_die(H->m * sizeof(kmer_count_pair_t));

    size_t i, j;
    kmer_count_pair_t *slot, *slot0;
    size_t slot_size;
    for (i = 0, j = 0; i < H->n; ++i) {
        slot0     = H->slots[i];
        slot_size = H->slot_sizes[i];
        for (slot = slot0; (size_t) (slot - slot0) < slot_size; ++slot, ++j) {
            dump[j].x     = slot->x;
            dump[j].count = slot->count;
        }
    }

    qsort(dump, H->m, sizeof(kmer_count_pair_t), kmer_count_pair_cmp);

    return dump;
}



