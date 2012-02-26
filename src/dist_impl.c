
/* WARNING: DO NOT COMPILE THIS FILE DIRECTLY. You should compile "dist.c".
 */




void dfun(init)(dist_t* D)
{
    D->update_delay = DISTSIZE * update_delay_factor;

    /* initialize to pseudocounts of 1 */
    size_t i;
    for (i = 0; i < DISTSIZE; ++i) D->xs[i].count = 1;
    dfun(update)(D);
}



void dfun(set)(dist_t* D, const uint16_t* cs)
{
    size_t i;
    for (i = 0; i < DISTSIZE; ++i) {
        D->xs[i].count = cs[i];
    }

    dfun(update)(D);
}


void dfun(update)(dist_t* D)
{
    size_t i;

    /* rescale when we have exceeded the maximum count */
    uint32_t z = 0;
    for (i = 0; i < DISTSIZE; ++i) z += D->xs[i].count;
    if (z > max_count) {
        z = 0;
        for (i = 0; i < DISTSIZE; ++i) {
            D->xs[i].count = D->xs[i].count / 2 + 1;
            z += D->xs[i].count;
        }
    }


    /* update frequencies */
    const uint32_t scale = 0x80000000U / z;
    const uint32_t shift = 31 - dist_length_shift;
    uint32_t c = 0;

    for (i = 0; i < DISTSIZE; ++i) {
        D->xs[i].freq = (uint16_t) ((scale * c) >> shift);
        c += D->xs[i].count;
    }

    D->update_delay = DISTSIZE * update_delay_factor;
}


static void dfun(encode2)(ac_t* ac, dist_t* D, symb_t x, uint8_t update_rate)
{
    uint32_t b0 = ac->b;

    uint32_t u;
    if (x == DISTSIZE - 1) {
        u = (uint32_t) D->xs[x].freq * (ac->l >> dist_length_shift);
        ac->b += u;
        ac->l -= u;
    }
    else {
        u = (uint32_t) D->xs[x].freq * (ac->l >>= dist_length_shift);
        ac->b += u;
        ac->l = (uint32_t) D->xs[x + 1].freq * ac->l - u;
    }

    if (b0 > ac->b)         ac_propogate_carry(ac);
    if (ac->l < min_length) ac_renormalize_encoder(ac);

    D->xs[x].count += update_rate;

    if(!--D->update_delay) dfun(update)(D);
}

void dfun(encode)(ac_t* ac, dist_t* D, symb_t x)
{
    dfun(encode2)(ac, D, x, 1);
}


static symb_t dfun(decode2)(ac_t* ac, dist_t* D, uint8_t update_rate)
{
    if (ac->init_state) {
        ac->bufavail = ac->reader(ac->reader_data, ac->buf, ac->buflen);

        if (ac->bufavail < 4) {
            fprintf(stderr, "Malformed compressed data encountered.");
            exit(EXIT_FAILURE);
        }

        ac->v = ((uint32_t) ac->buf[0] << 24) | ((uint32_t) ac->buf[1] << 16) |
                ((uint32_t) ac->buf[2] << 8)  | ((uint32_t) ac->buf[3]);

        ac->l = max_length;

        ac->bufpos = 4;
        ac->init_state = false;
    }


    symb_t sym;
    uint32_t low_val, hi_val;


#if DISTSIZE == 41

    /* linear search */

    hi_val = ac->l;
    ac->l >>= dist_length_shift;

    __m128i fs, fs_low, fs_high, cmp1, cmp2, cmp3;
    const __m128i low_mask = _mm_set1_epi32(0xffff);
    const __m128i ls       = _mm_set1_epi32(ac->l);
    const __m128i vs       = _mm_set1_epi32(ac->v);
    const __m128i vs_low   = _mm_and_si128(vs, low_mask);
    const __m128i vs_high  = _mm_srli_epi32(vs, 16);

    symb_t low_sym = 0;
    int cmpmask;

    while (low_sym + 4 < DISTSIZE) {
        /* load the next four frequencies */
        fs = _mm_load_si128((__m128i*) (D->xs + low_sym));

        /* mask out the counts */
        fs = _mm_and_si128(fs, low_mask);

        /* multiply by l */
        fs = _mm_mullo_epi32(fs, ls);

        /* compare to ac->v */

        /* SSE provides no unsigned comparison operations so we have to
           simulate one using signed comparisons on 16-bit halves. */
        fs_high = _mm_srli_epi32(fs, 16);
        fs_low  = _mm_and_si128(fs, low_mask);
        cmp1 = _mm_cmpgt_epi32(fs_high, vs_high);
        cmp2 = _mm_cmpeq_epi32(fs_high, vs_high);
        cmp3 = _mm_cmpgt_epi32(fs_low, vs_low);
        cmpmask = _mm_movemask_epi8(_mm_or_si128(cmp1, _mm_and_si128(cmp2, cmp3)));

        if (cmpmask != 0) break;

        low_sym += 4;
    }

    /* we now know that the correct value is one of the next four */
    sym = low_sym;

    uint32_t vals[4];

    if (sym + 3 < DISTSIZE && (vals[3] = D->xs[sym + 3].freq * ac->l) <= ac->v) {
        low_val = vals[3];
        sym = low_sym + 3;
    }
    else if (sym + 2 < DISTSIZE && (vals[2] = D->xs[sym + 2].freq * ac->l) <= ac->v) {
        hi_val  = vals[3];
        low_val = vals[2];
        sym = low_sym + 2;
    }
    else if (sym + 1 < DISTSIZE && (vals[1] = D->xs[sym + 1].freq * ac->l) <= ac->v) {
        hi_val  = vals[2];
        low_val = vals[1];
        sym = low_sym + 1;
    }
    else if (sym + 0 < DISTSIZE && (vals[0] = D->xs[sym + 0].freq * ac->l) <= ac->v) {
        hi_val  = vals[1];
        low_val = vals[0];
        sym = low_sym + 0;
    }
    else {
        hi_val  = vals[0];
        low_val = D->xs[sym - 1].freq * ac->l;
        sym = low_sym - 1;
    } 

#if 0
    sym = DISTSIZE;
    hi_val = ac->l;

    ac->l >>= dist_length_shift;

    do {
        low_val = ac->l * D->xs[sym - 1].freq;

        if (low_val <= ac->v || sym == 1) break;

        hi_val = low_val;
        --sym;
    } while(true);

    --sym;
#endif

#else

    /* binary search */

    symb_t low_sym, mid_sym, hi_sym;
    uint32_t mid_val;

    low_val = 0;
    hi_val  = ac->l;

    low_sym = 0;
    hi_sym  = DISTSIZE;
    mid_sym = hi_sym / 2;

    ac->l >>= dist_length_shift;

    do {
        mid_val = ac->l * D->xs[mid_sym].freq;
        if (mid_val > ac->v) {
            hi_sym = mid_sym;
            hi_val = mid_val;
        }
        else {
            low_sym = mid_sym;
            low_val = mid_val;
        }

        mid_sym = (low_sym + hi_sym) / 2;

    } while(mid_sym != low_sym);

    sym = low_sym;

#endif

    ac->v -= low_val;
    ac->l = hi_val - low_val;

    if (ac->l < min_length) ac_renormalize_decoder(ac);

    D->xs[sym].count += update_rate;

    if (--D->update_delay == 0) dfun(update)(D);

    return sym;
}


symb_t dfun(decode)(ac_t* ac, dist_t* D)
{
    return dfun(decode2)(ac, D, 1);
}


void cdfun(init) (cond_dist_t* D, size_t n)
{
    D->n = n;
    D->xss = malloc_aligned_or_die(n * sizeof(dist_t));
    D->update_rate = 1;

    D->xss[0].update_delay = DISTSIZE * update_delay_factor;

    size_t i;
    for (i = 0; i < DISTSIZE; ++i) D->xss[0].xs[i].count = 1;
    dfun(update)(D->xss);


    for (i = 1; i < n; ++i) {
        memcpy(D->xss + i, D->xss, sizeof(dist_t));
    }
}


void cdfun(free) (cond_dist_t* D)
{
    if (D == NULL) return;

    free_aligned(D->xss);
}


void cdfun(set_update_rate) (cond_dist_t* D, uint8_t update_rate)
{
    D->update_rate = update_rate;
}


void cdfun(setall) (cond_dist_t* D, const uint16_t* cs)
{
    size_t i;
    for (i = 0; i < DISTSIZE; ++i) {
        D->xss[0].xs[i].count = cs[i];
    }

    dfun(update)(D->xss);

    for (i = 1; i < D->n; ++i) {
        memcpy(D->xss + i, D->xss, sizeof(dist_t));
    }
}


void cdfun(setone) (cond_dist_t* D, const uint16_t* cs, size_t i)
{
    size_t j;
    for (j = 0; j < DISTSIZE; ++j) {
        D->xss[i].xs[j].count = cs[j];
    }

    dfun(update)(D->xss + i);
}


void cdfun(encode)(ac_t* ac, cond_dist_t* D, uint32_t y, symb_t x)
{
    dfun(encode2)(ac, D->xss + y, x, D->update_rate);
}


symb_t cdfun(decode)(ac_t* ac, cond_dist_t* D, uint32_t y)
{
    return dfun(decode2)(ac, D->xss + y, D->update_rate);
}


