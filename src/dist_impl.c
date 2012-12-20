
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
            D->xs[i].count /= 2;
            if (D->xs[i].count == 0) D->xs[i].count = 1;
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


void dfun(encode2)(ac_t* ac, dist_t* D, symb_t x, uint8_t update_rate)
{
    prefetch(D, 1, 0);

    uint32_t b0 = ac->b;

    uint32_t u;
    if (expect(x == DISTSIZE - 1, 0)) {
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
    prefetch(D, 1, 0);

    symb_t low_sym, mid_sym, hi_sym;
    uint32_t low_val, mid_val, hi_val;

    low_val = 0;
    hi_val  = ac->l;

    low_sym = 0;
    hi_sym  = DISTSIZE;
    mid_sym = hi_sym / 2;

    ac->l >>= dist_length_shift;

    /* binary search */
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

    ac->v -= low_val;
    ac->l = hi_val - low_val;

    if (ac->l < min_length) ac_renormalize_decoder(ac);

    D->xs[low_sym].count += update_rate;

    if (--D->update_delay == 0) dfun(update)(D);

    return low_sym;
}


symb_t dfun(decode)(ac_t* ac, dist_t* D)
{
    return dfun(decode2)(ac, D, 1);
}


void cdfun(init) (cond_dist_t* D, size_t n)
{
    if (n == 0) n = 1;

    D->n = n;
    D->xss   = malloc_or_die(n * sizeof(dist_t));
    D->update_rate = 1;

    /*D->xss[0].update_delay = DISTSIZE * update_delay_factor;*/
    D->xss[0].update_delay = 1;

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

    free(D->xss);
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


