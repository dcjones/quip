
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

    symb_t s, n, m;
    uint32_t z, x, y = ac->l;

    x = s = 0;
    ac->l >>= dist_length_shift;
    n = DISTSIZE;
    m = n >> 1;

    do {
        z = ac->l * D->xs[m].freq;
        if (z > ac->v) {
            n = m;
            y = z;
        }
        else {
            s = m;
            x = z;
        }

        m = (s + n) >> 1;

    } while(m != s);

    ac->v -= x;
    ac->l = y - x;

    if (ac->l < min_length) ac_renormalize_decoder(ac);

    D->xs[s].count += update_rate;

    if (--D->update_delay == 0) dfun(update)(D);

    return s;
}


symb_t dfun(decode)(ac_t* ac, dist_t* D)
{
    return dfun(decode2)(ac, D, 1);
}


void cdfun(init) (cond_dist_t* D, size_t n)
{
    D->n = n;
    D->xss   = malloc_or_die(n * sizeof(dist_t));
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


