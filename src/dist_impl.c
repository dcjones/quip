
/* WARNING: DO NOT COMPILE THIS FILE DIRECTLY. You should compile "dist.c".
 */




void dfun(init)(dist_t* D, bool decode)
{
    if (decode && DISTSIZE > 16) {
        D->dec = malloc_or_die((dec_size + 6) * sizeof(uint16_t));
    }
    else {
        D->dec = NULL;
    }

    D->update_delay = DISTSIZE * update_delay_factor;

    /* initialize to pseudocounts of 1 */
    size_t i;
    for (i = 0; i < DISTSIZE; ++i) D->xs[i].count = 1;
    dfun(update)(D);
}



void dfun(free)(dist_t* D)
{
    if (D == NULL) return;
    free(D->dec);
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
    uint32_t j, w, c = 0;

    for (i = 0; i < DISTSIZE; ++i) {
        D->xs[i].freq = (uint16_t) ((scale * c) >> shift);
        c += D->xs[i].count;
    }


    /* update decoder table */
    if (D->dec) {
        for (i = 0, j = 0; i < DISTSIZE; ++i) {
            w = D->xs[i].freq >> dec_shift;
            while (j < w) D->dec[++j] = i - 1;
        }
        D->dec[0] = 0;
        while (j < dec_size) D->dec[++j] = DISTSIZE - 1;
    }


    D->update_delay = DISTSIZE * update_delay_factor;
}



void dfun(encode)(ac_t* ac, dist_t* D, symb_t x)
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

    D->xs[x].count++;
    D->use_count++;

    if(--D->update_delay == 0) dfun(update)(D);
}


symb_t dfun(decode)(ac_t* ac, dist_t* D)
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

    if (D->dec) {
        uint32_t dv = ac->v / (ac->l >>= dist_length_shift);
        uint32_t t = dv >> dec_shift;

        s = D->dec[t];
        n = D->dec[t + 1] + 1;

        /* binary search in [s, n] */
        while (n > s + 1) {
            m = (s + n) >> 1;
            if (D->xs[m].freq > dv) n = m;
            else                    s = m;
        }

        x = D->xs[s].freq * ac->l;
        if (s != DISTSIZE - 1) y = D->xs[s + 1].freq * ac->l;
    }
    else {
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
    }

    ac->v -= x;
    ac->l = y - x;

    if (ac->l < min_length) ac_renormalize_decoder(ac);

    D->xs[s].count++;
    D->use_count++;

    if (--D->update_delay == 0) dfun(update)(D);

    return s;
}



void cdfun(init) (cond_dist_t* D, size_t n, bool decode)
{
    D->xss   = malloc_or_die(n * sizeof(dist_t));
    D->index = malloc_or_die(n * sizeof(uint32_t));
    D->ord   = malloc_or_die(n * sizeof(uint32_t));
    D->n = n;

    D->xss[0].update_delay = DISTSIZE * update_delay_factor;
    if (decode && DISTSIZE > 16) {
        /* allocate enough memory for everyone */
        D->xss[0].dec = malloc_or_die(n * (dec_size + 6) * sizeof(uint16_t));
    }
    else D->xss[0].dec = NULL;

    size_t i;
    for (i = 0; i < DISTSIZE; ++i) D->xss[0].xs[i].count = 1;
    dfun(update)(D->xss);

    D->index[0] = 0;
    D->ord[0]   = 0;

    for (i = 1; i < n; ++i) {
        memcpy(D->xss + i, D->xss, sizeof(dist_t));
        D->xss[i].dec = NULL;
        D->index[i] = i;
        D->ord[i]   = i;
    }

    /* properly initialize the dec tables */
    if (decode && DISTSIZE > 16) {
        for (i = 1; i < n; ++i) {
            D->xss[i].dec = D->xss[i - 1].dec + dec_size + 6;
            memcpy(D->xss[i].dec, D->xss[i - 1].dec, (dec_size + 6) * sizeof(uint16_t));
        }
    }

}


void cdfun(free) (cond_dist_t* D)
{
    if (D == NULL) return;

    free(D->xss[0].dec);
    free(D->xss);
    free(D->index);
    free(D->ord);
}



/* Sort the ord array in a cond_dist. */

#define getkey(i) (D->xss[D->index[D->ord[i]]].use_count)

static void cdfun(inssort_ord) (cond_dist_t* D, size_t i, size_t j)
{
    if (j - i + 1 < 2) return;

    uint32_t key, ord;
    dist_t val;
    size_t u, v;
    for (u = i + 1; u <= j; ++u) {
        key = getkey(u);
        ord = D->ord[u];
        memcpy(&val, &D->ord[u], sizeof(dist_t));

        for (v = u - 1; v < 0U - 1 && key < getkey(v); --v) {
            D->ord[v + 1] = D->ord[v];
            memcpy(&D->xss[v + 1], &D->xss[v], sizeof(dist_t));
        }

        D->ord[v + 1] = ord;
        memcpy(&D->ord[v + 1], &val, sizeof(dist_t));
    }
}


static void cdfun(quicksort_ord) (cond_dist_t* D, size_t i, size_t j)
{
    if (j - i + 1 < 7) {
        cdfun(inssort_ord) (D, i, j);
        return;
    }

    size_t pivot = i + rand() % (j - i + 1);

    uint32_t pivot_key = getkey(pivot);
    uint32_t pivot_ord = D->ord[pivot];
    dist_t   pivot_val;
    memcpy(&pivot_val, &D->xss[pivot], sizeof(dist_t));

    uint32_t key, ord;
    dist_t val;

    uint32_t min_left_key = 0U - 1;

    D->ord[pivot] = D->ord[j];
    D->ord[j]     = pivot_ord;

    memcpy(&D->xss[pivot], &D->xss[j], sizeof(dist_t));
    memcpy(&D->xss[j],     &pivot_val, sizeof(dist_t));

    size_t u, v;
    for (u = i, v = i; v < j; ++v) {
        key = getkey(v);
        if (key <= pivot_key) {

            ord = D->ord[v];
            D->ord[v] = D->ord[u];
            D->ord[u] = ord;

            memcpy(&val,       &D->xss[v], sizeof(dist_t));
            memcpy(&D->xss[v], &D->xss[u], sizeof(dist_t));
            memcpy(&D->xss[u], &val,       sizeof(dist_t));

            if (key < min_left_key) min_left_key = key;

            ++u;
        }
    }

    D->ord[j] = D->ord[u];
    D->ord[u] = pivot_ord;

    memcpy(&D->xss[j], &D->xss[u], sizeof(dist_t));
    memcpy(&D->xss[u], &pivot_val, sizeof(dist_t));

    if (i > u && min_left_key != pivot_key) cdfun(quicksort_ord) (D, i, u - 1);
    cdfun(quicksort_ord) (D, u + 1, j);
}




void cdfun(reorder) (cond_dist_t* D)
{
    /* sort ord based on use count */
    cdfun(quicksort_ord) (D, 0, D->n - 1);

    size_t i;
    for (i = 0; i < D->n; ++i) {
        D->index[D->ord[i]] = i;
    }
}



#undef getkey
#undef getval






void cdfun(encode)(ac_t* ac, cond_dist_t* D, uint32_t y, symb_t x)
{
    /* XXX */
    dfun(encode)(ac, D->xss + D->index[y], x);
    /*dfun(encode)(ac, D->xss + y, x);*/
}


symb_t cdfun(decode)(ac_t* ac, cond_dist_t* D, uint32_t y)
{
    /* XXX */
    return dfun(decode)(ac, D->xss + D->index[y]);
    /*return dfun(decode)(ac, D->xss + y);*/
}


