
#include "idenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>

/* does a particular character constitute a separator */
static const bool issep[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static const size_t max_group_len  = 50;
static const size_t max_offset     = 16;

/* method of encoding a field */
typedef enum {
    ID_GROUP_MATCH = 0,
    ID_GROUP_STR,
    ID_GROUP_NUM,
    ID_GROUP_NUM_OFF
} enc_t;

/* token type */
typedef enum {
    ID_TOK_STR = 0,
    ID_TOK_NUM
} tok_type_t;


/* groups that do match the previous id are encoded either as a numerical offset
 * or as text. */
typedef struct {
    tok_type_t  type;
    uint64_t    num;
    size_t      pos;
    size_t      len;
} tok_t;


/* tokenization of ids */
const uint8_t* next_id_token(const uint8_t* id, tok_type_t* t)
{
    const uint8_t* s = id;

    if (*s == '\0') {
        ++s;
        *t = ID_TOK_STR;
    }
    else if (issep[*s]) {
        ++s;
        while ((size_t) (s - id) < max_group_len && issep[*s]) ++s;

        if ((size_t) (s - id) < max_group_len && !('1' <= *s && *s <= '9')) {
            while ((size_t) (s - id) < max_group_len &&
                   *s != '\0' && !issep[*s]) ++s;
        }

        *t = ID_TOK_STR;
    }
    else if ('1' <= *s && *s <= '9') {
        ++s;
        while ((size_t) (s - id) < max_group_len && '0' <= *s && *s <= '9') ++s;
        *t = ID_TOK_NUM;
    }
    else {
        ++s;
        while ((size_t) (s - id) < max_group_len &&
               *s != '\0' && !issep[*s]) ++s;
        *t = ID_TOK_STR;
    }

    return s;
}



struct idenc_t_
{
    ac_t* ac;

    /* distribution over how the group is stored (numeric offset or text) */
    dist4_t* d_type;

    /* distribution over number of leading matches given group */
    dist50_t* d_str_matches;

    /* distribution over the length text given group */
    dist50_t* d_str_len;

    /* distribution over new characters given group */
    cond_dist128_t* d_str_char;

    /* distribution over numerical offset given group */
    dist16_t* d_off;

    /* distribution over byte values used when encoding numbers */
    cond_dist256_t* d_num;

    /* number of tokenized groups allowed for */
    size_t max_group_cnt;

    /* previous id string */
    uint8_t* lastid;
    size_t lastid_len, lastid_size;

    /* previous id token stream */
    tok_t* toks;
    size_t toks_len, toks_size;
};


static void idenc_init(idenc_t* E)
{
    E->max_group_cnt = 0;

    E->d_type        = NULL;
    E->d_str_matches = NULL;
    E->d_str_len     = NULL;
    E->d_str_char    = NULL;
    E->d_off         = NULL;
    E->d_num         = NULL;

    E->lastid = NULL;
    E->lastid_len  = 0;
    E->lastid_size = 0;

    E->toks = NULL;
    E->toks_len  = 0;
    E->toks_size = 0;
}



idenc_t* idenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    idenc_init(E);

    return E;
}



idenc_t* idenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_decoder(reader, reader_data);

    idenc_init(E);

    return E;
}



void idenc_free(idenc_t* E)
{
    free(E->lastid);
    free(E->toks);
    ac_free(E->ac);

    size_t i;
    for (i = 0; i < E->max_group_cnt; ++i){
        cond_dist128_free(&E->d_str_char[i]);
        cond_dist256_free(&E->d_num[i]);
    }

    free(E->d_type);
    free(E->d_str_matches);
    free(E->d_str_len);
    free(E->d_str_char);
    free(E->d_off);
    free(E->d_num);

    free(E);
}


static void idenc_add_group(idenc_t* E)
{
    E->max_group_cnt++;

    E->d_type        = realloc_or_die(E->d_type,        E->max_group_cnt * sizeof(dist4_t));
    E->d_str_matches = realloc_or_die(E->d_str_matches, E->max_group_cnt * sizeof(dist50_t));
    E->d_str_len     = realloc_or_die(E->d_str_len,     E->max_group_cnt * sizeof(dist50_t));
    E->d_str_char    = realloc_or_die(E->d_str_char,    E->max_group_cnt * sizeof(cond_dist128_t));
    E->d_off         = realloc_or_die(E->d_off,         E->max_group_cnt * sizeof(dist16_t));
    E->d_num         = realloc_or_die(E->d_num,         E->max_group_cnt * sizeof(cond_dist256_t));

    size_t i = E->max_group_cnt - 1;
    dist4_init  (&E->d_type[i]);
    dist50_init (&E->d_str_matches[i]);
    dist50_init (&E->d_str_len[i]);
    dist16_init (&E->d_off[i]);

    /* these allocate quite a bit, so we initialized them lazily */
    memset(&E->d_str_char[i], 0, sizeof(cond_dist128_t));
    memset(&E->d_num[i], 0, sizeof(cond_dist256_t));
}


static void encode_str(idenc_t* E, size_t i, const uint8_t* str, tok_t* tok)
{
    /* previous id string */
    size_t prev_tok_len = E->toks_len > i ? E->toks[i].len : 0;
    const uint8_t* t = E->toks_len > i ? E->lastid + E->toks[i].pos : NULL;
    const uint8_t* s = str + tok->pos;

    size_t j = 0; /* string offset */
    while (j < tok->len && j < prev_tok_len && s[j] == t[j]) ++j;

    /* lazy initialization of d_str_char */
    if (E->d_str_char[i].n == 0) {
        cond_dist128_init(&E->d_str_char[i], 128);
    }

    if (tok->len == prev_tok_len && j == tok->len) {
        dist4_encode(E->ac, &E->d_type[i], ID_GROUP_MATCH);
    }
    else {
        dist4_encode  (E->ac, &E->d_type[i], ID_GROUP_STR);
        dist50_encode (E->ac, &E->d_str_matches[i], j);
        dist50_encode (E->ac, &E->d_str_len[i], tok->len - j);

        char c = 0;
        for (; j < tok->len; ++j) {
            cond_dist128_encode(E->ac, &E->d_str_char[i], c, s[j]);
            c = s[j];
        }
    }
}


static void encode_num(idenc_t* E, size_t i, const uint8_t* str, tok_t* tok)
{
    const uint8_t* s = str + tok->pos;
    char* last;

    /* is this the same the number that was previously encoded:
        avoid relatively expensive to strtoull when we can */
    if (i < E->toks_len && E->toks[i].type == ID_TOK_NUM &&
        tok->len == E->toks[i].len && memcmp(str, E->lastid + E->toks[i].pos, tok->len) == 0) {

        dist4_encode(E->ac, &E->d_type[i], ID_GROUP_NUM_OFF);
        dist16_encode(E->ac, &E->d_off[i], 0);
        tok->num = E->toks[i].num;
        return;
    }

    uint64_t off, x = strtoull((char*) s, &last, 10);

    /* if conversion to 64-bit unsigned int was not possible,
       encode the number as as string. */
    if (x == ULLONG_MAX || last != (char*) s + tok->len) {
        tok->type = ID_TOK_STR;
        encode_str(E, i, str, tok);
        return;
    }

    tok->num = x;

    /* lazy initialization of d_num */
    if (E->d_num[i].n == 0) {
        cond_dist256_init(&E->d_num[i], 9 * 256);
    }

    /* Encode numbers as offsets from a previous number when possible. */
    if (i < E->toks_len && E->toks[i].type == ID_TOK_NUM &&
        x >= E->toks[i].num && (off = x - E->toks[i].num) < max_offset) {

        dist4_encode(E->ac, &E->d_type[i], ID_GROUP_NUM_OFF);
        dist16_encode(E->ac, &E->d_off[i], off);

    }
    else {
        dist4_encode(E->ac, &E->d_type[i], ID_GROUP_NUM);
        dist_encode_uint64(E->ac, &E->d_num[i], x);
    }
}



void idenc_encode(idenc_t* E, const short_read_t* seq)
{
    size_t i; /* token number */
    const uint8_t *s, *t;
    t = seq->id.s;
    tok_t tok;

    for (i = 0, s = seq->id.s; true; ++i, s = t) {
        if (i + 1 > E->max_group_cnt) idenc_add_group(E);

        t = next_id_token(s, &tok.type);
        tok.pos = s - seq->id.s;
        tok.len = t - s;

        if      (tok.type == ID_TOK_STR) encode_str(E, i, seq->id.s, &tok);
        else if (tok.type == ID_TOK_NUM) encode_num(E, i, seq->id.s, &tok);

        if (E->toks_size <= i) {
            E->toks_size += 1;
            E->toks = realloc_or_die(E->toks, E->toks_size * sizeof(tok_t));
        }
        memcpy(E->toks + i, &tok, sizeof(tok_t));

        if (*s == '\0') break;
    }


    if (E->lastid_size < seq->id.n + 1) {
        E->lastid_size = seq->id.n + 1;
        free(E->lastid);
        E->lastid = malloc_or_die((seq->id.n + 1) * sizeof(uint8_t));
    }
    memcpy(E->lastid, seq->id.s, (seq->id.n + 1) * sizeof(uint8_t));
    E->lastid_len = seq->id.n + 1;

    E->toks_len = i - 1;
}


void idenc_decode(idenc_t* E, short_read_t* seq)
{
    str_t* id = &seq->id; /* for convenience */
    size_t i; /* token number */
    size_t j; /* offset into id */
    uint32_t matches, mismatches;
    enc_t type;
    tok_t tok;
    uint64_t off;

    for (i = 0, j = 0; true; ++i) {
        if (i + 1 > E->max_group_cnt) idenc_add_group(E);

        type = dist4_decode(E->ac, &E->d_type[i]);

        if (type == ID_GROUP_MATCH) {
            str_reserve(id, j + E->toks[i].len + 1);
            memcpy(id->s + j, E->lastid + E->toks[i].pos, E->toks[i].len);

            tok.type = ID_TOK_STR;
            tok.pos  = j;
            tok.len  = E->toks[i].len;
            j += tok.len;
        }
        else if (type == ID_GROUP_STR) {
            matches    = dist50_decode(E->ac, &E->d_str_matches[i]);
            mismatches = dist50_decode(E->ac, &E->d_str_len[i]);

            tok.type = ID_TOK_STR;
            tok.pos  = j;
            tok.len  = matches + mismatches;;

            str_reserve(id, j + matches + mismatches + 1);
            if (matches > 0) {
                memcpy(id->s + j, E->lastid + E->toks[i].pos, matches);
                j += matches;
            }

            /* lazy initialization of d_str_char */
            if (E->d_str_char[i].n == 0) {
                cond_dist128_init(&E->d_str_char[i], 128);
            }

            char c = 0;
            while (mismatches--) {
                c = id->s[j++] = cond_dist128_decode(E->ac, &E->d_str_char[i], c);
            }
        }
        else if (type == ID_GROUP_NUM_OFF) {
            str_reserve(id, j + 20);

            off = dist16_decode(E->ac, &E->d_off[i]);

            tok.type = ID_TOK_NUM;
            tok.pos  = j;
            tok.num  = E->toks[i].num + off;
            tok.len  = snprintf((char*) id->s + j, 20, "%"PRIu64, tok.num);
            j += tok.len;
        }
        else if (type == ID_GROUP_NUM) {
            str_reserve(id, j + 20);

            /* lazy initialization of d_num */
            if (E->d_num[i].n == 0) {
                cond_dist256_init(&E->d_num[i], 9 * 256);
            }

            tok.type = ID_TOK_NUM;
            tok.pos  = j;
            tok.num  = dist_decode_uint64(E->ac, &E->d_num[i]);
            tok.len  = snprintf((char*) id->s + j, 20, "%"PRIu64, tok.num);
            j += tok.len;
        }

        if (E->toks_size <= i) {
            E->toks_size += 1;
            E->toks = realloc_or_die(E->toks, E->toks_size * sizeof(tok_t));
        }
        memcpy(E->toks + i, &tok, sizeof(tok_t));

        if (id->s[tok.pos] == '\0') break;
    }

    seq->id.n = j - 1;

    if (E->lastid_size < seq->id.n + 1) {
        E->lastid_size = seq->id.n + 1;
        free(E->lastid);
        E->lastid = malloc_or_die((seq->id.n + 1) * sizeof(char));
    }
    memcpy(E->lastid, seq->id.s, (seq->id.n + 1) * sizeof(char));
    E->lastid_len = seq->id.n + 1;

    E->toks_len = i - 1;
}


size_t idenc_finish(idenc_t* E)
{
    return ac_finish_encoder(E->ac);
}


void idenc_flush(idenc_t* E)
{
    ac_flush_encoder(E->ac);
}



void idenc_start_decoder(idenc_t* E)
{
    ac_start_decoder(E->ac);
}


void idenc_reset_decoder(idenc_t* E)
{
    ac_reset_decoder(E->ac);
}


