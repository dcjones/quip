
#include "idenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>

/* does a particular character constitute a separator */
static const bool issep[256] =
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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
    ID_GROUP_OFF,
    ID_GROUP_NUM
} enc_t;

/* token type */
typedef enum {
    ID_TOK_END = 0,
    ID_TOK_STR,
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
const char* next_id_token(const char* id, tok_type_t* t)
{
    const char* s = id;

    if (*s == '\0') {
        ++s;
        *t = ID_TOK_END;
    }
    else if (issep[(uint8_t) *s]) {
        ++s;
        while ((size_t) (s - id) < max_group_len && issep[(uint8_t) *s]) ++s;
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
               *s != '\0' && !issep[(uint8_t) *s]) ++s;
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
    char* lastid;
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


static void encode_str(idenc_t* E, size_t i, const char* str, tok_t* tok)
{
    /* previous id string */
    const char* t = E->toks_len > i ? E->lastid + E->toks[i].pos : "\0";
    const char* s = str + tok->pos;

    size_t j = 0; /* string offset */
    while (j < tok->len && s[j] != '\0' && t[j] != '\0' && s[j] == t[j]) ++j;

    /* lazy initialization of d_str_char */
    if (E->d_str_char[i].n == 0) {
        cond_dist128_init(&E->d_str_char[i], 128);
    }

    if (j == tok->len) {
        dist4_encode(E->ac, &E->d_type[i], ID_GROUP_MATCH);
    }
    else {
        dist4_encode  (E->ac, &E->d_type[i], ID_GROUP_STR);
        dist50_encode (E->ac, &E->d_str_matches[i], j);
        dist50_encode (E->ac, &E->d_str_len[i], tok->len - j);

        char c = 0;
        for (; j < tok->len; ++j) {
            cond_dist128_encode(E->ac, &E->d_str_char[i], c, s[i]);
            c = s[i];
        }
    }
}


static void encode_num(idenc_t* E, size_t i, const  char* str, tok_t* tok)
{
    const char* s = str + tok->pos;
    char* last;
    uint64_t off, x = strtoull(s, &last, 10);


    /* if conversion to 64-bit unsigned int was not possible,
       encode the number as as string. */
    if (x == ULLONG_MAX || last != s + tok->len) {
        tok->type = ID_TOK_STR;
        encode_str(E, i, str, tok);
    }

    tok->num = x;

    /* lazy initialization of d_num */
    if (E->d_num[i].n == 0) {
        cond_dist256_init(&E->d_num[i], 8 * 256);
    }

    /* Encode numbers as offsets from a previous number when possible. */
    if (i < E->toks_len && E->toks[i].type == ID_TOK_NUM &&
        x >= E->toks[i].num && (off = x - E->toks[i].num) < max_offset) {

        dist4_encode(E->ac, &E->d_type[i], ID_GROUP_OFF);
        dist16_encode(E->ac, &E->d_off[i], off);

    }
    else {
        dist4_encode(E->ac, &E->d_type[i], ID_GROUP_NUM);
        dist_encode_uint64(E->ac, &E->d_num[i], x);
    }
}



void idenc_encode(idenc_t* E, const seq_t* seq)
{
    size_t i; /* token number */
    const char *s, *t;
    t = seq->id1.s;
    tok_t tok;

    for (i = 0, s = seq->id1.s; true; ++i, s = t) {
        if (i + 1 > E->max_group_cnt) idenc_add_group(E);

        t = next_id_token(s, &tok.type);
        tok.pos = s - seq->id1.s;
        tok.len = t - s;

        if      (tok.type == ID_TOK_STR) encode_str(E, i, seq->id1.s, &tok);
        else if (tok.type == ID_TOK_NUM) encode_num(E, i, seq->id1.s, &tok);
        else    break;

        if (E->toks_size <= i) {
            E->toks_size += 1;
            E->toks = realloc_or_die(E->toks, E->toks_size * sizeof(tok_t));
        }
        memcpy(E->toks + i, &tok, sizeof(tok_t));
    }


    if (E->lastid_size < seq->id1.n + 1) {
        E->lastid_size = seq->id1.n + 1;
        free(E->lastid);
        E->lastid = malloc_or_die((seq->id1.n + 1) * sizeof(char));
    }
    memcpy(E->lastid, seq->id1.s, (seq->id1.n + 1) * sizeof(char));
    E->lastid_len = seq->id1.n + 1;

    E->toks_len = i - 1;
}



void idenc_flush(idenc_t* E)
{
    ac_flush_encoder(E->ac);
}


void idenc_decode(idenc_t* E, seq_t* seq)
{
    str_t* id = &seq->id1; /* for convenience */

    /* TODO */
}


void idenc_start_decoder(idenc_t* E)
{
    ac_start_decoder(E->ac);
}


void idenc_reset_decoder(idenc_t* E)
{
    ac_reset_decoder(E->ac);
}


