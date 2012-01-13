
#include "idenc.h"
#include "misc.h"

/*
 * TODO
 * All of these functions are placeholders for now.
 */

struct idenc_t_
{

};


idenc_t* idenc_alloc(quip_block_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    return E;
}

void idenc_free(idenc_t* E)
{
    free(E);
}

void idenc_encode(idenc_t* E, const seq_t* s)
{
}

void idenc_clear(idenc_t* E)
{
}




