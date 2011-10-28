
#include "bloom.h"
#include "bitset.h"


struct bloom_t_
{
    bitset_t* bits;

};


bloom_t* bloom_alloc(size_t m, size_t k)
{
    bloom_t* B = malloc(sizeof(bloom_t));
    B->bits = bitset_alloc(m);

    return B;
}


void bloom_free(bloom_t* B)
{
    bitset_free(B->bits);
    free(B);
}



