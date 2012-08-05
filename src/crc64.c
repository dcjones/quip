
#include "crc64.h"
#include "config.h"


/* 64-bit cyclic redundancy check. This was adapted
 * from the public domain implementation by
 * Lasse Collin in liblzma.
 */


#ifdef WORDS_BIGENDIAN
#include "crc64_table_be.h"
#else
#include "crc64_table_le.h"
#endif


#ifndef bswap16
#   define bswap16(num) \
        (((uint16_t)(num) << 8) | ((uint16_t)(num) >> 8))
#endif

#ifndef bswap32
#   define bswap32(num) \
        ( (((uint32_t)(num) << 24)                       ) \
        | (((uint32_t)(num) <<  8) & UINT32_C(0x00FF0000)) \
        | (((uint32_t)(num) >>  8) & UINT32_C(0x0000FF00)) \
        | (((uint32_t)(num) >> 24)                       ) )
#endif

#ifndef bswap64
#   define bswap64(num) \
        ( (((uint64_t)(num) << 56)                               ) \
        | (((uint64_t)(num) << 40) & UINT64_C(0x00FF000000000000)) \
        | (((uint64_t)(num) << 24) & UINT64_C(0x0000FF0000000000)) \
        | (((uint64_t)(num) <<  8) & UINT64_C(0x000000FF00000000)) \
        | (((uint64_t)(num) >>  8) & UINT64_C(0x00000000FF000000)) \
        | (((uint64_t)(num) >> 24) & UINT64_C(0x0000000000FF0000)) \
        | (((uint64_t)(num) >> 40) & UINT64_C(0x000000000000FF00)) \
        | (((uint64_t)(num) >> 56)                               ) )
#endif


#ifdef WORDS_BIGENDIAN
#   define A(x) ((x) >> 24)
#   define B(x) (((x) >> 16) & 0xFF)
#   define C(x) (((x) >> 8) & 0xFF)
#   define D(x) ((x) & 0xFF)

#   define S8(x) ((x) << 8)
#   define S32(x) ((x) << 32)

#   define A1(x) ((x) >> 56)

#else
#   define A(x) ((x) & 0xFF)
#   define B(x) (((x) >> 8) & 0xFF)
#   define C(x) (((x) >> 16) & 0xFF)
#   define D(x) ((x) >> 24)

#   define S8(x) ((x) >> 8)
#   define S32(x) ((x) >> 32)

#   define A1 A
#endif


uint64_t crc64_update(uint8_t* buf, size_t size, uint64_t crc)
{
    crc = ~crc;

#ifdef WORDS_BIGENDIAN
    crc = bswap64(crc);
#endif

    if (size > 4) {
        while ((uintptr_t)(buf) & 3) {
            crc = crc64_table[0][*buf++ ^ A1(crc)] ^ S8(crc);
            --size;
        }

        const uint8_t *const limit = buf + (size & ~(size_t)(3));
        size &= (size_t)(3);

        while (buf < limit) {
#ifdef WORDS_BIGENDIAN
            const uint32_t tmp = (crc >> 32)
                    ^ *(const uint32_t *)(buf);
#else
            const uint32_t tmp = crc ^ *(const uint32_t *)(buf);
#endif
            buf += 4;

            crc = crc64_table[3][A(tmp)]
                ^ crc64_table[2][B(tmp)]
                ^ S32(crc)
                ^ crc64_table[1][C(tmp)]
                ^ crc64_table[0][D(tmp)];
        }
    }

    while (size-- != 0)
        crc = crc64_table[0][*buf++ ^ A1(crc)] ^ S8(crc);

#ifdef WORDS_BIGENDIAN
    crc = bswap64(crc);
#endif

    return ~crc;

    return 0;
}

