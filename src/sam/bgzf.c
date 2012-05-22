/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/*
  2012-03-24 by dcjones: use quip streams for more general I/O.
  2009-06-29 by lh3: cache recent uncompressed blocks.
  2009-06-25 by lh3: optionally use my knetfile library to access file on a FTP.
  2009-06-12 by lh3: support a mode string like "wu" where 'u' for uncompressed output */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "bgzf.h"
#include "../misc.h"

#if defined(_WIN32) || defined(_MSC_VER)
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif

typedef int8_t bgzf_byte_t;

static const int DEFAULT_COMPRESS_LEVEL = 5;
static const int DEFAULT_BLOCK_SIZE = 64 * 1024;
static const int MAX_BLOCK_SIZE = 64 * 1024;

static const int BLOCK_HEADER_LENGTH = 18;
static const int BLOCK_FOOTER_LENGTH = 8;

static const int GZIP_ID1 = 31;
static const int GZIP_ID2 = 139;
static const int CM_DEFLATE = 8;
static const int FLG_FEXTRA = 4;
static const int OS_UNKNOWN = 255;
static const int BGZF_ID1 = 66; // 'B'
static const int BGZF_ID2 = 67; // 'C'
static const int BGZF_LEN = 2;
static const int BGZF_XLEN = 6; // BGZF_LEN+4

static const int GZIP_WINDOW_BITS = -15; // no zlib header
static const int Z_DEFAULT_MEM_LEVEL = 8;


void
packInt16(uint8_t* buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

int
unpackInt16(const uint8_t* buffer)
{
    return (buffer[0] | (buffer[1] << 8));
}

void
packInt32(uint8_t* buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

static inline
int
bgzf_min(int x, int y)
{
    return (x < y) ? x : y;
}

static
void
report_error(BGZF* fp, const char* message) {
    fp->error = message;
}


BGZF* bgzf_open_out(quip_writer_t writer, void* writer_data)
{
    const int compress_level = -1;

    BGZF* fp = malloc(sizeof(BGZF));

    fp->file_descriptor = -1;
    fp->open_mode = 'w';
    fp->owned_file = 0;
    fp->compress_level = compress_level < 0? DEFAULT_COMPRESS_LEVEL: compress_level; // Z_DEFAULT_COMPRESSION==-1
    if (fp->compress_level > 9) fp->compress_level = DEFAULT_COMPRESS_LEVEL;

    fp->writer = writer;
    fp->writer_data = writer_data;

    fp->uncompressed_block_size = DEFAULT_BLOCK_SIZE;
    fp->uncompressed_block = NULL;
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
    fp->block_offset = 0;
    fp->block_length = 0;

    return fp;

}


BGZF* bgzf_open_in (quip_reader_t reader, void* reader_data)
{
    BGZF* fp = malloc(sizeof(BGZF));

    fp->uncompressed_block_size = MAX_BLOCK_SIZE;
    fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
    fp->file_descriptor = -1;
    fp->open_mode = 'r';

    fp->reader      = reader;
    fp->reader_data = reader_data;
    
    fp->block_offset = 0;
    fp->block_length = 0;

    return fp;
}


static
int
deflate_block(BGZF* fp, int block_length)
{
    // Deflate the block in fp->uncompressed_block into fp->compressed_block.
    // Also adds an extra field that stores the compressed block length.

    bgzf_byte_t* buffer = fp->compressed_block;
    int buffer_size = fp->compressed_block_size;

    // Init gzip header
    buffer[0] = GZIP_ID1;
    buffer[1] = GZIP_ID2;
    buffer[2] = CM_DEFLATE;
    buffer[3] = FLG_FEXTRA;
    buffer[4] = 0; // mtime
    buffer[5] = 0;
    buffer[6] = 0;
    buffer[7] = 0;
    buffer[8] = 0;
    buffer[9] = OS_UNKNOWN;
    buffer[10] = BGZF_XLEN;
    buffer[11] = 0;
    buffer[12] = BGZF_ID1;
    buffer[13] = BGZF_ID2;
    buffer[14] = BGZF_LEN;
    buffer[15] = 0;
    buffer[16] = 0; // placeholder for block length
    buffer[17] = 0;

    // loop to retry for blocks that do not compress enough
    int input_length = block_length;
    int compressed_length = 0;
    while (1) {
        z_stream zs;
        zs.zalloc = NULL;
        zs.zfree = NULL;
        zs.next_in = fp->uncompressed_block;
        zs.avail_in = input_length;
        zs.next_out = (void*)&buffer[BLOCK_HEADER_LENGTH];
        zs.avail_out = buffer_size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

        int status = deflateInit2(&zs, fp->compress_level, Z_DEFLATED,
                                  GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
        if (status != Z_OK) {
            report_error(fp, "deflate init failed");
            return -1;
        }
        status = deflate(&zs, Z_FINISH);
        if (status != Z_STREAM_END) {
            deflateEnd(&zs);
            if (status == Z_OK) {
                // Not enough space in buffer.
                // Can happen in the rare case the input doesn't compress enough.
                // Reduce the amount of input until it fits.
                input_length -= 1024;
                if (input_length <= 0) {
                    // should never happen
                    report_error(fp, "input reduction failed");
                    return -1;
                }
                continue;
            }
            report_error(fp, "deflate failed");
            return -1;
        }
        status = deflateEnd(&zs);
        if (status != Z_OK) {
            report_error(fp, "deflate end failed");
            return -1;
        }
        compressed_length = zs.total_out;
        compressed_length += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
        if (compressed_length > MAX_BLOCK_SIZE) {
            // should never happen
            report_error(fp, "deflate overflow");
            return -1;
        }
        break;
    }

    packInt16((uint8_t*)&buffer[16], compressed_length-1);
    uint32_t crc = crc32(0L, NULL, 0L);
    crc = crc32(crc, fp->uncompressed_block, input_length);
    packInt32((uint8_t*)&buffer[compressed_length-8], crc);
    packInt32((uint8_t*)&buffer[compressed_length-4], input_length);

    int remaining = block_length - input_length;
    if (remaining > 0) {
        if (remaining > input_length) {
            // should never happen (check so we can use memcpy)
            report_error(fp, "remainder too large");
            return -1;
        }
        memcpy((uint8_t*) fp->uncompressed_block,
               (uint8_t*) fp->uncompressed_block + input_length,
               remaining);
    }
    fp->block_offset = remaining;
    return compressed_length;
}

static
int
inflate_block(BGZF* fp, int block_length)
{
    // Inflate the block in fp->compressed_block into fp->uncompressed_block

    z_stream zs;
	int status;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = (uint8_t*) fp->compressed_block + 18;
    zs.avail_in = block_length - 16;
    zs.next_out = fp->uncompressed_block;
    zs.avail_out = fp->uncompressed_block_size;

    status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if (status != Z_OK) {
        report_error(fp, "inflate init failed");
        return -1;
    }
    status = inflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END) {
        inflateEnd(&zs);
        report_error(fp, "inflate failed");
        return -1;
    }
    status = inflateEnd(&zs);
    if (status != Z_OK) {
        report_error(fp, "inflate failed");
        return -1;
    }
    return zs.total_out;
}

static
int
check_header(const bgzf_byte_t* header)
{
    return (header[0] == GZIP_ID1 &&
            header[1] == (bgzf_byte_t) GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & FLG_FEXTRA) != 0 &&
            unpackInt16((uint8_t*)&header[10]) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            unpackInt16((uint8_t*)&header[14]) == BGZF_LEN);
}


int
bgzf_read_block(BGZF* fp)
{
    bgzf_byte_t header[BLOCK_HEADER_LENGTH];
	int count, size = 0, block_length, remaining;

    count = fp->reader(fp->reader_data, (uint8_t*) header, sizeof(header));

    if (count == 0) {
        fp->block_length = 0;
        return 0;
    }
	size = count;
    if ((size_t) count != sizeof(header)) {
        report_error(fp, "read failed");
        return -1;
    }
    if (!check_header(header)) {
        report_error(fp, "invalid block header");
        return -1;
    }
    block_length = unpackInt16((uint8_t*)&header[16]) + 1;
    bgzf_byte_t* compressed_block = (bgzf_byte_t*) fp->compressed_block;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;

    count = fp->reader(fp->reader_data, (uint8_t*) &compressed_block[BLOCK_HEADER_LENGTH], remaining);

    if (count != remaining) {
        report_error(fp, "read failed");
        return -1;
    }
	size += count;
    count = inflate_block(fp, block_length);
    if (count < 0) return -1;
    if (fp->block_length != 0) {
        // Do not reset offset if this read follows a seek.
        fp->block_offset = 0;
    }
    fp->block_length = count;
    return 0;
}

int
bgzf_read(BGZF* fp, void* data, int length)
{
    if (length <= 0) {
        return 0;
    }
    if (fp->open_mode != 'r') {
        report_error(fp, "file not open for reading");
        return -1;
    }

    int bytes_read = 0;
    bgzf_byte_t* output = data;
    while (bytes_read < length) {
        int copy_length, available = fp->block_length - fp->block_offset;
		bgzf_byte_t *buffer;
        if (available <= 0) {
            if (bgzf_read_block(fp) != 0) {
                return -1;
            }
            available = fp->block_length - fp->block_offset;
            if (available <= 0) {
                break;
            }
        }
        copy_length = bgzf_min(length-bytes_read, available);
        buffer = fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;
    }
    if (fp->block_offset == fp->block_length) {
        fp->block_offset = 0;
        fp->block_length = 0;
    }
    return bytes_read;
}

int bgzf_flush(BGZF* fp)
{
    while (fp->block_offset > 0) {
        int block_length;
		block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) return -1;
        fp->writer(fp->writer_data, (uint8_t*) fp->compressed_block, block_length);
    }
    return 0;
}

int bgzf_flush_try(BGZF *fp, int size)
{
	if (fp->block_offset + size > fp->uncompressed_block_size)
		return bgzf_flush(fp);
	return -1;
}

int bgzf_write(BGZF* fp, const void* data, int length)
{
	const bgzf_byte_t *input = data;
	int block_length, bytes_written;
    if (fp->open_mode != 'w') {
        report_error(fp, "file not open for writing");
        return -1;
    }

    if (fp->uncompressed_block == NULL)
        fp->uncompressed_block = malloc(fp->uncompressed_block_size);

    input = data;
    block_length = fp->uncompressed_block_size;
    bytes_written = 0;
    while (bytes_written < length) {
        int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
        bgzf_byte_t* buffer = fp->uncompressed_block;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        bytes_written += copy_length;
        if (fp->block_offset == block_length) {
            if (bgzf_flush(fp) != 0) {
                break;
            }
        }
    }
    return bytes_written;
}

int bgzf_close(BGZF* fp)
{
    if (fp->open_mode == 'w') {
        if (bgzf_flush(fp) != 0) return -1;

        int block_length = deflate_block(fp, 0);
        fp->writer(fp->writer_data, (uint8_t*) fp->compressed_block, block_length);
    }

    free(fp->uncompressed_block);
    free(fp->compressed_block);
    free(fp);
    return 0;
}

int bgzf_check_EOF(ATTRIB_UNUSED BGZF *fp)
{
    /* quip I/O behaves like a pipe, so we cannot seek to check
     * if the file as the proper EOF magic. */
    errno = ESPIPE;
    return -1;
}

int64_t bgzf_seek(ATTRIB_UNUSED BGZF* fp, ATTRIB_UNUSED int64_t pos, ATTRIB_UNUSED int where)
{
    /* quip I/O behaves like a pipe, so we cannot seek. */
    errno = ESPIPE;
    return -1;
}
