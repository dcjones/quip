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

#ifndef __BGZF_H
#define __BGZF_H

#include "../quip.h"
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

//typedef int8_t bool;

typedef struct {
    int file_descriptor;
    char open_mode;  // 'r' or 'w'
    int16_t owned_file, compress_level;

    quip_reader_t reader;
    void*         reader_data;

    quip_writer_t writer;
    void*         writer_data;

    int uncompressed_block_size;
    int compressed_block_size;
    void* uncompressed_block;
    void* compressed_block;
    int block_length;
    int block_offset;
    const char* error;
} BGZF;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Open the specified file for reading or writing.
 * Mode must be either "r" or "w".
 * Returns null on error.
 */
BGZF* bgzf_open_out(quip_writer_t writer, void* writer_data);
BGZF* bgzf_open_in (quip_reader_t reader, void* reader_data);

/*
 * Close the BGZ file and free all associated resources.
 * Does not close the underlying file descriptor if created with bgzf_fdopen.
 * Returns zero on success, -1 on error.
 */
int bgzf_close(BGZF* fp);

/*
 * Read up to length bytes from the file storing into data.
 * Returns the number of bytes actually read.
 * Returns zero on end of file.
 * Returns -1 on error.
 */
int bgzf_read(BGZF* fp, void* data, int length);

/*
 * Write length bytes from data to the file.
 * Returns the number of bytes written.
 * Returns -1 on error.
 */
int bgzf_write(BGZF* fp, const void* data, int length);

int bgzf_check_EOF(BGZF *fp);
int bgzf_read_block(BGZF* fp);
int bgzf_flush(BGZF* fp);
int bgzf_flush_try(BGZF *fp, int size);

#ifdef __cplusplus
}
#endif

static inline int bgzf_getc(BGZF *fp)
{
	int c;
	if (fp->block_offset >= fp->block_length) {
		if (bgzf_read_block(fp) != 0) return -2; /* error */
		if (fp->block_length == 0) return -1; /* end-of-file */
	}
	c = ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    if (fp->block_offset == fp->block_length) {
        fp->block_offset = 0;
        fp->block_length = 0;
    }
	return c;
}

#endif
