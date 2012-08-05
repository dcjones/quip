#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include "sam.h"

void bam_init_header_hash(bam_header_t *header);

#define TYPE_BAM  1
#define TYPE_READ 2

bam_header_t *bam_header_dup(const bam_header_t *h0)
{
	bam_header_t *h;
	int i;
	h = bam_header_init();
	*h = *h0;
	h->hash = h->dict = h->rg2lib = 0;
	h->text = (char*)calloc(h->l_text + 1, 1);
	memcpy(h->text, h0->text, h->l_text);
	h->target_len = (uint32_t*)calloc(h->n_targets, 4);
	h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
	for (i = 0; i < h->n_targets; ++i) {
		h->target_len[i] = h0->target_len[i];
		h->target_name[i] = strdup(h0->target_name[i]);
	}

	bam_init_header_hash(h);

	return h;
}
static void append_header_text(bam_header_t *header, char* text, int len)
{
	int x = header->l_text + 1;
	int y = header->l_text + len + 1; // 1 byte null
	if (text == 0) return;
	kroundup32(x);
	kroundup32(y);
	if (x < y) header->text = (char*)realloc(header->text, y);
	strncpy(header->text + header->l_text, text, len); // we cannot use strcpy() here.
	header->l_text += len;
	header->text[header->l_text] = 0;
}


samfile_t* samopen_out(quip_writer_t writer, void* writer_data, bool binary, void* aux)
{
	samfile_t *fp;
	fp = (samfile_t*)calloc(1, sizeof(samfile_t));

	fp->header = bam_header_dup((const bam_header_t*)aux);

	if (binary) { // binary
		fp->type |= TYPE_BAM;
		fp->x.bam = bam_open_out(writer, writer_data);
		if (fp->x.bam == 0) goto open_err_ret;
		bam_header_write(fp->x.bam, fp->header);

	} else { // text
		// open file

		fp->x.tamw.writer      = writer;
		fp->x.tamw.writer_data = writer_data;
		fp->type |= BAM_OFDEC<<2;

		// write header
		if (aux) {
			int i;
			bam_header_t *alt;
			// parse the header text
			alt = bam_header_init();
			alt->l_text = fp->header->l_text; alt->text = fp->header->text;
			sam_header_parse(alt);
			alt->l_text = 0; alt->text = 0;
			// check if there are @SQ lines in the header
			writer(writer_data, (uint8_t*) fp->header->text, fp->header->l_text);

			if (alt->n_targets) { // then write the header text without dumping ->target_{name,len}
				if (alt->n_targets != fp->header->n_targets && bam_verbose >= 1)
					fprintf(stderr, "[samopen] inconsistent number of target sequences. Output the text header.\n");
			} else { // then dump ->target_{name,len}
				char target_len_str[11];
				for (i = 0; i < fp->header->n_targets; ++i) {
					writer(writer_data, (uint8_t*) "@SQ\tSN:", 7);
					writer(writer_data, (uint8_t*) fp->header->target_name[i], strlen(fp->header->target_name[i]));
					writer(writer_data, (uint8_t*) "\tLN:", 4);
					snprintf(target_len_str, sizeof(target_len_str), "%"PRIu32, fp->header->target_len[i]);
					writer(writer_data, (uint8_t*) target_len_str, strlen(target_len_str));
					writer(writer_data, (uint8_t*) "\n", 1);
				}
			}
			bam_header_destroy(alt);
		}
	}

	return fp;

open_err_ret:
	free(fp);
	return 0;
}



samfile_t* samopen_in(quip_reader_t reader, void* reader_data, bool binary, void* aux)
{
	samfile_t *fp;
	fp = (samfile_t*)calloc(1, sizeof(samfile_t));

	fp->type |= TYPE_READ;
	if (binary) { // binary
		fp->type |= TYPE_BAM;
		fp->x.bam = bam_open_in(reader, reader_data);
		if (fp->x.bam == 0) goto open_err_ret;
		fp->header = bam_header_read(fp->x.bam);
	} else { // text
		fp->x.tamr = sam_open_in(reader, reader_data);
		if (fp->x.tamr == 0) goto open_err_ret;
		fp->header = sam_header_read(fp->x.tamr);
		if (fp->header->n_targets == 0) { // no @SQ fields
			if (aux) { // check if aux is present
				bam_header_t *textheader = fp->header;
				fp->header = sam_header_read2((const char*)aux);
				if (fp->header == 0) goto open_err_ret;
				append_header_text(fp->header, textheader->text, textheader->l_text);
				bam_header_destroy(textheader);
			}
			if (fp->header->n_targets == 0 && bam_verbose >= 1)
				fprintf(stderr, "[samopen] no @SQ lines in the header.\n");
		} else if (bam_verbose >= 2) fprintf(stderr, "[samopen] SAM header is present: %d sequences.\n", fp->header->n_targets);
	}

	return fp;

open_err_ret:
	free(fp);
	return 0;
}


void samclose(samfile_t *fp)
{
	if (fp == 0) return;
	if (fp->header) bam_header_destroy(fp->header);
	if (fp->type & 1) bam_close(fp->x.bam);
	else if (fp->type == 2) sam_close(fp->x.tamr);
	free(fp);
}

int samread(samfile_t *fp, bam1_t *b)
{
	if (fp == 0 || !(fp->type & TYPE_READ)) return -1; // not open for reading
	if (fp->type & TYPE_BAM) return bam_read1(fp->x.bam, b);
	else return sam_read1(fp->x.tamr, fp->header, b);
}

int samwrite(samfile_t *fp, const bam1_t *b)
{
	if (fp == 0 || (fp->type & TYPE_READ)) return -1; // not open for writing
	if (fp->type & TYPE_BAM) return bam_write1(fp->x.bam, b);
	else {
		char *s = bam_format1_core(fp->header, b, fp->type>>2&3);
		int l = strlen(s);
		fp->x.tamw.writer(fp->x.tamw.writer_data, (uint8_t*) s, strlen(s));
		fp->x.tamw.writer(fp->x.tamw.writer_data, (uint8_t*) "\n", 1);
		free(s);
		return l + 1;
	}
}

