
AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libquip.a

libquip_a_SOURCES = \
          version.h \
          ac.h              ac.c \
          assembler.h       assembler.c \
          bloom.h           bloom.c \
          dist.h            dist.c \
          crc64.h           crc64.c \
          idenc.h           idenc.c \
          kmer.h            kmer.c \
          kmerhash.h        kmerhash.c \
          misc.h            misc.c \
          samopt.h          samopt.c \
		  samoptenc.h \
          seqenc.h          seqenc.c \
          seqmap.h          seqmap.c \
          strmap.h          strmap.c \
          twobit.h          twobit.c \
          qualenc.h         qualenc.c \
          quiplib.c         quip.h \
          fastqfmt.h        fastqfmt.c \
          samfmt.h          samfmt.c \
          quipfmt.h         quipfmt.c \
          sam/bam.h         sam/bam.c \
          sam/sam.h         sam/sam.c \
          sam/bgzf.h        sam/bgzf.c \
          sam/sam_header.h  sam/sam_header.c \
          sam/kstring.h     sam/kstring.c \
          sam/bam_endian.h \
          sam/bam_import.c \
          sam/bam_aux.c \
          sam/khash.h \
          sam/kseq.h

bin_PROGRAMS = quip

quip_SOURCES = quip.c

quip_LDADD = libquip.a

include_HEADERS = quip.h

noinst_HEADERS = dist_impl.h dist_impl.c \
                 dist_template_on.h dist_template_off.h \
                 crc64_table_be.h crc64_table_le.h

quiplinks = unquip quipcat

install-exec-hook :
	cd $(DESTDIR)$(bindir) && \
	target=`echo quip | sed '$(transform)'`$(EXEEXT) && \
	for name in $(quiplinks); do \
		link=`echo $$name | sed '$(transform)'`$(LN_EXEEXT) && \
		rm -f $$link && \
		$(LN_S) $$target $$link; \
	done

EXTRADIST = sam/COPYING

