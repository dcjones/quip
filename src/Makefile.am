
SUBDIRS = sam

lib_LIBRARIES = libquip.a

libquip_a_SOURCES = \
          version.h \
          ac.h            ac.c \
          assembler.h     assembler.c \
          bloom.h         bloom.c \
          dist.h          dist.c \
          crc64.h         crc64.c \
          idenc.h         idenc.c \
          kmer.h          kmer.c \
          kmerhash.h      kmerhash.c \
          misc.h          misc.c \
          seqenc.h        seqenc.c \
          seqenc_prior.h  seqenc_prior.c \
          seqset.h        seqset.c \
          twobit.h        twobit.c \
          qualenc.h       qualenc.c \
          quiplib.c       quip.h \
          fastqfmt.h      fastqfmt.c \
          samfmt.h        samfmt.c \
          quipfmt.h       quipfmt.c

libquip_a_LIBADD = sam/libsam.a

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


