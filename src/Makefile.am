

lib_LTLIBRARIES = libquip.la

libquip_la_SOURCES = \
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
          quip_parse.h    quip_parse.c \
          qualenc.h       qualenc.c \
          qualenc_prior.h qualenc_prior.c \
          quiplib.c       quip.h

bin_PROGRAMS = quip

quip_SOURCES = quip.c

quip_LDADD = libquip.la

include_HEADERS = quip.h quip_parse.h

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


