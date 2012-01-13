
SUBDIRS = lzma

bin_PROGRAMS=quip

quip_SOURCES = quip.c version.h \
			   ac.h          ac.c \
			   assembler.h   assembler.c \
			   bloom.h       bloom.c \
			   cumdist.h     cumdist.c \
			   idenc.h       idenc.c \
			   kmer.h        kmer.c \
			   kmerhash.h    kmerhash.c \
			   misc.h        misc.c \
			   seqenc.h      seqenc.c \
			   seqset.h      seqset.c \
			   sw.h          sw.c \
			   twobit.h      twobit.c \
			   parse.h       parse.c \
			   qualenc.h     qualenc.c \
			   quiplib.c     quip.h

quip_LDADD = lzma/liblzma.la

