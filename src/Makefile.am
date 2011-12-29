
bin_PROGRAMS=faze kmerhash_test

faze_SOURCES = faze.c version.h \
			   assembler.h   assembler.c \
			   bloom.h       bloom.c \
			   cumdist.h     cumdist.c \
			   kmer.h        kmer.c \
			   kmerhash.h    kmerhash.c \
			   misc.h        misc.c \
			   seqset.h      seqset.c \
			   sw.h          sw.c \
			   twobit.h      twobit.c \
			   parse.h       parse.c \
			   qual.h        qual.c

faze_LDADD = -lm


kmerhash_test_SOURCES=kmerhash_test.c \
					  misc.h misc.c \
					  kmer.h kmer.c \
					  kmerhash.h kmerhash.c

kmerhash_test_LDADD = -lm

