
bin_PROGRAMS = fastqmd5 bammd5
check_PROGRAMS = random_fastq

TESTS = test_fastq

random_fastq_SOURCES = random_fastq.c

fastqmd5_SOURCES = fastqmd5.c md5.h md5.c
fastqmd5_LDADD = ../src/libquip.a

bammd5_SOURCES = bammd5.c md5.h md5.c
bammd5_LDADD= ../src/libquip.a


