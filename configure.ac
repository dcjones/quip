
AC_INIT([quip], [0.0], [dcjones@cs.washington.ed])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AM_PROG_CC_C_O
AC_C_BIGENDIAN

AC_TYPE_UINT8_T
AC_TYPE_INT32_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_SIZE_T

AC_HEADER_STDBOOL

opt_CFLAGS="--std=c99 -Wall -Wextra -pedantic -g -O3"
dbg_CFLAGS="--std=c99 -Wall -Wextra -pedantic -g -O0"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debugging],
                              [enable debugging info (default is no)])],
              [], [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])


AC_DEFINE([_FILE_OFFSET_BITS], [64],
          [Do not crash on >4GB files on 32bit machines.])

# LZMA configuration
AC_DEFINE([HAVE_ENCODER_LZMA1], 1, [LIBLZMA: LZMA1 Encoder])
AC_DEFINE([HAVE_ENCODER_LZMA2], 1, [LIBLZMA: LZMA2 Encoder])
AC_DEFINE([HAVE_DECODER_LZMA1], 1, [LIBLZMA: LZMA1 Decoder])
AC_DEFINE([HAVE_DECODER_LZMA2], 1, [LIBLZMA: LZMA2 Decoder])
AC_DEFINE([HAVE_CHECK_CRC32],   1, [LIBLZMA: CRC32 Check])
AC_DEFINE([HAVE_CHECK_CRC64],   1, [LIBLZMA: CRC64 Check])
AC_DEFINE([HAVE_CHECK_SHA256],  1, [LIBLZMA: SHA256 CHeck])
AC_DEFINE([HAVE_MF_BT2], 1, [LIBLZMA: BT2 Match Finder])
AC_DEFINE([HAVE_MF_BT3], 1, [LIBLZMA: BT3 Match Finder])
AC_DEFINE([HAVE_MF_BT4], 1, [LIBLZMA: BT4 Match Finder])
AC_DEFINE([HAVE_MF_HC3], 1, [LIBLZMA: HC3 Match Finder])
AC_DEFINE([HAVE_MF_HC4], 1, [LIBLZMA: HC4 Match Finder])

AC_PROG_LIBTOOL

AC_CHECK_FUNC(fileno, ,
              AC_MSG_ERROR([The 'fileno' function is missing.]))

# check getopt
AC_CHECK_HEADER(getopt.h, ,
                AC_MSG_ERROR([The posix getopt.h header is needed.]))


CXXFLAGS=$CFLAGS

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/lzma/Makefile
                 src/version.h])

AM_CONFIG_HEADER([src/config.h])

AC_OUTPUT

