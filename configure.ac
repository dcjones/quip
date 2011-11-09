
AC_INIT( [fazer], [0.1], [dcjones@cs.washington.ed] )
AM_INIT_AUTOMAKE( [foreign -Wall -Werror] )
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AM_PROG_CC_C_O
AC_C_BIGENDIAN

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

AC_PROG_LIBTOOL

AC_CHECK_FUNC(fileno, ,
              AC_MSG_ERROR([The 'fileno' function is missing.]))

# check getopt
AC_CHECK_HEADER(getopt.h, ,
                AC_MSG_ERROR([The posix getopt.h header is needed.]))


# TODO: check for uint64_t


CXXFLAGS=$CFLAGS

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/version.h])

AM_CONFIG_HEADER([src/config.h])

AC_OUTPUT

