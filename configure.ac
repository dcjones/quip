
AC_INIT([quip],
        m4_esyscmd([./git-version-gen | tr -d '\n']),
        [dcjones@cs.washington.ed])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_LN_S
AC_PROG_RANLIB
AM_PROG_AR

AC_C_BIGENDIAN

AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_INT32_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_SIZE_T

AC_HEADER_STDBOOL

# TODO: if zlib is not found, compiile without BAM support
AC_CHECK_LIB(z, compress, , AC_MSG_ERROR([zlib is required but not was not found.]))

# Check if the compiler has a builtin prefetch function
AC_LANG(C)
AC_MSG_CHECKING([for __builtin_prefetch])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
    [[#include<stdlib.h>]],
    [[__builtin_prefetch(NULL, 0, 0);]])],
  [
    AC_MSG_RESULT([yes])
    AC_DEFINE_UNQUOTED([HAVE_PREFETCH], 1, [Define to 1 if you have the `__builtin_prefetch' function.] ) ],
  [
    AC_MSG_RESULT([no])
    AC_DEFINE_UNQUOTED([HAVE_PREFETCH], 0, [Define to 1 if you have the `__builtin_prefetch' function.] ) ],
  ])

opt_CFLAGS="--std=gnu99 -Wall -Wextra -pedantic -g -O3 -DNDEBUG"
dbg_CFLAGS="--std=gnu99 -Wall -Wextra -pedantic -g -O0"

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging info (default is no)])],
              [], [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes],
      [CFLAGS="$dbg_CFLAGS"],
      [CFLAGS="$opt_CFLAGS"])

ACX_PTHREAD()
LIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"

AC_DEFINE([_FILE_OFFSET_BITS], [64],
          [Do not crash on >4GB files on 32bit machines.])

# These are POSIX standards, but we check just to be careful.
AC_CHECK_FUNC(open, ,
              AC_MSG_ERROR([The 'open' function is missing.]))
AC_CHECK_FUNC(fdopen, ,
              AC_MSG_ERROR([The 'fdopen' function is missing.]))
AC_CHECK_FUNC(fileno, ,
              AC_MSG_ERROR([The 'fileno' function is missing.]))
AC_CHECK_FUNC(unlink, ,
              AC_MSG_ERROR([The 'unlink' function is missing.]))
AC_CHECK_FUNC(isatty, ,
              AC_MSG_ERROR([The 'isatty' function is missing.]))

AC_CHECK_HEADER(getopt.h, ,
                AC_MSG_ERROR([The posix getopt.h header is needed.]))


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/version.h])

AM_CONFIG_HEADER([src/config.h])

AC_OUTPUT

