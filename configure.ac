AC_INIT(longPhase,1.0)
AC_CONFIG_HEADERS([config.h])

m4_include([m4/hts_with_curses.m4])
m4_include([m4/ax_with_htslib.m4])

AC_PROG_AWK
AC_PROG_CC

AC_SYS_LARGEFILE

if test "$ax_cv_htslib" != yes; then
  AC_MSG_ERROR([HTSlib development files not found])
fi

AC_CONFIG_FILES([config.mk])
AC_OUTPUT

