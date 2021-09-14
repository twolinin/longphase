AC_INIT(longPhase,1.0)
AC_CONFIG_HEADERS([config.h])

m4_include([m4/hts_with_curses.m4])
m4_include([m4/ax_with_htslib.m4])

AC_PROG_AWK
AC_PROG_CC

AC_SYS_LARGEFILE

AX_WITH_HTSLIB
if test "$ax_cv_htslib" != yes; then
  AC_MSG_ERROR([HTSlib development files not found])
fi

Hsource=
Hinstall="#"
AC_SUBST([Hsource])
AC_SUBST([Hinstall])

HTS_WITH_CURSES

AS_IF([test "x$hts_late_cflags" != x],[CFLAGS="$CFLAGS $hts_late_cflags"])

AC_CONFIG_FILES([config.mk])
#AC_CONFIG_FILES([Makefile])
AC_OUTPUT

