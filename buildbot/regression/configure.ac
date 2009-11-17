AC_INIT([regression], [1], [nbock@lanl.gov])
AM_INIT_AUTOMAKE

AC_PROG_CC

if test "${prefix}" = "NONE"; then
  AC_DEFINE_UNQUOTED([EXEPATH], ["${ac_default_prefix}/bin"], [The installation path])
else
  AC_DEFINE_UNQUOTED([EXEPATH], ["${prefix}/bin"], [The installation path])
fi

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
master/Makefile
slave_1/Makefile
slave_2/Makefile
slave_3/Makefile
])

AC_OUTPUT
