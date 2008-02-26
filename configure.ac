AC_INIT([MondoSCF], [1.0-alpha-latest], [mchalla@lanl.gov])
AC_COPYRIGHT([
This code is part of the MondoSCF suite of programs for linear scaling 
electronic structure theory and ab initio molecular dynamics.

Copyright (2004). The Regents of the University of California. This 
material was produced under U.S. Government contract W-7405-ENG-36 
for Los Alamos National Laboratory, which is operated by the University 
of California for the U.S. Department of Energy. The U.S. Government has 
rights to use, reproduce, and distribute this software.  NEITHER THE 
GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED, 
OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  

This program is free software; you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by the 
Free Software Foundation; either version 2 of the License, or (at your 
option) any later version. Accordingly, this program is distributed in 
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. See the GNU General Public License at www.gnu.org for details.

While you may do as you like with this software, the GNU license requires
that you clearly mark derivative software.  In addition, you are encouraged 
to return derivative works to the MondoSCF group for review, and possible 
disemination in future releases.
])

AM_INIT_AUTOMAKE

if test -f localversion; then
  . localversion
else
  localversion="not-set"
fi

AC_DEFINE_UNQUOTED([LOCAL_VERSION], ["${localversion}"], [The local version of Mondo])
AC_MSG_NOTICE([configuring local version ${localversion}])

# Fix up the header definitions for hdf.
AC_OUTPUT_COMMANDS([
  echo "creating H5config.h"
  sed 's/#define /#define H5_/' < config.h | \
    sed 's/#undef /#undef H5_/' > tempconf
  if test ! -f H5config.h; then
    mv -f tempconf H5config.h
  elif (diff tempconf H5config.h > /dev/null); then
    rm -f tempconf
    echo "H5config.h is unchanged"
  else
    mv -f tempconf H5config.h
  fi


  # We don't want inline defined for C++ compilers
  # Don't worry about the C++ ifdef wrappers in the H5pubconf file, since
  # 'H5_inline' isn't a C++ keyword.
  echo "adding inline statements to H5config.h"
  cat >> H5config.h <<EOF

#if defined(__cplusplus) && defined(inline)
#undef inline
#endif
EOF

  echo "creating MondoConfig.h"
  grep "#define" config.h > tempconf
  if test ! -f MondoConf.h; then
    mv -f tempconf MondoConfig.h
  elif (diff tempconf MondoConfig.h > /dev/null); then
    rm -f tempconf
    echo "MondoConfig.h is unchanged"
  else
    mv -f tempconf MondoConf.h
  fi
])

# Set the optimization flags.
if test -z "${CFLAGS}"; then
  AC_MSG_NOTICE([Setting default optimization: CFLAGS = ${CFLAGS:=-O2}])
fi
if test -z "${FFLAGS}"; then
  AC_MSG_NOTICE([Setting default optimization: FFLAGS = ${FFLAGS:=-O2}])
fi
if test -z "${FCFLAGS}"; then
  AC_MSG_NOTICE([Setting default optimization: FCFLAGS = ${FCFLAGS:=-O2}])
fi

# Find compilers.
AC_PROG_CC
AC_PROG_FC
AC_PROG_F77

# Put the C-compiler into standard C mode.
AC_PROG_CC_STDC
AC_GNU_SOURCE

# Set line-length limits.
AC_MSG_CHECKING([how to set line length limits of f90 compiler])
if test "`basename ${FC}`" = "gfortran"; then
  AC_MSG_RESULT([gfortran detected, setting unlimited line-length limit (-ffree-line-length-none)])
  FCFLAGS="${FCFLAGS} -ffree-line-length-none"
else
  AC_MSG_RESULT([not needed])
fi

# Check for ranlib.
#AC_PROG_RANLIB

# Check for libtool.
AC_PROG_LIBTOOL

# Set $(MAKE).
AC_PROG_MAKE_SET

# Check the build system.
AC_CANONICAL_BUILD

# Check for some standard headers.
AC_LANG(C)
AC_HEADER_STDC
AC_CHECK_HEADERS(sys/resource.h sys/time.h unistd.h sys/ioctl.h sys/types.h \
  math.h sys/wait.h)

if test "${ac_cv_header_sys_wait_h}" = "yes"; then

  AC_MSG_CHECKING([for WSTOPSIG macro])
  AC_LINK_IFELSE([
#include <sys/wait.h>

int
main ()
{
  int status = 0;

  if (WSTOPSIG(status) == 0) { return 0; }
  else { return -1; }
}],
  [AC_DEFINE(HAVE_WSTOPSIG, 1, [If WIFSTOPPED is true of status, this macro
   returns the signal number of the signal that caused the child process to
   stop.])
  AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])])

  AC_MSG_CHECKING([for WCOREDUMP macro])
  AC_LINK_IFELSE([
#include <sys/wait.h>

int
main ()
{
  int status = 0;

  if (WCOREDUMP(status) == 0) { return 0; }
  else { return -1; }
}],
  [AC_DEFINE(HAVE_WCOREDUMP, 1, [This macro returns a nonzero value if the child
   process terminated and produced a core dump.])
  AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])])

  AC_MSG_CHECKING([for WEXITSTATUS macro])
  AC_LINK_IFELSE([
#include <sys/wait.h>

int
main ()
{
  int status = 0;

  if (WEXITSTATUS(status) == 0) { return 0; }
  else { return -1; }
}],
  [AC_DEFINE(HAVE_WEXITSTATUS, 1, [If WIFEXITED is true of status, this macro
   returns the low-order 8 bits of the exit status value from the child
   process.])
  AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])])

fi

AC_MSG_CHECKING([for rand48])
AC_LINK_IFELSE([
#include <stdlib.h>

int
main ()
{
  srand48(0);
  return 0;
}],
  [AC_DEFINE(HAVE_SRAND48, 1, [These functions generate pseudo-random numbers
   using the linear congruential algorithm and 48-bit integer arithmetic.])
  AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])])

# Check for phipac.
AC_MSG_CHECKING([whether to use PHiPAC])
AC_ARG_ENABLE([phipac],
  [AS_HELP_STRING([--enable-phipac],
  [Enable the use of phipac optimized dgemms.])],
  [AC_DEFINE(PHIPAC, 1, [Use PHiPAC])
  AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])])

# Say something about PHiPAC/BLAS.
AC_MSG_NOTICE([*** NOTICE ***
The only part of MondoSCF affected by this flag is MondoMods/SMWrap.C. If you
have built MondoSCF already and you are simply changing from using PHiPAC to
using BLAS or vice versa, you don't have to rebuild the whole code. Simply
delete MondoMods/SMWrap.o and run make again.  This will rebuild all that is
necessary.])

# Test some stuff for C-code.
AC_LANG(C)

AC_CHECK_LIB([m], [round], [],
  [AC_MSG_FAILURE([we need the round() function])])

# Switch tests to fortran.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM

# Find whether the compiler is gfortran and too old...
if test "`basename ${FC}`" = "gfortran"; then
  AC_MSG_CHECKING([for possible gfortran compiler problems])
  AC_LINK_IFELSE([
PROGRAM test_compiler

  TYPE TestType
    INTEGER, ALLOCATABLE, DIMENSION(:) :: I
  END TYPE TestType

END PROGRAM test_compiler
    ],
    [AC_MSG_RESULT([good, none found])],
    [AC_MSG_RESULT([failed test])
    AC_MSG_FAILURE([the compiler does not understand the ALLOCATABLE attribute inside a type definition])])
fi

# Find flags for static linking.
AC_ARG_ENABLE([static],
  [AS_HELP_STRING([--disable-static],
  [Disable static linking of MondoSCF executables.])],
  [enable_static_linking="${enableval}"],
  [enable_static_linking="yes"])

static_linking_flags=""

if test "${enable_static_linking}" = "yes"; then

  if test "`basename ${FC}`" = "gfortran"; then
    static_linking_flags="-static"
    AC_MSG_CHECKING([identified gfortran compiler, trying ${static_linking_flags} flag])
  elif test "`basename ${FC}`" = "ifort"; then
    static_linking_flags="-static"
    AC_MSG_CHECKING([identified intel ifort compiler, trying ${static_linking_flags} flag])
  else
    AC_MSG_NOTICE([unknown compiler, no static linking flags found])
  fi

  if test -n "${static_linking_flags}"; then

    old_LDFLAGS="${LDFLAGS}"
    if test -z "${LDFLAGS}"; then
      LDFLAGS="${static_linking_flags}"
    else
      LDFLAGS="${LDFLAGS} ${static_linking_flags}"
    fi

    # Test the linking flags.
    AC_LINK_IFELSE([
      program link_test

        write(*,*) "hmmmm"

      end program link_test
      ],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no, not linking statically])
      LDFLAGS="${old_LDFLAGS}"])

    unset old_LDFLAGS

  fi

else
  AC_MSG_NOTICE([static linking disabled])
fi

# Check for internal LAPACK/BLAS.
AC_MSG_CHECKING([for internal lapack/blas])
AC_ARG_ENABLE([internal-lapack],
  [AS_HELP_STRING([--disable-internal-lapack],
  [Disable the internal lapack/blas (3.1.1) and use your own. You might have to
  use the options --with-lapack=DIR or --with-blas=DIR to specify exactly where
  your library is.])],
  [use_internal_lapack="${enableval}"],
  [use_internal_lapack="yes"])

AC_ARG_ENABLE([internal-lapack-checks],
  [AS_HELP_STRING([--enable-internal-lapack-checks],
  [Enable the checks in the internal lapack. Lapack can be checked by running
   "make check".])],
  [compile_internal_lapack_checks="${enableval}"],
  [compile_internal_lapack_checks="no"])

AM_CONDITIONAL(INTERNAL_LAPACK, [test "${use_internal_lapack}" = "yes"])
AM_CONDITIONAL(INTERNAL_LAPACK_CHECKS, [test "${compile_internal_lapack_checks}" = "yes"])

LAPACK_TIMER="NONE"

if test "${use_internal_lapack}" = "yes"; then
  AC_MSG_RESULT([yes])
  INTERNAL_LAPACK_LIBS="-lMondoLAPACK -lMondoBLAS"

  # Test for no-optimization compiler flag.
  AC_MSG_CHECKING([for no-optmiziation flags])
  if test "`basename ${FC}`" = "gfortran"; then
    NOOPT="-O0"
    LAPACK_TIMER="INT_ETIME"
    AC_MSG_RESULT([identified gfortran compiler, using ${NOOPT} no-optimization flag])
  elif test "`basename ${FC}`" = "ifort"; then
    NOOPT="-O0 -fltconsistency -fp-port"
    LAPACK_TIMER="EXT_ETIME"
    AC_MSG_RESULT([identified intel ifort compiler, using ${NOOPT} no-optimization flag])
  else
    AC_MSG_RESULT([unknown compiler, no no-optimization flags found])
  fi
  AC_SUBST(NOOPT, [${NOOPT}])

else
  AC_MSG_RESULT([no])
  INTERNAL_LAPACK_LIBS=""

  AC_ARG_WITH([lapack],
    [AS_HELP_STRING([--with-lapack=DIR],
    [Specify non-standard directory for the BLAS/LAPACK libraries.])],
    [LDFLAGS="${LDFLAGS} -L${withval}/lib"
    AC_MSG_NOTICE([using -L${withval}/lib])],
    [])

  AC_ARG_WITH([lapacklibs],
    [AS_HELP_STRING([--with-lapacklibs=-lblas -llapack...],
    [Specify the blas/lapack libraries as specified on the linker command line (default is -lblas -llapack)])],
    [more_LIBS="${withval}"
    AC_MSG_NOTICE([using ${withval}])],
    [more_LIBS="-lblas -llapack"])

  if test -z "${LIBS}"; then
    LIBS="${more_LIBS}"
  else
    LIBS="${LIBS} ${more_LIBS}"
  fi
  unset more_LIBS

  # Check for BLAS.
  AC_MSG_CHECKING([for BLAS (dgemm)])
  AC_LINK_IFELSE([
  program test

    call dgemm()

  end program test
  ],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_FAILURE([no])])

  # Check for LAPACK.
  AC_MSG_CHECKING([for LAPACK (dsyev)])
  AC_LINK_IFELSE([
  program test

    call dsyev()

  end program test
  ],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_FAILURE([no])])
fi

AC_SUBST(LAPACK_TIMER, [${LAPACK_TIMER}])

# Check for internal hdf5.
AC_MSG_CHECKING([for internal hdf5])
AC_ARG_ENABLE([internal-hdf5],
  [AS_HELP_STRING([--disable-internal-hdf5],
  [Disable the internal hdf5 (1.4.5-post2) and use your own. You might have to
   use the option --with-hdf5=DIR to specify where exactly your hdf5 library
   is.])],
  [use_internal_hdf5="${enableval}"],
  [use_internal_hdf5="yes"])

AM_CONDITIONAL(INTERNAL_HDF5, [test "${use_internal_hdf5}" = "yes"])

if test "${use_internal_hdf5}" = "yes"; then
  AC_MSG_RESULT([yes])
  INTERNAL_HDF5_LIBS="-lMondoHDF5"
  AC_DEFINE(HAVE_INTERNAL_HDF5, 1, [Use the inernal hdf5 library])

  HDF5_CPPFLAGS="-UH5_DEBUG_API -DNDEBUG"

  case "${build_os}" in
    *linux*)

    HDF5_CPPFLAGS="${HDF5_CPPFLAGS} -D_POSIX_SOURCE -D_BSD_SOURCE"
    ;;
  esac

  AC_SUBST(HDF5_CPPFLAGS, [${HDF5_CPPFLAGS}])

  # Store old CPPFLAGS.
  old_CPPFLAGS=${CPPFLAGS}
  CPPFLAGS="${old_CPPFLAGS} ${HDF5_CPPFLAGS}"

  # Set language to C.
  AC_LANG(C)

  # Check for system header files.
  AC_HEADER_TIME

  AC_CHECK_HEADERS(sys/ioctl.h sys/stat.h)
  AC_CHECK_HEADERS(sys/socket.h)
  AC_CHECK_HEADERS(stddef.h setjmp.h)
  AC_CHECK_HEADERS(stdint.h, C9x=yes)

  # Data types and their sizes.
  AC_TYPE_OFF_T
  AC_CHECK_TYPE(size_t, unsigned long)
  AC_CHECK_TYPE(ssize_t, long)
  AC_C_BIGENDIAN
  AC_CHECK_SIZEOF(char, 1)
  AC_CHECK_SIZEOF(short, 2)
  AC_CHECK_SIZEOF(int, 4)
  AC_CHECK_SIZEOF(long, 4)
  AC_CHECK_SIZEOF(long long, 8)
  AC_CHECK_SIZEOF(__int64, 8)
  AC_CHECK_SIZEOF(float, 4)
  AC_CHECK_SIZEOF(double, 8)
  AC_CHECK_SIZEOF(long double, 8)

  AC_CHECK_SIZEOF(        int8_t, 1)
  AC_CHECK_SIZEOF(       uint8_t, 1)
  AC_CHECK_SIZEOF(  int_least8_t, 1)
  AC_CHECK_SIZEOF( uint_least8_t, 1)
  AC_CHECK_SIZEOF(   int_fast8_t, 1)
  AC_CHECK_SIZEOF(  uint_fast8_t, 1)

  AC_CHECK_SIZEOF(       int16_t, 2)
  AC_CHECK_SIZEOF(      uint16_t, 2)
  AC_CHECK_SIZEOF( int_least16_t, 2)
  AC_CHECK_SIZEOF(uint_least16_t, 2)
  AC_CHECK_SIZEOF(  int_fast16_t, 2)
  AC_CHECK_SIZEOF( uint_fast16_t, 2)

  AC_CHECK_SIZEOF(       int32_t, 4)
  AC_CHECK_SIZEOF(      uint32_t, 4)
  AC_CHECK_SIZEOF( int_least32_t, 4)
  AC_CHECK_SIZEOF(uint_least32_t, 4)
  AC_CHECK_SIZEOF(  int_fast32_t, 4)
  AC_CHECK_SIZEOF( uint_fast32_t, 4)

  AC_CHECK_SIZEOF(       int64_t, 8)
  AC_CHECK_SIZEOF(      uint64_t, 8)
  AC_CHECK_SIZEOF( int_least64_t, 8)
  AC_CHECK_SIZEOF(uint_least64_t, 8)
  AC_CHECK_SIZEOF(  int_fast64_t, 8)
  AC_CHECK_SIZEOF( uint_fast64_t, 8)

  AC_CHECK_SIZEOF(size_t, 4)
  AC_CHECK_SIZEOF(ssize_t, 4)

  AC_CHECK_SIZEOF(off_t, 4)

  # ----------------------------------------------------------------------
  # Try to figure out how to print `long long'.  Some machines use `%lld'
  # and others use `%qd'.  There may be more!  The final `l' is a
  # default in case none of the others work.
  # Need to patch up LD_LIBRARY_PATH so that the execution can find all
  # the dynamic library.  The correct way to do it should be updating
  # LD_LIBRARY_PATH along with LDFLAGS or do it with the AC_TRY_RUN macro.
  AC_MSG_CHECKING(how to print long long)
  AC_CACHE_VAL(hdf5_cv_printf_ll,
    LD_LIBRARY_PATH="$LD_LIBRARY_PATH`echo $LDFLAGS | sed -e 's/-L/:/g' -e 's/ //g'`"
    export LD_LIBRARY_PATH

    for hdf5_cv_printf_ll in l L q ll unknown; do
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
	char *s = malloc(128);
	long long x = (long long)1048576 * (long long)1048576;
	sprintf(s,"%${hdf5_cv_printf_ll}d",x);
	exit(strcmp(s,"1099511627776"));
}
      ], break,,continue)
    done)

  AC_MSG_RESULT(%${hdf5_cv_printf_ll}d and %${hdf5_cv_printf_ll}u)
  AC_DEFINE_UNQUOTED(PRINTF_LL_WIDTH, "$hdf5_cv_printf_ll",
    [Width for printf() for type 'long long' or '__int64', us. 'll'])

  AC_CHECK_HEADERS(pthread.h)
  AC_CHECK_LIB(pthread, pthread_create,,)

  # ----------------------------------------------------------------------
  # How does one figure out the local time zone?  Anyone know of a
  # Posix way to do this?
  #

  # First check if `struct tm' has a 'tm_gmtoff' member.
  AC_MSG_CHECKING(for tm_gmtoff in struct tm)
  AC_TRY_COMPILE([
#include <sys/time.h>
#include <time.h>
  ],
  [struct tm tm; tm.tm_gmtoff=0;],
  AC_DEFINE(HAVE_TM_GMTOFF, 1,
    [Define if 'tm_gmtoff' is a member of 'struct tm'])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  # check if `struct tm' has a `__tm_gmtoff' member.
  AC_MSG_CHECKING(for __tm_gmtoff in struct tm)
  AC_TRY_COMPILE([
#include <sys/time.h>
#include <time.h>],[struct tm tm; tm.__tm_gmtoff=0;
  ],
  AC_DEFINE(HAVE___TM_GMTOFF, 1,
    [Define if '__tm_gmtoff' is a member of 'struct tm'])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  # Check whether the global variable 'timezone' is defined.
  AC_MSG_CHECKING(for global timezone variable)
  AC_TRY_LINK([
#include <sys/time.h>
#include <time.h>], [timezone=0;
  ],
  AC_DEFINE(HAVE_TIMEZONE, 1,
    [Define if 'timezone' is a global variable])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  # Check whether `struct timezone' is defined.
  AC_STRUCT_TIMEZONE
  AC_MSG_CHECKING(for struct timezone)
  AC_TRY_COMPILE([
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
  ],
  [struct timezone tz; tz.tz_minuteswest=0;],
  AC_DEFINE(HAVE_STRUCT_TIMEZONE, 1,
    [Define if 'struct timezone' is defined])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  # ----------------------------------------------------------------------
  # How do we figure out the width of a tty in characters?
  #
  AC_CHECK_FUNCS(_getvideoconfig gettextinfo GetConsoleScreenBufferInfo)
  AC_CHECK_FUNCS(_scrsize ioctl)

  AC_MSG_CHECKING(for struct videoconfig)
  AC_TRY_COMPILE(,[struct videoconfig w; w.numtextcols=0;],
  AC_DEFINE(HAVE_STRUCT_VIDEOCONFIG, 1,
    [Define if 'struct videoconfig' is defined])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  AC_MSG_CHECKING(for struct text_info)
  AC_TRY_COMPILE(,[struct text_info w; w.screenwidth=0;],
  AC_DEFINE(HAVE_STRUCT_TEXT_INFO, 1,
    [Define if 'struct text_info' is defined])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  AC_MSG_CHECKING(for TIOCGWINSZ)
  AC_TRY_COMPILE([#include <sys/ioctl.h>],[int w=TIOCGWINSZ;],
  AC_DEFINE(HAVE_TIOCGWINSZ, 1,
    [Define if the ioctl TIOGWINSZ is defined])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  AC_MSG_CHECKING(for TIOCGGETD)
  AC_TRY_COMPILE([#include <sys/ioctl.h>],[int w=TIOCGETD;],
  AC_DEFINE(HAVE_TIOCGETD, 1,
    [Define if the ioctl TIOCGETD is defined])
  AC_MSG_RESULT(yes),
  AC_MSG_RESULT(no))

  AC_CHECK_FUNCS(difftime fork gethostname getpwuid getrusage)
  AC_CHECK_FUNCS(gettimeofday BSDgettimeofday longjmp setsysinfo sigaction)
  AC_CHECK_FUNCS(signal snprintf vsnprintf strdup system waitpid)

  # ----------------------------------------------------------------------
  # Check compiler characteristics
  #
  AC_C_CONST
  AC_C_INLINE

  AC_MSG_CHECKING(for __attribute__ extension)
  AC_TRY_COMPILE(,[int __attribute__((unused)) x],
    AC_DEFINE(HAVE_ATTRIBUTE, 1,
    [Define if the __attribute__(()) extension is present])
    AC_MSG_RESULT(yes),
    AC_MSG_RESULT(no))

  AC_MSG_CHECKING(for __FUNCTION__ extension)
  AC_TRY_COMPILE(,[(void)__FUNCTION__],
    AC_DEFINE(HAVE_FUNCTION, 1,
    [Define if the compiler understand the __FUNCTION__ keyword])
    AC_MSG_RESULT(yes),
    AC_MSG_RESULT(no))

  # Reset the CPPFLAGS.
  CPPFLAGS=${old_CPPFLAGS}

else

  AC_MSG_RESULT([no])
  INTERNAL_HDF5_LIBS=""

  AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5=DIR],
    [Specify non-standard directory for the hdf5 libraries.])],
    [LDFLAGS="${LDFLAGS} -L${withval}/lib"
    CPPFLAGS="${CPPFLAGS} -I${withval}/include"
    AC_MSG_NOTICE([using ${withval}/lib and ${withval}/include])],
    [])

  # Set language to C.
  AC_LANG(C)

  # Check for math library (ceil).
  AC_CHECK_LIB([m], [ceil],
    [],
    [AC_MSG_NOTICE([The external hdf5 library might require ceil()])])

  # Check for zlib compression.
  AC_CHECK_LIB([z], [zlibVersion],
    [],
    [AC_MSG_NOTICE([The external hdf5 library might require the zlib library.])])

  # Check for the pthread library.
  AC_CHECK_HEADERS(pthread.h)
  AC_CHECK_LIB([pthread], [pthread_create],
    [],
    [AC_MSG_NOTICE([The external hdf5 library might require the pthread library.])])

  # Check for HDF5 libraries.
  AC_CHECK_HEADER([hdf5.h],
    [],
    [AC_MSG_FAILURE([I need the hdf5.h header file])])
  AC_CHECK_LIB([hdf5], [H5get_libversion],
    [],
    [AC_MSG_FAILURE([I need the hdf5 library])])

fi

# Set the preprocessor.
CPP=cpp

# Set some preprocessor command line options.
FORTRAN_CPPFLAGS="-P -C -xc -traditional -I../MondoMods"

if test -z "${LIBS}"; then
  LIBS="${INTERNAL_LAPACK_LIBS} ${INTERNAL_HDF5_LIBS}"
else
  LIBS="${LIBS} ${INTERNAL_LAPACK_LIBS} ${INTERNAL_HDF5_LIBS}"
fi

# Set this stuff in the makefiles.
AC_SUBST(CPP, ${CPP})
AC_SUBST(AM_CPPFLAGS, ${AM_CPPFLAGS})
AC_SUBST(FORTRAN_CPPFLAGS, ${FORTRAN_CPPFLAGS})

AC_MSG_CHECKING([whether we should remove preprocessed files])
AC_ARG_ENABLE([delete-preprocessed-files],
  [AS_HELP_STRING([--disable-delete-preprocessed-files],
  [Disable the deleting of preprocessed files. This is useful for debugging...])],
  [remove_preprocessed_files="${enableval}"],
  [remove_preprocessed_files="yes"])
AC_MSG_RESULT([${remove_preprocessed_files}])
AC_SUBST(REMOVE_PREPROCESSED_FILES, ${remove_preprocessed_files})

# Some mondo variables.
AC_DEFINE_UNQUOTED(HAVE_MONDO_HOME, ["${prefix}"], [Define the main installation directory.])
AC_DEFINE_UNQUOTED(HAVE_MONDO_PLAT, ["${build}"], [Define build platform.])
AC_DEFINE_UNQUOTED(HAVE_MONDO_HOST, ["`hostname`"], [Define hostname.])
AC_DEFINE_UNQUOTED(HAVE_MONDO_MACH, ["${build_cpu}"], [Define build cpu.])
AC_DEFINE_UNQUOTED(HAVE_MONDO_SYST, ["${build_os}"], [Define build system (os).])
AC_DEFINE_UNQUOTED(HAVE_MONDO_VRSN, ["`uname -r`"], [Define build os version.])

# Enable fix for gfortran.
if test "${ac_cv_fc_compiler_gnu}" = "yes"; then
  AC_DEFINE(GFORTRAN, 1, [Enable fix for gfortran])
fi

# Output all of it....
AC_CONFIG_FILES([
Makefile
MondoSCF/Makefile
MondoMods/Makefile
MondoMods/lapack/Makefile
MondoMods/lapack/install/Makefile
MondoMods/lapack/blas/Makefile
MondoMods/lapack/lapack/Makefile
MondoMods/lapack/blas_testing/Makefile
MondoMods/lapack/lapack_testing/Makefile
MondoMods/lapack/lapack_testing/lin/Makefile
MondoMods/lapack/lapack_testing/eig/Makefile
MondoMods/lapack/lapack_testing/matgen/Makefile
MondoMods/hdf5/Makefile
SCFeqs/Makefile
OneE/Makefile
QCTC/Makefile
QCTC/LibQCTC/Makefile
HiCu/Makefile
TwoE/Makefile
TwoE/GRD/Makefile
TwoE/HRR/Makefile
TwoE/INT/Makefile
TwoE/MSC/Makefile
TwoE/VRR/Makefile
ONX/Makefile
])

AC_CONFIG_HEADERS([config.h])

# Configuration date
AC_DEFINE_UNQUOTED(CONFIG_DATE, ["`date`"], [The date this configuration script was run])

# User doing the configuration
AC_DEFINE_UNQUOTED(CONFIG_USER, ["`whoami`@`hostname`"], [The user running the configure script])
if test -n "$ORGANIZATION"; then
  AC_DEFINE_UNQUOTED(CONFIG_USER, ["$CONFIG_USER at $ORGANIZATION"], [The user running the configure script])
fi

AC_MSG_NOTICE([***   Summary   ***])
AC_MSG_NOTICE([CC            = ${CC}])
AC_MSG_NOTICE([FC            = ${FC}])
AC_MSG_NOTICE([F77           = ${F77}])
AC_MSG_NOTICE([CFLAGS        = ${CFLAGS}])
AC_MSG_NOTICE([FCFLAGS       = ${FCFLAGS}])
AC_MSG_NOTICE([FFLAGS        = ${FFLAGS}])
AC_MSG_NOTICE([LDFLAGS       = ${LDFLAGS}])
AC_MSG_NOTICE([LIBS          = ${LIBS}])
AC_MSG_NOTICE([CPPFLAGS      = ${CPPFLAGS}])
AC_MSG_NOTICE([HDF5_CPPFLAGS = ${HDF5_CPPFLAGS}])
AC_MSG_NOTICE([])
if test "${use_internal_lapack}" = "yes"; then
  AC_MSG_NOTICE([Using internal lapack])
else
  AC_MSG_NOTICE([Using external lapack library])
fi
if test "${use_internal_hdf5}" = "yes"; then
  AC_MSG_NOTICE([Using internal hdf5])
else
  AC_MSG_NOTICE([Using external hdf5 library])
fi
if test "${remove_preprocessed_files}" = "yes"; then
  AC_MSG_NOTICE([I will remove preprocessed files])
else
  AC_MSG_NOTICE([I will _not_ remove preprocessed files])
fi
AC_MSG_NOTICE([*** End Summary ***])

AC_OUTPUT
