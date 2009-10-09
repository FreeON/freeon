include $(top_srcdir)/make.inc

# No linking of anything...
LIBS =

if USE_DYNAMIC_LINKING
lib_LTLIBRARIES = libfreeonlapackextra.la
else
noinst_LTLIBRARIES = libfreeonlapackextra.la
endif

libfreeonlapackextra_la_SOURCES = \
  ilaver.f \
	lsame.f \
	slamch.f \
	second_@LAPACK_TIMER@.f \
	dlamch.f \
	dsecnd_@LAPACK_TIMER@.f

noinst_PROGRAMS = testlsame testslamch testdlamch testsecond testdsecnd testieee testversion

testlsame_SOURCES = lsametst.f
testlsame_LDADD = libfreeonlapackextra.la

testslamch_SOURCES = slamchtst.f
testslamch_LDADD = libfreeonlapackextra.la

testdlamch_SOURCES = dlamchtst.f
testdlamch_LDADD = libfreeonlapackextra.la

testsecond_SOURCES = secondtst.f
testsecond_LDADD = libfreeonlapackextra.la

testdsecnd_SOURCES = dsecndtst.f
testdsecnd_LDADD = libfreeonlapackextra.la

testieee_SOURCES = tstiee.f

testversion_SOURCES = LAPACK_version.f
testversion_LDADD = libfreeonlapackextra.la

EXTRA_DIST = \
  dsecnd_EXT_ETIME.f \
  dsecnd_EXT_ETIME_.f \
  dsecnd_INT_CPU_TIME.f \
  dsecnd_INT_ETIME.f \
  dsecnd_NONE.f \
  second_EXT_ETIME.f \
  second_EXT_ETIME_.f \
  second_INT_CPU_TIME.f \
  second_INT_ETIME.f \
  second_NONE.f

slamch.lo: slamch.f
	@echo "compiling $< with no-optimization"
	@$(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77) $(AM_FFLAGS) $(FFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(EXTRA_INCLUDES) @NOOPT@ -c -o $@ $< > /dev/null 2>&1 && exit 0; \
	  echo "The following compiler command failed:" 1>&2; \
	  echo '$(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77) $(AM_FFLAGS) $(FFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(EXTRA_INCLUDES) @NOOPT@ -c -o $@ $<' 1>&2; \
	  $(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77) $(AM_FFLAGS) $(FFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(EXTRA_INCLUDES) @NOOPT@ -c -o $@ $< > /dev/null 2>&1

dlamch.lo: dlamch.f
	@echo "compiling $< with no-optimization"
	@$(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77) $(AM_FFLAGS) $(FFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(EXTRA_INCLUDES) @NOOPT@ -c -o $@ $< > /dev/null 2>&1 && exit 0; \
	  echo "The following compiler command failed:" 1>&2; \
	  echo '$(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77) $(AM_FFLAGS) $(FFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(EXTRA_INCLUDES) @NOOPT@ -c -o $@ $<' 1>&2; \
	  $(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77) $(AM_FFLAGS) $(FFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(EXTRA_INCLUDES) @NOOPT@ -c -o $@ $< > /dev/null 2>&1
