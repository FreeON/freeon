include ../../../make.inc

# No linking of anything...
LIBS =

noinst_LTLIBRARIES = libMondoLAPACKEXTRA.la

libMondoLAPACKEXTRA_la_SOURCES = \
  ilaver.f \
	lsame.f \
	slamch.f \
	second_@LAPACK_TIMER@.f \
	dlamch.f \
	dsecnd_@LAPACK_TIMER@.f

noinst_PROGRAMS = testlsame testslamch testdlamch testsecond testdsecnd testieee testversion

testlsame_SOURCES = lsametst.f
testlsame_LDADD = libMondoLAPACKEXTRA.la

testslamch_SOURCES = slamchtst.f
testslamch_LDADD = libMondoLAPACKEXTRA.la

testdlamch_SOURCES = dlamchtst.f
testdlamch_LDADD = libMondoLAPACKEXTRA.la

testsecond_SOURCES = secondtst.f
testsecond_LDADD = libMondoLAPACKEXTRA.la

testdsecnd_SOURCES = dsecndtst.f
testdsecnd_LDADD = libMondoLAPACKEXTRA.la

testieee_SOURCES = tstiee.f

testversion_SOURCES = LAPACK_version.f
testversion_LDADD = libMondoLAPACKEXTRA.la

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

slamch.o: slamch.f
	$(F77_COMMAND) @NOOPT@ -c -o $@ $<

#slamch.o: slamch.f
#	@echo "compiling $< with no-optimizations"
#	@$(F77_COMMAND) @NOOPT@ -c -o $@ $< && exit 0; \
#	  echo "The following comand failed:" 1>&2; \
#		echo "$(F77_COMMAND) @NOOPT@ -c -o $@ $<" 1>&2; \
#		$(F77_COMMAND) @NOOPT@ -c -o $@ $< > /dev/null 2>&1

dlamch.o: dlamch.f
	$(F77_COMMAND) @NOOPT@ -c -o $@ $<

#dlamch.o: dlamch.f
#	@echo "compiling $< with no-optimizations"
#	@$(F77_COMMAND) @NOOPT@ -c -o $@ $< && exit 0; \
#	  echo "The following comand failed:" 1>&2; \
#		echo "$(F77_COMMAND) @NOOPT@ -c -o $@ $<" 1>&2; \
#		$(F77_COMMAND) @NOOPT@ -c -o $@ $< > /dev/null 2>&1
