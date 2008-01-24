include ../../../make.inc

# No linking of anything...
LIBS =

noinst_PROGRAMS = testlsame testslamch testdlamch testsecond testdsecnd testieee testversion

testlsame_SOURCES = lsame.f lsametst.f
testslamch_SOURCES = slamch.f lsame.f slamchtst.f
testdlamch_SOURCES = dlamch.f lsame.f dlamchtst.f
testsecond_SOURCES = second_@LAPACK_TIMER@.f secondtst.f
testdsecnd_SOURCES = dsecnd_@LAPACK_TIMER@.f dsecndtst.f
testieee_SOURCES = tstiee.f
testversion_SOURCES = ilaver.f LAPACK_version.f

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
