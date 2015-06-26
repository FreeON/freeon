include $(top_srcdir)/make.inc

lib_LTLIBRARIES = libbcsr.la
libbcsr_la_SOURCES = \
  bcsr.cc \
  bcsr.h

include_HEADERS = bcsr.h

bin_SCRIPTS = \
  diffusionCoefficient \
  MDO2xyz \
  getFrame \
  get_energies \
  get_fMax \
  parseMDO \
  velocityAutocorrelation

bin_PROGRAMS = \
	convertBCSR \
	plotDecay \
	printBCSR \
	standaloneBCSR

printBCSR_SOURCES = printBCSR.F90
printBCSR_LDADD = $(top_builddir)/Modules/libfreeonmodules.la $(LAPACK_LIBS) $(BLAS_LIBS)

standaloneBCSR_SOURCES = standaloneBCSR.F90
standaloneBCSR_LDADD = $(top_builddir)/Modules/libfreeonmodules.la $(LAPACK_LIBS) $(BLAS_LIBS)

plotDecay_SOURCES = plotDecay.F90
plotDecay_LDADD = $(top_builddir)/Modules/libfreeonmodules.la $(LAPACK_LIBS) $(BLAS_LIBS)

convertBCSR_SOURCES = \
  convertBCSR.cc \
  index.h \
  lapack_interface.h \
  logger.cc \
  logger.h

convertBCSR_LDADD = libbcsr.la $(LAPACK_LIBS) $(BLAS_LIBS)

EXTRA_DIST = $(bin_SCRIPTS)

TAGS_FILES = *F90 *.h *.cc ../Modules/*F90
CTAGSFLAGS = --sort=foldcase
