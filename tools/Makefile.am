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

bin_PROGRAMS = printBCSR convertBCSR

printBCSR_SOURCES = printBCSR.F90
printBCSR_LDADD = $(top_builddir)/Modules/libfreeonmodules.la

convertBCSR_SOURCES = \
  convertBCSR.cc \
  index.h \
  lapack_interface.h \
  logger.cc \
  logger.h

convertBCSR_LDADD = libbcsr.la

EXTRA_DIST = $(bin_SCRIPTS)
