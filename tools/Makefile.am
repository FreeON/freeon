include $(top_srcdir)/make.inc

lib_LTLIBRARIES = libbcsr.la
libbcsr_la_SOURCES = \
  bcsr.cc \
  bcsr.h

include_HEADERS = bcsr.h

bin_SCRIPTS = \
  diffusionCoefficient.py \
  MDO2xyz.py \
  getFrame.py \
  get_energies.py \
  get_fMax.py \
  parseMDO.py \
  velocityAutocorrelation.py

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
