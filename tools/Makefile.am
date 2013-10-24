include $(top_srcdir)/make.inc

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
  bcsr.cc \
  bcsr.h \
  convertBCSR.cc \
  index.h \
  lapack_interface.h \
  logger.cc \
  logger.h

EXTRA_DIST = $(bin_SCRIPTS)
