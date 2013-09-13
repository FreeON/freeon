include $(top_srcdir)/make.inc

bin_SCRIPTS = \
  diffusionCoefficient.py \
  MDO2xyz.py \
  getFrame.py \
  get_energies.py \
  get_fMax.py \
  parseMDO.py \
  velocityAutocorrelation.py

bin_PROGRAMS = BCSRInfo

BCSRInfo_SOURCES = BCSRInfo.F90
BCSRInfo_LDADD = $(top_builddir)/Modules/libfreeonmodules.la

EXTRA_DIST = $(bin_SCRIPTS)
