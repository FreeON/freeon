include $(top_srcdir)/make.inc

noinst_LTLIBRARIES = libMSC.la

libMSC_la_SOURCES = \
  ShellPairStruct.F90 \
  ShellPrint.F90 \
  Scratch.F90 \
  GetAtomPair.F90 \
  InitDbl.F

# Dependencies.
GetAtomPair.lo : \
  ShellPairStruct.lo
