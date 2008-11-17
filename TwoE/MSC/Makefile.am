include $(top_srcdir)/make.inc

noinst_LTLIBRARIES = libTwoEMisc.la

libTwoEMisc_la_SOURCES = \
  ShellPairStruct.f90 \
  ShellPrint.f90 \
  Scratch.f90 \
  GetAtomPair.f90 \
  InitDbl.f

libTwoEMisc_la_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES)

# Dependencies.
GetAtomPair.lo : \
  ShellPairStruct.lo
