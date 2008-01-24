include ../../make.inc

noinst_LIBRARIES = libTwoEMisc.a

libTwoEMisc_a_SOURCES = \
  ShellPairStruct.f90 \
  ShellPrint.f90 \
  Scratch.f90 \
  GetAtomPair.f90 \
  InitDbl.f

libTwoEMisc_a_DEPENDENCIES = $(foreach i, $(MondoMods_DEPENDENCY_FILES), ../$(i))

EXTRA_INCLUDES = -I../../MondoMods

# Dependencies.
GetAtomPair.o : \
  ShellPairStruct.o
