include $(top_srcdir)/make.inc

noinst_LTLIBRARIES = libktrax.la libqctc.la

libktrax_la_SOURCES = KTrax.F

if AUTOMATIC_FFLAGS
# Force moderate compiler optimization flags on KTrax.F
libktrax_la_FFLAGS = $(shell echo $(AM_FFLAGS) | sed -e 's:\s*-ftree-vectorize\s*: :')
else
libktrax_la_FFLAGS = $(AM_FFLAGS)
endif

libqctc_la_LIBADD = libktrax.la
libqctc_la_SOURCES = \
  AA95c.F \
  AA95.F \
  BlokTrPdJ.F90 \
  CTraX10.F \
  CTraX12.F \
  CTraX14.F \
  CTraX4.F \
  CTraX6.F \
  CTraX8.F \
  CTraX.F \
  DEI.F \
  Density.F90 \
  DivPFFTensor.F90 \
  Gamma.F \
  Globals.F90 \
  Indexing.F90 \
  JGen.F90 \
  MD3TRR.F \
  MondoPoles77.F \
  MondoPoles.F90 \
  NukE.F90 \
  ParaQ.F90 \
  PBCFarField.F90 \
  PlaneWise.F90 \
  PFFTensor.F90 \
  PoleNodeType.F90 \
  PoleTree.F90 \
  QSortI.F \
  TreeWalk.F90 \
  XLate.F

bin_PROGRAMS = QCTC JForce MakePFFT

QCTC_SOURCES = \
  QCTC.F90 \
  Thresholds.F90

JForce_SOURCES = \
  JForce.F90 \
  Thresholds.F90

MakePFFT_SOURCES = \
  MakePFFT.F90 \
  Thresholds.F90

QCTC_DEPENDENCIES     = $(LIBRARY_DEPENDENCY_FILES)
JForce_DEPENDENCIES   = $(LIBRARY_DEPENDENCY_FILES)
MakePFFT_DEPENDENCIES = $(LIBRARY_DEPENDENCY_FILES)

SOURCES_DEPCHECK = \
  $(QCTC_SOURCES) \
  $(JForce_SOURCES) \
  $(MakePFFT_SOURCES)

EXTRA_DIST = \
  CTraX.Inc \
  HGToSP_Density.Inc \
  HGToSP.Inc \
  HGToSP_PoleNode.Inc \
  MACErrBnd4.Inc \
  MACErrBnd6.Inc \
  PBCTensor/Majik_Kubic_WS1.Inc \
  PBCTensor/Majik_Kubic_WS2.Inc

EXTRA_INCLUDES = -I$(top_builddir)/OneE -I$(top_srcdir)/Modules/MMA/LookUpTables_800_6x
LDADD = $(LDADD_MODULES) libqctc.la

# In case we have ctags installed, make the necessary tags files.
.PHONY : tags
tags :
	ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.{f,f90} ../Modules/*.{c,f,f90}

# Dependencies.
BlokTrPdJ.lo : \
  PoleTree.lo \
  TreeWalk.lo \
  Thresholds.o \
  PBCFarField.lo

DivPFFTensor.lo : \
  MondoPoles.lo \
  PFFTensor.lo

JForce.o : \
  $(top_builddir)/MondoConfig.h \
  libqctc.la

JGen.lo : \
  Globals.lo \
  PoleTree.lo \
  TreeWalk.lo \
  PBCFarField.lo \
  PlaneWise.lo

MakePFFT.o : \
  libqctc.la

NukE.lo : \
  Globals.lo \
  PoleTree.lo \
  TreeWalk.lo \
  PBCFarField.lo

PBCFarField.lo : \
  PFFTensor.lo \
  PoleTree.lo \
  MondoPoles.lo \
  Globals.lo \
  PFFTensor.lo

PFFTensor.lo : \
  MondoPoles.lo \
  PlaneWise.lo

ParaQ.lo : \
  Globals.lo \
  PoleTree.lo

PoleTree.lo : \
  Thresholds.o \
  MondoPoles.lo \
  Globals.lo \
  Density.lo

QCTC.o : \
  $(top_builddir)/MondoConfig.h \
  libqctc.la

TreeWalk.lo : \
  Globals.lo \
  PoleTree.lo

Density.lo : \
  Thresholds.o \
  Globals.lo

Thresholds.o : \
  PoleNodeType.lo

Indexing.lo : \
  PoleNodeType.lo

MondoPoles.lo : \
  Indexing.lo \
  Globals.lo

PlaneWise.lo : \
  MondoPoles.lo
