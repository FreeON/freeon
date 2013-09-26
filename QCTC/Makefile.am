include $(top_srcdir)/make.inc

noinst_LTLIBRARIES = libktrax.la

libktrax_la_SOURCES = KTrax.F

if AUTOMATIC_FFLAGS
# Force moderate compiler optimization flags on KTrax.F
libktrax_la_FFLAGS = -O2 -g -pipe
else
libktrax_la_FFLAGS = $(AM_FFLAGS)
endif

bin_PROGRAMS = QCTC JForce MakePFFT

QCTC_SOURCES = \
  AA95c.F \
  AA95.F \
  CTraX10.F \
  CTraX12.F \
  CTraX14.F \
  CTraX4.F \
  CTraX6.F \
  CTraX8.F \
  CTraX.F \
  Density.F90 \
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
  PoleNodeType.F90 \
  PoleTree.F90 \
  QCTC.F90 \
  Thresholds.F90 \
  TreeWalk.F90 \
  XLate.F

JForce_SOURCES = \
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
  Density.F90 \
  Gamma.F \
  Globals.F90 \
  Indexing.F90 \
  JForce.F90 \
  MD3TRR.F \
  MondoPoles77.F \
  MondoPoles.F90 \
  ParaQ.F90 \
  PBCFarField.F90 \
  PoleNodeType.F90 \
  PoleTree.F90 \
  Thresholds.F90 \
  TreeWalk.F90 \
  XLate.F

MakePFFT_SOURCES = \
  XLate.F \
  CTraX.F \
  Globals.F90 \
  PoleNodeType.F90 \
  Thresholds.F90 \
  Indexing.F90 \
  QSortI.F \
  MondoPoles.F90 \
  DEI.F \
  PlaneWise.F90 \
  PFFTensor.F90 \
  DivPFFTensor.F90 \
  ParaQ.F90 \
  MakePFFT.F90


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
LDADD = $(LDADD_MODULES) libktrax.la

# In case we have ctags installed, make the necessary tags files.
.PHONY : tags
tags :
	ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.{f,f90} ../Modules/*.{c,f,f90}

# Dependencies.
BlokTrPdJ.o : \
  PoleTree.o \
  TreeWalk.o \
  Thresholds.o \
  PBCFarField.o

DivPFFTensor.o : \
  MondoPoles.o \
  PFFTensor.o

JForce.o : \
  $(top_builddir)/MondoConfig.h \
  libktrax.la \
  Globals.o \
  PoleTree.o \
  PBCFarField.o \
  BlokTrPdJ.o

JGen.o : \
  Globals.o \
  PoleTree.o \
  TreeWalk.o \
  PBCFarField.o \
  PlaneWise.o

MakePFFT.o : \
  libktrax.la \
  PFFTensor.o \
  DivPFFTensor.o

NukE.o : \
  Globals.o \
  PoleTree.o \
  TreeWalk.o \
  PBCFarField.o

PBCFarField.o : \
  PFFTensor.o \
  PoleTree.o \
  MondoPoles.o \
  Globals.o \
  PFFTensor.o

PFFTensor.o : \
  MondoPoles.o \
  PlaneWise.o

ParaQ.o : \
  Globals.o \
  PoleTree.o

PoleTree.o : \
  Thresholds.o \
  MondoPoles.o \
  Globals.o \
  Density.o

QCTC.o : \
  $(top_builddir)/MondoConfig.h \
  libktrax.la \
  PoleTree.o \
  Globals.o \
  PBCFarField.o \
  JGen.o \
  NukE.o

TreeWalk.o : \
  Globals.o \
  PoleTree.o

Density.o : \
  Thresholds.o \
  Globals.o

Thresholds.o : \
  PoleNodeType.o

Indexing.o : \
  PoleNodeType.o

MondoPoles.o : \
  Indexing.o \
  Globals.o

PlaneWise.o : \
  MondoPoles.o
