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
  CTraX10.F \
  CTraX12.F \
  CTraX14.F \
  CTraX4.F \
  CTraX6.F \
  CTraX8.F \
  CTraX.F \
  Density.F90 \
  Gamma.F \
  Indexing.F90 \
  MondoPoles77.F \
  MondoPoles.F90 \
  PoleNodeType.F90 \
  PoleTree.F90 \
  TreeWalk.F90 \
  XLate.F

bin_PROGRAMS = QCTC JForce MakePFFT

QCTC_SOURCES = \
  Globals.F90 \
  JGen.F90 \
  MD3TRR.F \
  NukE.F90 \
  ParaQ.F90 \
  PBCFarField.F90 \
  QCTC.F90 \
  Thresholds.F90

JForce_SOURCES = \
  BlokTrPdJ.F90 \
  Globals.F90 \
  JForce.F90 \
  MD3TRR.F \
  ParaQ.F90 \
  PBCFarField.F90 \
  Thresholds.F90

MakePFFT_SOURCES = \
  DEI.F \
  DivPFFTensor.F90 \
  Globals.F90 \
  MakePFFT.F90 \
  ParaQ.F90 \
  PFFTensor.F90 \
  PlaneWise.F90 \
  QSortI.F \
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
  libqctc.la \
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
  libqctc.la \
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
  libqctc.la \
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
