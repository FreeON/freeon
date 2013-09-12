include $(top_srcdir)/make.inc

bin_PROGRAMS = QCTC JForce MakePFFT

QCTC_SOURCES = \
  AA95c.f \
  AA95.f \
  CTraX10.f \
  CTraX12.f \
  CTraX14.f \
  CTraX4.f \
  CTraX6.f \
  CTraX8.f \
  CTraX.f \
  Density.f90 \
  Gamma.f \
  Globals.f90 \
  Indexing.f90 \
  JGen.f90 \
  KTrax.f \
  MD3TRR.f \
  MondoPoles77.f \
  MondoPoles.f90 \
  NukE.f90 \
  ParaQ.f90 \
  PBCFarField.f90 \
  PoleNodeType.f90 \
  PoleTree.f90 \
  QCTC.f90 \
  Thresholds.f90 \
  TreeWalk.f90 \
  XLate.f

JForce_SOURCES = \
  AA95c.f \
  AA95.f \
  BlokTrPdJ.f90 \
  CTraX10.f \
  CTraX12.f \
  CTraX14.f \
  CTraX4.f \
  CTraX6.f \
  CTraX8.f \
  CTraX.f \
  Density.f90 \
  Gamma.f \
  Globals.f90 \
  Indexing.f90 \
  JForce.f90 \
  KTrax.f \
  MD3TRR.f \
  MondoPoles77.f \
  MondoPoles.f90 \
  ParaQ.f90 \
  PBCFarField.f90 \
  PoleNodeType.f90 \
  PoleTree.f90 \
  Thresholds.f90 \
  TreeWalk.f90 \
  XLate.f

MakePFFT_SOURCES = \
  XLate.f \
  CTraX.f \
  Globals.f90 \
  PoleNodeType.f90 \
  Thresholds.f90 \
  Indexing.f90 \
  QSortI.f \
  MondoPoles.f90 \
  DEI.f \
  PlaneWise.f90 \
  PFFTensor.f90 \
  DivPFFTensor.f90 \
  ParaQ.f90 \
  MakePFFT.f90


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
LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

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
  Globals.o \
  PoleTree.o \
  PBCFarField.o \
  BlokTrPdJ.o \
  $(top_builddir)/MondoConfig.h

JGen.o : \
  Globals.o \
  PoleTree.o \
  TreeWalk.o \
  PBCFarField.o \
  PlaneWise.o

MakePFFT.o : \
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
