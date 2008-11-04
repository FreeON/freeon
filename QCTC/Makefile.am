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
  DEI.f PlaneWise.f90 \
  PFFTensor.f90 \
  DivPFFTensor.f90 \
  ParaQ.f90 \
  MakePFFT.f90


QCTC_DEPENDENCIES     = $(MondoMods_DEPENDENCY_FILES)
JForce_DEPENDENCIES   = $(MondoMods_DEPENDENCY_FILES)
MakePFFT_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES)

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
  PBCTensor/Majik_Kubic_WS1.Inc \
  PBCTensor/Majik_Kubic_WS2.Inc

EXTRA_INCLUDES = -I$(top_srcdir)/MondoMods -I$(top_srcdir)/OneE -I$(top_srcdir)/MondoMods/MMA/LookUpTables_800_6x
LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

# Dependencies.

BlokTrPdJ.o: Globals.o PBCFarField.o PoleTree.o Thresholding.o TreeWalk.o
Density.o: Thresholds.o
DivPFFTensor.o: MondoPoles.o PFFTensor.o Thresholding.o
Globals.o:
Indexing.o: PoleNodeType.o
JForce.o: BlokTrPdJ.o Globals.o PBCFarField.o PoleTree.o
JGen.o: Globals.o PBCFarField.o PoleTree.o Thresholds.o TreeWalk.o
MakePFFT.o: DivPFFTensor.o PFFTensor.o
MondoPoles.o: Globals.o Indexing.o PoleNodeType.o Thresholding.o
NukE.o: Globals.o PBCFarField.o PoleTree.o Thresholding.o Thresholds.o \
	TreeWalk.o
ParaQ.o: Globals.o PoleNodeType.o PoleTree.o
PBCFarField.o: Globals.o MondoPoles.o PoleTree.o Thresholding.o
PFFTensor.o: MondoPoles.o PlaneWise.o Thresholding.o
PlaneWise.o: MondoPoles.o
PoleNodeType.o:
PoleTree.o: Density.o Globals.o Indexing.o MondoPoles.o PoleNodeType.o \
	Thresholds.o
QCTC.o: Density.o Globals.o JGen.o NukE.o PBCFarField.o PoleTree.o \
	Thresholds.o TreeWalk.o
Thresholding.o:
Thresholds.o: PoleNodeType.o Thresholding.o
TreeWalk.o: Globals.o Indexing.o MondoPoles.o PoleNodeType.o PoleTree.o \
	Thresholding.o Thresholds.o
