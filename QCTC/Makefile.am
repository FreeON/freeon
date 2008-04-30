include $(top_srcdir)/make.inc

bin_PROGRAMS = QCTC JForce MakePFFT

QCTC_SOURCES = \
PoleNodeType.f90 Indexing.f90 XLate.f  CTraX.f \
CTraX4.f CTraX6.f CTraX8.f CTraX10.f  CTraX12.f  CTraX14.f  \
Globals.f90 Gamma.f KTrax.f MD3TRR.f AA95.f AA95c.f   \
MondoPoles77.f MondoPoles.f90 Thresholds.f90    \
Density.f90 PoleTree.f90 TreeWalk.f90 PBCFarField.f90   ParaQ.f90 \
  JGen.f90 \
  NukE.f90 \
  QCTC.f90

JForce_SOURCES = \
PoleNodeType.f90 Indexing.f90 XLate.f  CTraX.f \
CTraX4.f CTraX6.f CTraX8.f CTraX10.f  CTraX12.f  CTraX14.f  \
Globals.f90 Gamma.f KTrax.f MD3TRR.f AA95.f AA95c.f   \
MondoPoles77.f MondoPoles.f90 Thresholds.f90    \
Density.f90 PoleTree.f90 TreeWalk.f90 PBCFarField.f90  ParaQ.f90 \
  BlokTrPdJ.f90 \
  JForce.f90


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
CTraX.Inc  HGToSP_Density.Inc  HGToSP.Inc  HGToSP_PoleNode.Inc MACErrBnd4.Inc  \
  PBCTensor/Majik_Kubic_WS1.Inc \
  PBCTensor/Majik_Kubic_WS2.Inc

EXTRA_INCLUDES = -I../MondoMods -I../OneE -I../MondoMods/MMA/LookUpTables_800_6x
LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5) 

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
	ParaQ.o\
	Globals.o \
	PoleTree.o \
  PBCFarField.o \
	BlokTrPdJ.o

JGen.o : \
  Globals.o \
	PoleTree.o \
	TreeWalk.o \
	PBCFarField.o

MakePFFT.o : \
  PFFTensor.o \
  PlaneWise.o \
  DEI.o 

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
	Globals.o

QCTC.o : \
  PoleTree.o \
	Globals.o \
	ParaQ.o \
	PBCFarField.o \
	JGen.o \
	NukE.o

TreeWalk.o : \
  Globals.o \
  PoleTree.o
