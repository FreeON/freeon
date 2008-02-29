include ../make.inc

SUBDIRS = LibQCTC

bin_PROGRAMS = QCTC JForce MakePFFT

QCTC_SOURCES = \
  XLate.f \
  CTraX.f \
  Globals.f90 \
  MondoPoles.f90 \
  Thresholds.f90 \
  PoleTree.f90 \
  PFFTensor.f90 \
  DivPFFTensor.f90 \
  ParaQ.f90 \
  PBCFarField.f90 \
  TreeWalk.f90 \
  JGen.f90 \
  NukE.f90 \
  QCTC.f90

JForce_SOURCES = \
  XLate.f \
  CTraX.f \
  Globals.f90 \
  MondoPoles.f90 \
  Thresholds.f90 \
  PoleTree.f90 \
  PFFTensor.f90 \
  DivPFFTensor.f90 \
  ParaQ.f90 \
  PBCFarField.f90 \
  TreeWalk.f90 \
  BlokTrPdJ.f90 \
	JForce.f90

MakePFFT_SOURCES = \
  XLate.f \
  CTraX.f \
  Globals.f90 \
  MondoPoles.f90 \
  Thresholds.f90 \
  PoleTree.f90 \
  PFFTensor.f90 \
  DivPFFTensor.f90 \
  ParaQ.f90 \
  PBCFarField.f90 \
  TreeWalk.f90 \
	MakePFFT.f90

LibQCTC_DEPENDENCY_FILES = \
  LibQCTC/libQCTC.a \
  LibQCTC/eriglobals.mod \
  LibQCTC/polenodetype.mod \
  LibQCTC/poleglobals.mod

QCTC_DEPENDENCIES     = $(MondoMods_DEPENDENCY_FILES) $(LibQCTC_DEPENDENCY_FILES)
JForce_DEPENDENCIES   = $(MondoMods_DEPENDENCY_FILES) $(LibQCTC_DEPENDENCY_FILES)
MakePFFT_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES) $(LibQCTC_DEPENDENCY_FILES)

SOURCES_DEPCHECK = \
  $(QCTC_SOURCES) \
  $(JForce_SOURCES) \
  $(MakePFFT_SOURCES)

EXTRA_DIST = \
  AuxInt.Inc \
  CTraX.Inc \
  HGToSP.Inc \
  HGTraX.Inc \
  IrRegulars.Inc \
  PBCTensor/Majik_Kubic_WS1.Inc \
  PBCTensor/Majik_Kubic_WS2.Inc

EXTRA_INCLUDES = -I../MondoMods -I../OneE -ILibQCTC
LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5) LibQCTC/libQCTC.a
CPPMISC = -DEXPLICIT_SOURCE

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
	BlokTrPdJ.o

JGen.o : \
  Globals.o \
	PoleTree.o \
	TreeWalk.o \
	PBCFarField.o

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
  MondoPoles.o

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
	PBCFarField.o \
	JGen.o \
	NukE.o

TreeWalk.o : \
  Globals.o \
  PoleTree.o
