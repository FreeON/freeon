include $(top_srcdir)/make.inc

if INTERNAL_LAPACK
  INTERNAL_LAPACK_DIRS = lapack
endif

if INTERNAL_HDF5
  INTERNAL_HDF5_DIRS = hdf5
endif

SUBDIRS = $(INTERNAL_HDF5_DIRS) $(INTERNAL_LAPACK_DIRS)

noinst_LTLIBRARIES = libMondoMods.la

libMondoMods_la_SOURCES = \
  AInv.f90 \
  Anneal.f90 \
  AtomPairs.f90 \
  BasisSetParameters.f90 \
  BoundingBox.f90 \
  BraBloks.f90 \
  CholFact.f90 \
  Clock.f90 \
  D1Mach.f \
  DCSEvl.f \
  DErf.f \
  DErfc.f \
  DFRepository.f \
  DMach.f \
  DerivedTypes.f90 \
  FastMat.f90 \
  FileCopy.c \
  CWrappers.f90 \
  Functionals.f90 \
  GetAddress.c \
  GetHostname.c \
  GetLimits.c \
  GetRUsage.c \
  GetPWD.c \
  GlobalCharacters.f90 \
  GlobalObjects.f90 \
  GlobalScalars.f90 \
  HDF5Wrap.c \
  InOut.f90 \
  Indexing.f90 \
  InitDS.f \
  LinAlg.f90 \
  LinAlg77.f \
  MPIInclude.f90 \
  Macros.f90 \
  MatFunk.f90 \
  McMurch.f90 \
  Mechanics.f90 \
  MemMan.f90 \
  MondoLogger.f90 \
  MondoLogger.h \
  MondoMPI.f90 \
  Order.f90 \
  Order77.f \
  PBC.f90 \
  Parse.f90 \
  ParsingConstants.f90 \
  PrettyPrint.f90 \
  ProcessControl.f90 \
  Random.c \
  SMWrap.c \
  Seconds.c \
  SetEq.f90 \
  SetEq77.f \
  Slatec.f90 \
  SpecFun.f90 \
  Thresholding.f90 \
  TimerWrapper.c \
  Trap.c \
  UnlimitStack.c \
  Utilities.f90 \
  mm_double_NN_1.c \
  mm_double_NN_c.c \
  mm_double_NT_1.c \
  mm_double_NT_c.c \
  mm_double_TN_1.c \
  mm_double_TN_c.c \
  mm_double_TT_1.c \
  mm_double_TT_c.c \
  phipac_dgemm.c

EXTRA_DIST = \
  Final_TravO.Inc \
  MMA/LookUpTables/Erf.Inc \
  MMA/LookUpTables/ErrorFunction.Inc \
  MMA/LookUpTables/Exp.Inc \
  MMA/LookUpTables/F0.Inc \
  MMA/LookUpTables/F1.Inc \
  MMA/LookUpTables/F10.Inc \
  MMA/LookUpTables/F11.Inc \
  MMA/LookUpTables/F12.Inc \
  MMA/LookUpTables/F13.Inc \
  MMA/LookUpTables/F14.Inc \
  MMA/LookUpTables/F15.Inc \
  MMA/LookUpTables/F2.Inc \
  MMA/LookUpTables/F3.Inc \
  MMA/LookUpTables/F4.Inc \
  MMA/LookUpTables/F5.Inc \
  MMA/LookUpTables/F6.Inc \
  MMA/LookUpTables/F7.Inc \
  MMA/LookUpTables/F8.Inc \
  MMA/LookUpTables/F9.Inc \
  MMA/LookUpTables/GammaAssymptotics.Inc \
  MMA/LookUpTables/GammaDimensions.Inc \
  MMA/LookUpTables/GammaFunctions.Inc \
  MMA/LookUpTables/GammaGrid.Inc \
  MMA/LookUpTables/InverseExp.Inc \
  MMA/Functionals/B88x.Inc \
  MMA/Functionals/LYPc.Inc \
  MMA/Functionals/PBEc.Inc \
  MMA/Functionals/PBEx.Inc \
  MMA/Functionals/PW91c.Inc \
  MMA/Functionals/PW91x.Inc \
  MMA/Functionals/SDx.Inc \
  MMA/Functionals/VWN3c.Inc \
  MMA/Functionals/VWN5c.Inc \
  MMA/Functionals/XAx.Inc \
  MMA/Functionals/XG04x.Inc \
  MMA/PAC/PFunk.Inc

EXTRA_INCLUDES = \
  -IMMA/LookUpTables \
  -IMMA/Functionals  \
  -IMMA/PAC

SOURCES_DEPCHECK = $(libMondoMods_la_SOURCES)

# Dependencies.
InOut.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  MemMan.o \
  Indexing.o \
  Parse.o \
  SetEq.o \
  MondoLogger.o

InOut.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  MemMan.lo \
  Indexing.lo \
  Parse.lo \
  SetEq.lo \
  MondoLogger.lo

DerivedTypes.o : \
  GlobalScalars.o \
  GlobalCharacters.o

DerivedTypes.lo : \
  GlobalScalars.lo \
  GlobalCharacters.lo

MemMan.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  ProcessControl.o \
  MondoLogger.o

MemMan.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  ProcessControl.lo \
  MondoLogger.lo

Parse.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  ParsingConstants.o \
  ProcessControl.o \
  MemMan.o

Parse.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  ParsingConstants.lo \
  ProcessControl.lo \
  MemMan.lo

Clock.o : DerivedTypes.o
Clock.lo : DerivedTypes.lo

GlobalObjects.o : DerivedTypes.o
GlobalObjects.lo : DerivedTypes.lo

BasisSetParameters.o : GlobalCharacters.o
BasisSetParameters.lo : GlobalCharacters.lo

MondoMPI.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  ProcessControl.o \
  MPIInclude.o \
  MemMan.o \
  Parse.o

MondoMPI.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  ProcessControl.lo \
  MPIInclude.lo \
  MemMan.lo \
  Parse.lo

SetEq.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  MemMan.o

SetEq.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  MemMan.lo

Order.o : \
  DerivedTypes.o \
  GlobalCharacters.o \
  GlobalScalars.o \
  GlobalObjects.o \
  ProcessControl.o \
  MemMan.o \
  Anneal.o \
  Parse.o

Order.lo : \
  DerivedTypes.lo \
  GlobalCharacters.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  ProcessControl.lo \
  MemMan.lo \
  Anneal.lo \
  Parse.lo

Anneal.o : MemMan.o PrettyPrint.o
Anneal.lo : MemMan.lo PrettyPrint.lo

Indexing.o : DerivedTypes.o
Indexing.lo : DerivedTypes.lo

PrettyPrint.o : \
  DerivedTypes.o \
  GlobalCharacters.o \
  GlobalScalars.o \
  GlobalObjects.o \
  BasisSetParameters.o \
  Clock.o \
  Parse.o \
  SetEq.o \
  InOut.o \
  MondoLogger.o

PrettyPrint.lo : \
  DerivedTypes.lo \
  GlobalCharacters.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  BasisSetParameters.lo \
  Clock.lo \
  Parse.lo \
  SetEq.lo \
  InOut.lo \
  MondoLogger.lo

LinAlg.o : \
  GlobalCharacters.o \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  ProcessControl.o \
  PrettyPrint.o \
  MemMan.o \
  SetEq.o \
  InOut.o \
  Order.o \
  Clock.o \
  Thresholding.o

LinAlg.lo : \
  GlobalCharacters.lo \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  ProcessControl.lo \
  PrettyPrint.lo \
  MemMan.lo \
  SetEq.lo \
  InOut.lo \
  Order.lo \
  Clock.lo \
  Thresholding.lo

GlobalCharacters.o : GlobalScalars.o
GlobalCharacters.lo : GlobalScalars.lo

McMurch.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  MemMan.o

McMurch.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  MemMan.lo

DerivedTypes.o : \
  GlobalScalars.o \
  GlobalCharacters.o

DerivedTypes.lo : \
  GlobalScalars.lo \
  GlobalCharacters.lo

MatFunk.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  PrettyPrint.o

MatFunk.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  PrettyPrint.lo

Thresholding.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  InOut.o \
  SpecFun.o \
  Parse.o \
  McMurch.o

Thresholding.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  InOut.lo \
  SpecFun.lo \
  Parse.lo \
  McMurch.lo

SpecFun.o : \
  DerivedTypes.o \
	Parse.o

SpecFun.lo : \
  DerivedTypes.lo \
	Parse.lo

AtomPairs.o : DerivedTypes.o
AtomPairs.lo : DerivedTypes.lo

PBC.o : DerivedTypes.o MondoLogger.o
PBC.lo : DerivedTypes.lo MondoLogger.lo

Functional.o : DerivedTypes.o
Functional.lo : DerivedTypes.lo

BraBloks.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  PrettyPrint.o \
  McMurch.o \
  Thresholding.o

BraBloks.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  PrettyPrint.lo \
  McMurch.lo \
  Thresholding.lo

Mechanics.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
	InOut.o

Mechanics.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
	InOut.lo

Slatec.o : GlobalScalars.o ProcessControl.o
Slatec.lo : GlobalScalars.lo ProcessControl.lo

AInv.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  MemMan.o \
  LinAlg.o \
  MatFunk.o

AInv.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  MemMan.lo \
  LinAlg.lo \
  MatFunk.lo

ParsingConstants.o : GlobalScalars.o
ParsingConstants.lo : GlobalScalars.lo

BoundingBox.o : DerivedTypes.o
BoundingBox.lo : DerivedTypes.lo

Functionals.o : \
  DerivedTypes.o \
  ProcessControl.o \
  Parse.o

Functionals.lo : \
  DerivedTypes.lo \
  ProcessControl.lo \
  Parse.lo

FastMat.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  MemMan.o \
  Order.o \
  LinAlg.o \
  Thresholding.o

FastMat.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  MemMan.lo \
  Order.lo \
  LinAlg.lo \
  Thresholding.lo

CholFact.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  GlobalCharacters.o \
  InOut.o \
  MemMan.o \
  PrettyPrint.o \
  ParsingConstants.o

CholFact.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  GlobalCharacters.lo \
  InOut.lo \
  MemMan.lo \
  PrettyPrint.lo \
  ParsingConstants.lo

Macros.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  GlobalCharacters.o \
  ParsingConstants.o \
  PrettyPrint.o \
  Parse.o \
  Clock.o \
  InOut.o \
  Functionals.o \
  AtomPairs.o \
  Mechanics.o

Macros.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  GlobalCharacters.lo \
  ParsingConstants.lo \
  PrettyPrint.lo \
  Parse.lo \
  Clock.lo \
  InOut.lo \
  Functionals.lo \
  AtomPairs.lo \
  Mechanics.lo

PBC.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  PrettyPrint.o \
  Thresholding.o \
  Order.o \
  MemMan.o

PBC.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  PrettyPrint.lo \
  Thresholding.lo \
  Order.lo \
  MemMan.lo

AtomPairs.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  BoundingBox.o \
  Thresholding.o \
  PBC.o \
  LinAlg.o

AtomPairs.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  BoundingBox.lo \
  Thresholding.lo \
  PBC.lo \
  LinAlg.lo

EllToInclude.o : DerivedTypes.o GlobalScalars.o
EllToInclude.lo : DerivedTypes.lo GlobalScalars.lo

FTest.o : \
  DerivedTypes.o \
	GlobalScalars.o \
	GlobalCharacters.o \
	InOut.o \
	PrettyPrint.o \
	Parse.o \
	MemMan.o \
	SpecFun.o

FTest.lo : \
  DerivedTypes.lo \
	GlobalScalars.lo \
	GlobalCharacters.lo \
	InOut.lo \
	PrettyPrint.lo \
	Parse.lo \
	MemMan.lo \
	SpecFun.lo

SMTest.o : \
  DerivedTypes.o \
	GlobalScalars.o \
	GlobalCharacters.o \
	InOut.o \
	PrettyPrint.o \
	Parse.o \
	MemMan.o \
	Order.o

SMTest.lo : \
  DerivedTypes.lo \
	GlobalScalars.lo \
	GlobalCharacters.lo \
	InOut.lo \
	PrettyPrint.lo \
	Parse.lo \
	MemMan.lo \
	Order.lo

ProcessControl.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  ParsingConstants.o \
  MondoLogger.o

ProcessControl.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  ParsingConstants.lo \
  MondoLogger.lo

MondoLogger.o : \
  GlobalObjects.o \
  GlobalCharacters.o \
  ParsingConstants.o

MondoLogger.lo : \
  GlobalObjects.lo \
  GlobalCharacters.lo \
  ParsingConstants.lo

CWrappers.o : \
  ParsingConstants.o \
  MondoLogger.o

CWrappers.lo : \
  ParsingConstants.lo \
  MondoLogger.lo
