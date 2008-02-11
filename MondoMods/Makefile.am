include ../make.inc

if INTERNAL_LAPACK
  INTERNAL_LAPACK_DIRS = lapack
endif

if INTERNAL_HDF5
  INTERNAL_HDF5_DIRS = hdf5
endif

SUBDIRS = $(INTERNAL_LAPACK_DIRS) $(INTERNAL_HDF5_DIRS)

noinst_LIBRARIES = libMondoMods.a

libMondoMods_a_SOURCES = \
  DFRepository.f \
  GlobalScalars.f90 \
  GlobalCharacters.f90 \
  DerivedTypes.f90 \
  GlobalObjects.f90 \
  BasisSetParameters.f90 \
  MPIInclude.f90 \
  ProcessControl.f90 \
  Trap.c \
  MemMan.f90 \
  GetAddress.c \
  Clock.f90 \
  Seconds.c \
  TimerWrapper.c \
  ParsingConstants.f90 \
  Parse.f90 \
  MondoMPI.f90 \
  Indexing.f90 \
  SetEq.f90 \
  SetEq77.f \
  InOut.f90 \
  HDF5Wrap.c \
  Anneal.f90 \
  Order.f90 \
  Order77.f \
  Random.c \
  PrettyPrint.f90 \
  MondoLogger.h \
  MondoLogger.f90 \
  McMurch.f90 \
  LinAlg.f90 \
  LinAlg77.f \
  MatFunk.f90 \
  AInv.f90 \
  GetLimits.c \
  GetRUsage.c \
  SMWrap.c \
  SpecFun.f90 \
  D1Mach.f \
  DCSEvl.f \
  DErf.f \
  DErfc.f \
  DMach.f \
  InitDS.f \
  Mechanics.f90 \
  Thresholding.f90 \
  Functionals.f90 \
  PBC.f90 \
  AtomPairs.f90 \
  BoundingBox.f90 \
  BraBloks.f90 \
  Macros.f90 \
  FastMat.f90 \
  CholFact.f90 \
  Slatec.f90 \
  mm_double_NN_1.c \
  mm_double_NT_1.c \
  mm_double_TN_1.c \
  mm_double_TT_1.c \
  mm_double_NN_c.c \
  mm_double_NT_c.c \
  mm_double_TN_c.c \
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

EXTRA_INCLUDES = -IMMA/LookUpTables \
                 -IMMA/Functionals  \
                 -IMMA/PAC

SOURCES_DEPCHECK = $(libMondoMods_a_SOURCES)

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

DerivedTypes.o : \
  GlobalScalars.o \
  GlobalCharacters.o

MemMan.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  ProcessControl.o

Parse.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  ParsingConstants.o \
  ProcessControl.o \
  MemMan.o

Clock.o : \
  DerivedTypes.o

GlobalObjects.o : \
  DerivedTypes.o

BasisSetParameters.o : \
	GlobalCharacters.o

MondoMPI.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  ProcessControl.o \
  MPIInclude.o \
  MemMan.o \
  Parse.o

SetEq.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  MemMan.o

Order.o : \
  DerivedTypes.o \
  GlobalCharacters.o \
  GlobalScalars.o \
  GlobalObjects.o \
  ProcessControl.o \
  MemMan.o \
  Anneal.o \
  Parse.o

Anneal.o : \
  MemMan.o \
  PrettyPrint.o

Indexing.o : \
  DerivedTypes.o

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

GlobalCharacters.o : \
  GlobalScalars.o

McMurch.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  MemMan.o

DerivedTypes.o : \
  GlobalScalars.o \
  GlobalCharacters.o

MatFunk.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  PrettyPrint.o

Thresholding.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  InOut.o \
  SpecFun.o \
  Parse.o \
  McMurch.o

SpecFun.o : \
  DerivedTypes.o \
	Parse.o

AtomPairs.o : \
  DerivedTypes.o

PBC.o : \
  DerivedTypes.o \
  MondoLogger.o

Functional.o : \
  DerivedTypes.o

BraBloks.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  PrettyPrint.o \
  McMurch.o \
  Thresholding.o

Mechanics.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
	InOut.o

Slatec.o : \
  GlobalScalars.o \
  ProcessControl.o

AInv.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  MemMan.o \
  LinAlg.o \
  MatFunk.o

ParsingConstants.o : \
  GlobalScalars.o

BoundingBox.o : \
  DerivedTypes.o

Functionals.o : \
  DerivedTypes.o \
  ProcessControl.o \
  Parse.o

FastMat.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  MemMan.o \
  Order.o \
  LinAlg.o \
  Thresholding.o

CholFact.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalObjects.o \
  GlobalCharacters.o \
  InOut.o \
  MemMan.o \
  PrettyPrint.o \
  ParsingConstants.o

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

PBC.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  PrettyPrint.o \
  Thresholding.o \
  Order.o \
  MemMan.o

AtomPairs.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  GlobalObjects.o \
  BoundingBox.o \
  Thresholding.o \
  PBC.o \
  LinAlg.o

EllToInclude.o : \
  DerivedTypes.o \
	GlobalScalars.o

FTest.o : \
  DerivedTypes.o \
	GlobalScalars.o \
	GlobalCharacters.o \
	InOut.o \
	PrettyPrint.o \
	Parse.o \
	MemMan.o \
	SpecFun.o

SMTest.o : \
  DerivedTypes.o \
	GlobalScalars.o \
	GlobalCharacters.o \
	InOut.o \
	PrettyPrint.o \
	Parse.o \
	MemMan.o \
	Order.o

ProcessControl.o : \
  DerivedTypes.o \
  GlobalScalars.o \
  GlobalCharacters.o \
  ParsingConstants.o \
  MondoLogger.o

MondoLogger.o : \
  GlobalObjects.o \
  GlobalCharacters.o \
  ParsingConstants.o
