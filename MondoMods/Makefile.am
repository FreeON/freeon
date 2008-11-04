include $(top_srcdir)/make.inc

if INTERNAL_LAPACK
INTERNAL_LAPACK_DIRS = lapack
endif

if INTERNAL_HDF5
INTERNAL_HDF5_DIRS = hdf5
endif

SUBDIRS = $(INTERNAL_HDF5_DIRS) $(INTERNAL_LAPACK_DIRS)

if USE_DYNAMIC_LINKING
lib_LTLIBRARIES = libMondoMods.la
else
noinst_LTLIBRARIES = libMondoMods.la
endif

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

AInv.lo: DerivedTypes.lo GlobalCharacters.lo GlobalScalars.lo LinAlg.lo MatFunk.lo \
	MemMan.lo
Anneal.lo: MemMan.lo PrettyPrint.lo
AtomPairs.lo: BoundingBox.lo DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo \
	GlobalScalars.lo LinAlg.lo PBC.lo Thresholding.lo
BasisSetParameters.lo: GlobalCharacters.lo
BoundingBox.lo: DerivedTypes.lo
BraBloks.lo: DerivedTypes.lo GlobalScalars.lo McMurch.lo Mechanics.lo \
	PrettyPrint.lo Thresholding.lo
CholFact.lo: DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo GlobalScalars.lo \
	InOut.lo MemMan.lo MondoLogger.lo ParsingConstants.lo PrettyPrint.lo
Clock.lo: DerivedTypes.lo MondoMPI.lo
CWrappers.lo: MondoLogger.lo ParsingConstants.lo
DerivedTypes.lo: GlobalCharacters.lo GlobalScalars.lo
FastMat.lo: DerivedTypes.lo GlobalObjects.lo GlobalScalars.lo LinAlg.lo MemMan.lo \
	MondoMPI.lo Order.lo Thresholding.lo
Functionals.lo: DerivedTypes.lo Parse.lo ProcessControl.lo
GlobalCharacters.lo: GlobalScalars.lo
GlobalObjects.lo: DerivedTypes.lo
Indexing.lo: DerivedTypes.lo
InOut.lo: DerivedTypes.lo GlobalCharacters.lo GlobalScalars.lo Indexing.lo \
	MemMan.lo MondoMPI.lo Parse.lo ProcessControl.lo SetEq.lo
LinAlg.lo: Clock.lo DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo \
	GlobalScalars.lo InOut.lo LinAlg.lo MemMan.lo MondoMPI.lo Order.lo \
	PrettyPrint.lo ProcessControl.lo SetEq.lo Thresholding.lo
Macros.lo: AtomPairs.lo Clock.lo DerivedTypes.lo Functionals.lo GlobalCharacters.lo \
	GlobalObjects.lo GlobalScalars.lo InOut.lo Mechanics.lo MondoMPI.lo \
	Parse.lo ParsingConstants.lo PrettyPrint.lo
MatFunk.lo: DerivedTypes.lo GlobalCharacters.lo GlobalScalars.lo PrettyPrint.lo
McMurch.lo: DerivedTypes.lo GlobalScalars.lo MemMan.lo
Mechanics.lo: DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo \
	GlobalScalars.lo InOut.lo
MemMan.lo: DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo GlobalScalars.lo \
	ProcessControl.lo
MondoLogger.lo: GlobalCharacters.lo GlobalObjects.lo ParsingConstants.lo
MondoMPI.lo: DerivedTypes.lo GlobalObjects.lo GlobalScalars.lo MPIInclude.lo \
	MemMan.lo Parse.lo ProcessControl.lo
Order.lo: Anneal.lo DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo \
	GlobalScalars.lo MemMan.lo Parse.lo ProcessControl.lo
Parse.lo: DerivedTypes.lo GlobalCharacters.lo GlobalScalars.lo MemMan.lo \
	ParsingConstants.lo ProcessControl.lo
ParsingConstants.lo: GlobalScalars.lo
PBC.lo: DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo GlobalScalars.lo \
	LinAlg.lo MemMan.lo Order.lo PrettyPrint.lo Thresholding.lo
PrettyPrint.lo: BasisSetParameters.lo Clock.lo DerivedTypes.lo GlobalCharacters.lo \
	GlobalObjects.lo GlobalScalars.lo InOut.lo MondoMPI.lo Parse.lo \
	ParsingConstants.lo SetEq.lo
ProcessControl.lo: DerivedTypes.lo GlobalCharacters.lo GlobalScalars.lo \
	MPIInclude.lo MondoLogger.lo ParsingConstants.lo
SetEq.lo: DerivedTypes.lo GlobalObjects.lo GlobalScalars.lo MemMan.lo MondoMPI.lo
Slatec.lo: GlobalScalars.lo ProcessControl.lo
SpecFun.lo: DerivedTypes.lo Parse.lo ProcessControl.lo SpecFun.lo
Thresholding.lo: DerivedTypes.lo GlobalCharacters.lo GlobalObjects.lo \
	GlobalScalars.lo InOut.lo McMurch.lo Mechanics.lo Parse.lo SpecFun.lo

