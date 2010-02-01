include $(top_srcdir)/make.inc

if INTERNAL_LAPACK
INTERNAL_LAPACK_DIRS = lapack
endif

if INTERNAL_HDF5
INTERNAL_HDF5_DIRS = hdf5-1.8.3
HDF5_INCLUDES = -I$(srcdir)/hdf5-1.8.3/src -I$(builddir)/hdf5-1.8.3/src
endif

DIST_SUBDIRS = lapack hdf5-1.8.3
SUBDIRS      = $(INTERNAL_HDF5_DIRS) $(INTERNAL_LAPACK_DIRS)

if USE_DYNAMIC_LINKING
lib_LTLIBRARIES = libfreeonmodules.la
else
noinst_LTLIBRARIES = libfreeonmodules.la
endif

libfreeonmodules_la_SOURCES = \
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
  FreeONSleep.c \
  Functionals.f90 \
  GetAddress.c \
  GetHostname.c \
  GetLimits.c \
  GetPID.c \
  GetPWD.c \
  GetRUsage.c \
  GlobalCharacters.f90 \
  GlobalObjects.f90 \
  GlobalScalars.f90 \
  HDF5Wrap.c \
  InOut.f90 \
  Indexing.f90 \
  InitDS.f \
  LinAlg.f90 \
  LinAlg77.f \
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
  TempDirectory.c \
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
  MMA/LookUpTables_800_6x/F0_77.Inc \
  MMA/LookUpTables_800_6x/F0.Inc \
  MMA/LookUpTables_800_6x/F10_77.Inc \
  MMA/LookUpTables_800_6x/F10.Inc \
  MMA/LookUpTables_800_6x/F11_77.Inc \
  MMA/LookUpTables_800_6x/F11.Inc \
  MMA/LookUpTables_800_6x/F12_77.Inc \
  MMA/LookUpTables_800_6x/F12.Inc \
  MMA/LookUpTables_800_6x/F13_77.Inc \
  MMA/LookUpTables_800_6x/F13.Inc \
  MMA/LookUpTables_800_6x/F14_77.Inc \
  MMA/LookUpTables_800_6x/F14.Inc \
  MMA/LookUpTables_800_6x/F15_77.Inc \
  MMA/LookUpTables_800_6x/F15.Inc \
  MMA/LookUpTables_800_6x/F1_77.Inc \
  MMA/LookUpTables_800_6x/F1.Inc \
  MMA/LookUpTables_800_6x/F2_77.Inc \
  MMA/LookUpTables_800_6x/F2.Inc \
  MMA/LookUpTables_800_6x/F3_77.Inc \
  MMA/LookUpTables_800_6x/F3.Inc \
  MMA/LookUpTables_800_6x/F4_77.Inc \
  MMA/LookUpTables_800_6x/F4.Inc \
  MMA/LookUpTables_800_6x/F5_77.Inc \
  MMA/LookUpTables_800_6x/F5.Inc \
  MMA/LookUpTables_800_6x/F6_77.Inc \
  MMA/LookUpTables_800_6x/F6.Inc \
  MMA/LookUpTables_800_6x/F7_77.Inc \
  MMA/LookUpTables_800_6x/F7.Inc \
  MMA/LookUpTables_800_6x/F8_77.Inc \
  MMA/LookUpTables_800_6x/F8.Inc \
  MMA/LookUpTables_800_6x/F9_77.Inc \
  MMA/LookUpTables_800_6x/F9.Inc \
  MMA/LookUpTables_800_6x/GammaAssymptotics.Inc \
  MMA/LookUpTables_800_6x/GammaDimensions_77.Inc \
  MMA/LookUpTables_800_6x/GammaDimensions.Inc \
  MMA/LookUpTables_800_6x/GammaGrid_77.Inc \
  MMA/LookUpTables_800_6x/GammaGrid.Inc \
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
  $(HDF5_INCLUDES) \
  -I$(srcdir)/MMA/LookUpTables \
  -I$(srcdir)/MMA/Functionals  \
  -I$(srcdir)/MMA/PAC

SOURCES_DEPCHECK = $(libfreeonmodules_la_SOURCES)

# In case we have ctags installed, make the necessary tags files.
.PHONY : tags
tags :
	ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.{c,f,f90}

# Dependencies.
InOut.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  MemMan.lo \
  Indexing.lo \
  Parse.lo \
  SetEq.lo \
  MondoLogger.lo \
  Utilities.lo \
  MondoMPI.lo \
  $(top_builddir)/MondoConfig.h

DerivedTypes.lo : \
  GlobalScalars.lo \
  GlobalCharacters.lo

MemMan.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  ProcessControl.lo \
  MondoLogger.lo \
  $(top_builddir)/MondoConfig.h

Parse.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  ParsingConstants.lo \
  ProcessControl.lo \
  MemMan.lo \
  MondoLogger.lo

Clock.lo : \
  DerivedTypes.lo \
  $(top_builddir)/MondoConfig.h

GlobalObjects.lo : \
  DerivedTypes.lo \
  $(top_builddir)/MondoConfig.h

BasisSetParameters.lo : GlobalCharacters.lo

MondoMPI.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  ProcessControl.lo \
  MemMan.lo \
  Parse.lo \
  MondoLogger.lo \
  Utilities.lo \
  $(top_builddir)/MondoConfig.h

SetEq.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  MemMan.lo \
  MondoMPI.lo \
  $(top_builddir)/MondoConfig.h

Order.lo : \
  DerivedTypes.lo \
  GlobalCharacters.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  ProcessControl.lo \
  MemMan.lo \
  Anneal.lo \
  Parse.lo

Anneal.lo : MemMan.lo PrettyPrint.lo

Indexing.lo : DerivedTypes.lo

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
  MondoLogger.lo \
  MondoMPI.lo \
  $(top_builddir)/MondoConfig.h

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
  Thresholding.lo \
  MondoMPI.lo \
  $(top_builddir)/MondoConfig.h

GlobalCharacters.lo : GlobalScalars.lo

McMurch.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  MemMan.lo

DerivedTypes.lo : \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  $(top_builddir)/MondoConfig.h

MatFunk.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  PrettyPrint.lo

Thresholding.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  InOut.lo \
  SpecFun.lo \
  Parse.lo \
  McMurch.lo

SpecFun.lo : \
  DerivedTypes.lo \
  Parse.lo

AtomPairs.lo : DerivedTypes.lo

PBC.lo : \
  DerivedTypes.lo \
  MondoLogger.lo \
  LinAlg.lo

Functional.lo : DerivedTypes.lo

BraBloks.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  PrettyPrint.lo \
  McMurch.lo \
  Thresholding.lo

Mechanics.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  InOut.lo

Slatec.lo : GlobalScalars.lo ProcessControl.lo

AInv.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  MemMan.lo \
  LinAlg.lo \
  MatFunk.lo \
  $(top_builddir)/MondoConfig.h

ParsingConstants.lo : GlobalScalars.lo

BoundingBox.lo : DerivedTypes.lo

Functionals.lo : \
  DerivedTypes.lo \
  ProcessControl.lo \
  Parse.lo

FastMat.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  MemMan.lo \
  Order.lo \
  LinAlg.lo \
  Thresholding.lo \
  MondoMPI.lo \
  $(top_builddir)/MondoConfig.h

CholFact.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalObjects.lo \
  GlobalCharacters.lo \
  InOut.lo \
  MemMan.lo \
  PrettyPrint.lo \
  ParsingConstants.lo

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
  Mechanics.lo \
  MondoMPI.lo \
  MondoLogger.lo \
  $(top_builddir)/MondoConfig.h

PBC.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  PrettyPrint.lo \
  Thresholding.lo \
  Order.lo \
  MemMan.lo

AtomPairs.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  GlobalObjects.lo \
  BoundingBox.lo \
  Thresholding.lo \
  PBC.lo \
  LinAlg.lo

EllToInclude.lo : DerivedTypes.lo GlobalScalars.lo

FTest.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  InOut.lo \
  PrettyPrint.lo \
  Parse.lo \
  MemMan.lo \
  SpecFun.lo

SMTest.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  InOut.lo \
  PrettyPrint.lo \
  Parse.lo \
  MemMan.lo \
  Order.lo

ProcessControl.lo : \
  DerivedTypes.lo \
  GlobalScalars.lo \
  GlobalCharacters.lo \
  ParsingConstants.lo \
  MondoLogger.lo \
  $(top_builddir)/MondoConfig.h

MondoLogger.lo : \
  GlobalObjects.lo \
  GlobalCharacters.lo \
  ParsingConstants.lo

Utilities.lo : \
  MondoLogger.lo \
  ParsingConstants.lo \
  Parse.lo \
  GlobalScalars.lo \
  FreeONSleep.lo
