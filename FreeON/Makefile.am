include $(top_srcdir)/make.inc

bin_PROGRAMS = FreeON

FreeON_SOURCES = \
  ConflictDetect.f90 \
  ControlStructures.f90 \
  DynamicsKeys.f90 \
  FreeON.f90 \
  GeometryKeys.f90 \
  GeomOptKeys.f90 \
  GetPID.c \
  Hessian.f90 \
  InCoords.f90 \
  ls_rmsd.f90 \
  Massage.f90 \
  MDynamics.f90 \
  MLP.f90 \
  MonteCarlo.f90 \
  NEB.f90 \
  Numerics.f90 \
  Optimizer.f90 \
  OptionKeys.f90 \
  Overlay.f90 \
  ParallelKeys.f90 \
  ParseBasis.f90 \
  ParseCommands.f90 \
  ParseDynamics.f90 \
  ParseExtraCoords.f90 \
  ParseGeometries.f90 \
  ParseGeomOpt.f90 \
  ParseInput.f90 \
  ParseOptions.f90 \
  ParseParallel.f90 \
  ParsePeriodic.f90 \
  ParseProperties.f90 \
  PeriodicKeys.f90 \
  PrintParsed.f90 \
  PunchHDF.f90 \
  QUICCA.f90 \
  Response.f90 \
  RQI.f90 \
  SCFKeys.f90 \
  SCFs.f90 \
  Spawn.c \
  ZippyQuote.f90

FreeON_DEPENDENCIES = $(LIBRARY_DEPENDENCY_FILES)

SOURCES_DEPCHECK = $(FreeON_SOURCES)

FreeON_LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

CPPMISC = -DNEB2

# Link rule for FreeON.
FreeON_LINK = $(FCLINK)

# Dependencies.
ParseCommands.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/inout.mod \
  ControlStructures.o \
  $(top_builddir)/MondoMods/globalcharacters.mod \
  $(top_builddir)/MondoMods/utilities.mod

ParseOptions.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/memman.mod \
  Numerics.o \
  ConflictDetect.o \
  OptionKeys.o \
  $(top_builddir)/MondoMods/functionals.mod \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  $(top_builddir)/MondoMods/globalobjects.mod \
  $(top_builddir)/MondoMods/processcontrol.mod \
  $(top_builddir)/MondoMods/globalcharacters.mod \
  $(top_builddir)/MondoMods/parsingconstants.mod \
  ControlStructures.o \
  $(top_builddir)/MondoMods/mondologger.mod

ParseGeometries.o : \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/inout.mod \
  OptionKeys.o \
  $(top_builddir)/MondoMods/prettyprint.mod \
  DynamicsKeys.o \
  GeometryKeys.o \
  ControlStructures.o \
  ls_rmsd.o \
  NEB.o \
  ConflictDetect.o \
  $(top_builddir)/MondoMods/mondologger.mod

ParseExtraCoords.o : \
  ControlStructures.o \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  OptionKeys.o \
  GeometryKeys.o \
  ParseGeometries.o \
  Massage.o \
  InCoords.o \
  $(top_builddir)/MondoMods/setxyz.mod \
  $(top_builddir)/MondoMods/memman.mod

ParseGeomOpt.o : \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  GeomOptKeys.o \
  $(top_builddir)/MondoMods/macros.mod \
  ControlStructures.o \
  OptionKeys.o \
  $(top_builddir)/MondoMods/mondologger.mod

ParseProperties.o : \
  $(top_builddir)/MondoMods/inout.mod \
  ControlStructures.o \
  $(top_builddir)/MondoMods/mondologger.mod

ParseInput.o : \
  Massage.o \
  ParseBasis.o \
  ParseOptions.o \
  ParseCommands.o \
  ParseDynamics.o \
  ParsePeriodic.o \
  ParseGeometries.o \
  ControlStructures.o \
  $(top_builddir)/MondoMods/prettyprint.mod \
  ParseParallel.o \
  ParseGeomOpt.o \
  ParseExtraCoords.o \
  ParseProperties.o \
  PrintParsed.o \
  $(top_builddir)/MondoMods/mondologger.mod

ParsePeriodic.o : \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/pbc.mod \
  $(top_builddir)/MondoMods/atompairs.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  OptionKeys.o \
  PeriodicKeys.o \
  ControlStructures.o

ParseParallel.o : \
  $(top_builddir)/MondoMods/inout.mod \
  ParallelKeys.o \
  ControlStructures.o \
  $(top_builddir)/MondoMods/mondologger.mod

ParseBasis.o : \
  $(top_builddir)/MondoMods/order.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  ControlStructures.o \
  $(top_builddir)/MondoMods/basissetparameters.mod \
  $(top_builddir)/MondoMods/mondologger.mod

ParseDynamics.o : \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  DynamicsKeys.o \
  ControlStructures.o \
  OptionKeys.o \
  $(top_builddir)/MondoMods/mondologger.mod

ConflictDetect.o : \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/functionals.mod \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/processcontrol.mod \
  ControlStructures.o

Overlay.o : \
  $(top_builddir)/MondoMods/inout.mod \
  PunchHDF.o \
  ControlStructures.o \
  $(top_builddir)/MondoMods/utilities.mod

MDynamics.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoMods/order.mod \
  SCFs.o \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/matfunk.mod \
  Numerics.o \
  $(top_builddir)/MondoMods/atompairs.mod \
  ControlStructures.o \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/MondoMods/globalscalars.mod \
  $(top_builddir)/MondoMods/setxyz.mod \
  DynamicsKeys.o \
  $(top_builddir)/MondoMods/utilities.mod \
  $(top_builddir)/MondoMods/mondologger.mod

Response.o : \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/linalg.mod \
  $(top_builddir)/MondoMods/globalobjects.mod \
  SCFKeys.o \
  Overlay.o \
  PunchHDF.o \
  Numerics.o \
  OptionKeys.o \
  ControlStructures.o \
  NEB.o \
  SCFs.o \
  $(top_builddir)/MondoMods/setxyz.mod

Hessian.o : \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  $(top_builddir)/MondoMods/globalobjects.mod \
  $(top_builddir)/MondoMods/globalcharacters.mod \
  $(top_builddir)/MondoMods/inout.mod \
  InCoords.o \
  $(top_builddir)/MondoMods/memman.mod \
  $(top_builddir)/MondoMods/setxyz.mod \
  $(top_builddir)/MondoMods/processcontrol.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  $(top_builddir)/MondoMods/parsingconstants.mod \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/MondoMods/linalg.mod \
  $(top_builddir)/MondoMods/ainv.mod \
  $(top_builddir)/MondoMods/cholfactor.mod \
  ControlStructures.o

MonteCarlo.o : \
  $(top_builddir)/MondoMods/order.mod \
  SCFs.o \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/matfunk.mod \
  Numerics.o \
  $(top_builddir)/MondoMods/atompairs.mod \
  ControlStructures.o \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/MondoMods/globalscalars.mod \
  $(top_builddir)/MondoMods/setxyz.mod \
  DynamicsKeys.o \
  MDynamics.o

Optimizer.o : \
  SCFs.o \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/matfunk.mod \
  Numerics.o \
  $(top_builddir)/MondoMods/atompairs.mod \
  ControlStructures.o \
  InCoords.o \
  QUICCA.o \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/MondoMods/globalscalars.mod \
  Hessian.o \
  $(top_builddir)/MondoMods/setxyz.mod \
  $(top_builddir)/MondoMods/utilities.mod \
  $(top_builddir)/MondoMods/globalscalars.mod

FreeON.o : \
  $(top_builddir)/MondoConfig.h \
  SCFs.o \
  $(top_builddir)/MondoMods/macros.mod \
  Response.o \
  PunchHDF.o \
  Optimizer.o \
  ParseInput.o \
  ZippyQuote.o \
  MDynamics.o \
  MonteCarlo.o \
  $(top_builddir)/MondoMods/mondologger.mod \
  RQI.o

QUICCA.o : \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  $(top_builddir)/MondoMods/globalobjects.mod \
  $(top_builddir)/MondoMods/globalcharacters.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/memman.mod \
  $(top_builddir)/MondoMods/linalg.mod \
  $(top_builddir)/MondoMods/ainv.mod \
  $(top_builddir)/MondoMods/cholfactor.mod \
  ControlStructures.o \
  PunchHDF.o \
  $(top_builddir)/MondoMods/setxyz.mod \
  InCoords.o \
  Hessian.o \
  $(top_builddir)/MondoMods/slatec.mod

ZippyQuote.o : \
  $(top_builddir)/MondoMods/order.mod \
  ControlStructures.o

NEB.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoMods/inout.mod \
  ls_rmsd.o \
  $(top_builddir)/MondoMods/prettyprint.mod \
  ControlStructures.o \
  $(top_builddir)/MondoMods/order.mod \
  $(top_builddir)/MondoMods/mondologger.mod

Massage.o : \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/pbc.mod \
  $(top_builddir)/MondoMods/atompairs.mod \
  OptionKeys.o \
  DynamicsKeys.o \
  GeometryKeys.o \
  ControlStructures.o

SCFs.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoMods/parse.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/linalg.mod \
  $(top_builddir)/MondoMods/globalobjects.mod \
  SCFKeys.o \
  Overlay.o \
  DynamicsKeys.o \
  PunchHDF.o \
  Numerics.o \
  OptionKeys.o \
  ControlStructures.o \
  NEB.o \
  $(top_builddir)/MondoMods/setxyz.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  $(top_builddir)/MondoMods/mondologger.mod \
  $(top_builddir)/MondoMods/utilities.mod

InCoords.o : \
  ControlStructures.o \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  $(top_builddir)/MondoMods/globalobjects.mod \
  $(top_builddir)/MondoMods/globalcharacters.mod \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/memman.mod \
  $(top_builddir)/MondoMods/setxyz.mod \
  $(top_builddir)/MondoMods/processcontrol.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  $(top_builddir)/MondoMods/parsingconstants.mod \
  $(top_builddir)/MondoMods/linalg.mod \
  $(top_builddir)/MondoMods/ainv.mod \
  $(top_builddir)/MondoMods/cholfactor.mod \
  $(top_builddir)/MondoMods/pbc.mod \
  $(top_builddir)/MondoMods/atompairs.mod \
  $(top_builddir)/MondoMods/mondologger.mod

ParseGeomOpt.o : \
  $(top_builddir)/MondoMods/derivedtypes.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  GeomOptKeys.o \
  $(top_builddir)/MondoMods/macros.mod \
  ControlStructures.o \
  OptionKeys.o \
  $(top_builddir)/MondoMods/mondologger.mod

PunchHDF.o : \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/memman.mod \
  $(top_builddir)/MondoMods/pbc.mod \
  $(top_builddir)/MondoMods/indexing.mod \
  $(top_builddir)/MondoMods/atompairs.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  ControlStructures.o \
  OptionKeys.o

PrintParsed.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoLocalVersion.h \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/prettyprint.mod \
  ControlStructures.o \
  $(top_builddir)/MondoMods/utilities.mod

RQI.o : \
  $(top_builddir)/MondoMods/inout.mod \
  $(top_builddir)/MondoMods/macros.mod \
  $(top_builddir)/MondoMods/memman.mod \
  Overlay.o \
  $(top_builddir)/MondoMods/indexing.mod \
  $(top_builddir)/MondoMods/atompairs.mod \
  PunchHDF.o \
  $(top_builddir)/MondoMods/prettyprint.mod \
  $(top_builddir)/MondoMods/globalcharacters.mod \
  $(top_builddir)/MondoMods/globalscalars.mod \
  ControlStructures.o \
  OptionKeys.o \
  Response.o \
  $(top_builddir)/MondoMods/mondologger.mod
