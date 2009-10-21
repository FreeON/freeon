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

FreeON_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

CPPMISC = -DNEB2

# Link rule for FreeON.
FreeON_LINK = $(FCLINK)

# In case we have ctags installed, make the necessary tags files.
.PHONY : tags
tags :
	ctags --Fortran-kinds=+i+L *.{c,f90} ../Modules/*.{c,f,f90}

# Dependencies.
ParseCommands.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/inout.mod \
  ControlStructures.o \
  $(top_builddir)/Modules/globalcharacters.mod \
  $(top_builddir)/Modules/utilities.mod

ParseOptions.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/memman.mod \
  Numerics.o \
  ConflictDetect.o \
  OptionKeys.o \
  $(top_builddir)/Modules/functionals.mod \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  $(top_builddir)/Modules/globalobjects.mod \
  $(top_builddir)/Modules/processcontrol.mod \
  $(top_builddir)/Modules/globalcharacters.mod \
  $(top_builddir)/Modules/parsingconstants.mod \
  ControlStructures.o \
  $(top_builddir)/Modules/mondologger.mod

ParseGeometries.o : \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/inout.mod \
  OptionKeys.o \
  $(top_builddir)/Modules/prettyprint.mod \
  DynamicsKeys.o \
  GeometryKeys.o \
  ControlStructures.o \
  ls_rmsd.o \
  NEB.o \
  ConflictDetect.o \
  $(top_builddir)/Modules/mondologger.mod

ParseExtraCoords.o : \
  ControlStructures.o \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  OptionKeys.o \
  GeometryKeys.o \
  ParseGeometries.o \
  Massage.o \
  InCoords.o \
  $(top_builddir)/Modules/setxyz.mod \
  $(top_builddir)/Modules/memman.mod

ParseGeomOpt.o : \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  GeomOptKeys.o \
  $(top_builddir)/Modules/macros.mod \
  ControlStructures.o \
  OptionKeys.o \
  $(top_builddir)/Modules/mondologger.mod

ParseProperties.o : \
  $(top_builddir)/Modules/inout.mod \
  ControlStructures.o \
  $(top_builddir)/Modules/mondologger.mod

ParseInput.o : \
  Massage.o \
  ParseBasis.o \
  ParseOptions.o \
  ParseCommands.o \
  ParseDynamics.o \
  ParsePeriodic.o \
  ParseGeometries.o \
  ControlStructures.o \
  $(top_builddir)/Modules/prettyprint.mod \
  ParseParallel.o \
  ParseGeomOpt.o \
  ParseExtraCoords.o \
  ParseProperties.o \
  PrintParsed.o \
  $(top_builddir)/Modules/mondologger.mod

ParsePeriodic.o : \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/pbc.mod \
  $(top_builddir)/Modules/atompairs.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  OptionKeys.o \
  PeriodicKeys.o \
  ControlStructures.o

ParseParallel.o : \
  $(top_builddir)/Modules/inout.mod \
  ParallelKeys.o \
  ControlStructures.o \
  $(top_builddir)/Modules/mondologger.mod

ParseBasis.o : \
  $(top_builddir)/Modules/order.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  ControlStructures.o \
  $(top_builddir)/Modules/basissetparameters.mod \
  $(top_builddir)/Modules/mondologger.mod

ParseDynamics.o : \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  DynamicsKeys.o \
  ControlStructures.o \
  OptionKeys.o \
  $(top_builddir)/Modules/mondologger.mod

ConflictDetect.o : \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/functionals.mod \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/processcontrol.mod \
  ControlStructures.o

Overlay.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/Modules/inout.mod \
  PunchHDF.o \
  ControlStructures.o \
  $(top_builddir)/Modules/utilities.mod

MDynamics.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/Modules/order.mod \
  SCFs.o \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/matfunk.mod \
  Numerics.o \
  $(top_builddir)/Modules/atompairs.mod \
  ControlStructures.o \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/Modules/globalscalars.mod \
  $(top_builddir)/Modules/setxyz.mod \
  DynamicsKeys.o \
  $(top_builddir)/Modules/utilities.mod \
  $(top_builddir)/Modules/mondologger.mod

Response.o : \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/linalg.mod \
  $(top_builddir)/Modules/globalobjects.mod \
  SCFKeys.o \
  Overlay.o \
  PunchHDF.o \
  Numerics.o \
  OptionKeys.o \
  ControlStructures.o \
  NEB.o \
  SCFs.o \
  $(top_builddir)/Modules/setxyz.mod

Hessian.o : \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  $(top_builddir)/Modules/globalobjects.mod \
  $(top_builddir)/Modules/globalcharacters.mod \
  $(top_builddir)/Modules/inout.mod \
  InCoords.o \
  $(top_builddir)/Modules/memman.mod \
  $(top_builddir)/Modules/setxyz.mod \
  $(top_builddir)/Modules/processcontrol.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  $(top_builddir)/Modules/parsingconstants.mod \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/Modules/linalg.mod \
  $(top_builddir)/Modules/ainv.mod \
  $(top_builddir)/Modules/cholfactor.mod \
  ControlStructures.o

MonteCarlo.o : \
  $(top_builddir)/Modules/order.mod \
  SCFs.o \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/matfunk.mod \
  Numerics.o \
  $(top_builddir)/Modules/atompairs.mod \
  ControlStructures.o \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/Modules/globalscalars.mod \
  $(top_builddir)/Modules/setxyz.mod \
  DynamicsKeys.o \
  MDynamics.o

Optimizer.o : \
  SCFs.o \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/matfunk.mod \
  Numerics.o \
  $(top_builddir)/Modules/atompairs.mod \
  ControlStructures.o \
  InCoords.o \
  QUICCA.o \
  GeomOptKeys.o \
  PunchHDF.o \
  $(top_builddir)/Modules/globalscalars.mod \
  Hessian.o \
  $(top_builddir)/Modules/setxyz.mod \
  $(top_builddir)/Modules/utilities.mod \
  $(top_builddir)/Modules/globalscalars.mod

FreeON.o : \
  $(top_builddir)/MondoConfig.h \
  SCFs.o \
  $(top_builddir)/Modules/macros.mod \
  Response.o \
  PunchHDF.o \
  Optimizer.o \
  ParseInput.o \
  ZippyQuote.o \
  MDynamics.o \
  MonteCarlo.o \
  $(top_builddir)/Modules/mondologger.mod \
  RQI.o

QUICCA.o : \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  $(top_builddir)/Modules/globalobjects.mod \
  $(top_builddir)/Modules/globalcharacters.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/memman.mod \
  $(top_builddir)/Modules/linalg.mod \
  $(top_builddir)/Modules/ainv.mod \
  $(top_builddir)/Modules/cholfactor.mod \
  ControlStructures.o \
  PunchHDF.o \
  $(top_builddir)/Modules/setxyz.mod \
  InCoords.o \
  Hessian.o \
  $(top_builddir)/Modules/slatec.mod

ZippyQuote.o : \
  $(top_builddir)/Modules/order.mod \
  ControlStructures.o

NEB.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/Modules/inout.mod \
  ls_rmsd.o \
  $(top_builddir)/Modules/prettyprint.mod \
  ControlStructures.o \
  $(top_builddir)/Modules/order.mod \
  $(top_builddir)/Modules/mondologger.mod

Massage.o : \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/pbc.mod \
  $(top_builddir)/Modules/atompairs.mod \
  OptionKeys.o \
  DynamicsKeys.o \
  GeometryKeys.o \
  ControlStructures.o

SCFs.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/linalg.mod \
  $(top_builddir)/Modules/globalobjects.mod \
  SCFKeys.o \
  Overlay.o \
  DynamicsKeys.o \
  PunchHDF.o \
  Numerics.o \
  OptionKeys.o \
  ControlStructures.o \
  NEB.o \
  $(top_builddir)/Modules/setxyz.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  $(top_builddir)/Modules/mondologger.mod \
  $(top_builddir)/Modules/utilities.mod

InCoords.o : \
  ControlStructures.o \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  $(top_builddir)/Modules/globalobjects.mod \
  $(top_builddir)/Modules/globalcharacters.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/memman.mod \
  $(top_builddir)/Modules/setxyz.mod \
  $(top_builddir)/Modules/processcontrol.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  $(top_builddir)/Modules/parsingconstants.mod \
  $(top_builddir)/Modules/linalg.mod \
  $(top_builddir)/Modules/ainv.mod \
  $(top_builddir)/Modules/cholfactor.mod \
  $(top_builddir)/Modules/pbc.mod \
  $(top_builddir)/Modules/atompairs.mod \
  $(top_builddir)/Modules/mondologger.mod

ParseGeomOpt.o : \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  GeomOptKeys.o \
  $(top_builddir)/Modules/macros.mod \
  ControlStructures.o \
  OptionKeys.o \
  $(top_builddir)/Modules/mondologger.mod

PunchHDF.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/memman.mod \
  $(top_builddir)/Modules/pbc.mod \
  $(top_builddir)/Modules/indexing.mod \
  $(top_builddir)/Modules/atompairs.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  ControlStructures.o \
  OptionKeys.o

PrintParsed.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoLocalVersion.h \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  ControlStructures.o \
  $(top_builddir)/Modules/utilities.mod

RQI.o : \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/macros.mod \
  $(top_builddir)/Modules/memman.mod \
  Overlay.o \
  $(top_builddir)/Modules/indexing.mod \
  $(top_builddir)/Modules/atompairs.mod \
  PunchHDF.o \
  $(top_builddir)/Modules/prettyprint.mod \
  $(top_builddir)/Modules/globalcharacters.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  ControlStructures.o \
  OptionKeys.o \
  Response.o \
  $(top_builddir)/Modules/mondologger.mod

ControlStructures.o : \
  $(top_builddir)/MondoConfig.h

Numerics.o : \
  $(top_builddir)/MondoConfig.h
