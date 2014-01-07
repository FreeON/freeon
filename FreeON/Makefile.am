include $(top_srcdir)/make.inc

bin_PROGRAMS = FreeON

FreeON_SOURCES = \
  ConflictDetect.F90 \
  ControlStructures.F90 \
  DynamicsKeys.F90 \
  FreeON.F90 \
  GeometryKeys.F90 \
  GeomOptKeys.F90 \
  Harmonic.F90 \
  Hessian.F90 \
  InCoords.F90 \
  ls_rmsd.F90 \
  LennardJones.F90 \
  Massage.F90 \
  MDynamics.F90 \
  MLP.F90 \
  MonteCarlo.F90 \
  NEB.F90 \
  Numerics.F90 \
  Optimizer.F90 \
  OptionKeys.F90 \
  Overlay.F90 \
  ParallelKeys.F90 \
  ParseBasis.F90 \
  ParseCommands.F90 \
  ParseDynamics.F90 \
  ParseExtraCoords.F90 \
  ParseGeometries.F90 \
  ParseGeomOpt.F90 \
  ParseInput.F90 \
  ParseOptions.F90 \
  ParseParallel.F90 \
  ParsePeriodic.F90 \
  ParseProperties.F90 \
  PeriodicKeys.F90 \
  PrintParsed.F90 \
  PunchHDF.F90 \
  QUICCA.F90 \
  Response.F90 \
  RQI.F90 \
  SCFKeys.F90 \
  SCFs.F90 \
  Spawn.c \
  ZippyQuote.F90

FreeON_DEPENDENCIES = $(LIBRARY_DEPENDENCY_FILES)

SOURCES_DEPCHECK = $(FreeON_SOURCES)

FreeON_LDADD = $(LDADD_MODULES)

AM_CPPFLAGS += \
  -DFREEON_BASISSETS=\"$(datadir)/$(PACKAGE)/BasisSets\" \
  -DFREEON_HOME=\"${prefix}\"

CPPMISC = -DNEB2

# Link rule for FreeON.
#FreeON_LINK = $(FCLINK)

TAGS_FILES = \
	$(top_srcdir)/Modules/*.c \
	$(top_srcdir)/Modules/*.F \
	$(top_srcdir)/Modules/*.F90

AM_CTAGSFLAGS = --Fortran-kinds=+i+L

# Dependencies.
ParseCommands.o : \
  $(top_builddir)/MondoConfig.h \
  ControlStructures.o

ParseOptions.o : \
  $(top_builddir)/MondoConfig.h \
  Numerics.o \
  ConflictDetect.o \
  OptionKeys.o \
  ControlStructures.o

ParseGeometries.o : \
  OptionKeys.o \
  DynamicsKeys.o \
  GeometryKeys.o \
  ControlStructures.o \
  ls_rmsd.o \
  NEB.o \
  ConflictDetect.o

ParseExtraCoords.o : \
  ControlStructures.o \
  OptionKeys.o \
  GeometryKeys.o \
  ParseGeometries.o \
  Massage.o \
  InCoords.o

ParseGeomOpt.o : \
  GeomOptKeys.o \
  ControlStructures.o \
  OptionKeys.o

ParseProperties.o : \
  ControlStructures.o

ParseInput.o : \
  Massage.o \
  ParseBasis.o \
  ParseOptions.o \
  ParseCommands.o \
  ParseDynamics.o \
  ParsePeriodic.o \
  ParseGeometries.o \
  ControlStructures.o \
  ParseParallel.o \
  ParseGeomOpt.o \
  ParseExtraCoords.o \
  ParseProperties.o \
  PrintParsed.o

ParsePeriodic.o : \
  OptionKeys.o \
  PeriodicKeys.o \
  ControlStructures.o

ParseParallel.o : \
  ParallelKeys.o \
  ControlStructures.o \
  $(top_builddir)/MondoConfig.h

ParseBasis.o : \
  ControlStructures.o

ParseDynamics.o : \
  DynamicsKeys.o \
  ControlStructures.o \
  OptionKeys.o

ConflictDetect.o : \
  ControlStructures.o

Overlay.o : \
  $(top_builddir)/MondoConfig.h \
  PunchHDF.o \
  ControlStructures.o

MDynamics.o : \
  $(top_builddir)/MondoConfig.h \
  SCFs.o \
  Numerics.o \
  ControlStructures.o \
  GeomOptKeys.o \
  PunchHDF.o \
  DynamicsKeys.o

Response.o : \
  SCFKeys.o \
  Overlay.o \
  PunchHDF.o \
  Numerics.o \
  OptionKeys.o \
  ControlStructures.o \
  NEB.o \
  SCFs.o

Hessian.o : \
  InCoords.o \
  GeomOptKeys.o \
  PunchHDF.o \
  ControlStructures.o

MonteCarlo.o : \
  SCFs.o \
  Numerics.o \
  ControlStructures.o \
  GeomOptKeys.o \
  PunchHDF.o \
  DynamicsKeys.o \
  MDynamics.o

Optimizer.o : \
  SCFs.o \
  Numerics.o \
  ControlStructures.o \
  InCoords.o \
  QUICCA.o \
  GeomOptKeys.o \
  PunchHDF.o \
  Hessian.o

FreeON.o : \
  $(top_builddir)/MondoConfig.h \
  SCFs.o \
  Response.o \
  PunchHDF.o \
  Optimizer.o \
  ParseInput.o \
  ZippyQuote.o \
  MDynamics.o \
  MonteCarlo.o \
  RQI.o

QUICCA.o : \
  ControlStructures.o \
  PunchHDF.o \
  InCoords.o \
  Hessian.o

ZippyQuote.o : \
  ControlStructures.o

NEB.o : \
  $(top_builddir)/MondoConfig.h \
  ls_rmsd.o \
  ControlStructures.o

Massage.o : \
  OptionKeys.o \
  DynamicsKeys.o \
  GeometryKeys.o \
  ControlStructures.o

SCFs.o : \
  $(top_builddir)/MondoConfig.h \
  SCFKeys.o \
  Overlay.o \
  DynamicsKeys.o \
  PunchHDF.o \
  Numerics.o \
  OptionKeys.o \
  ControlStructures.o \
  NEB.o \
  LennardJones.o \
  Harmonic.o

InCoords.o : \
  ControlStructures.o

ParseGeomOpt.o : \
  GeomOptKeys.o \
  ControlStructures.o \
  OptionKeys.o

PunchHDF.o : \
  $(top_builddir)/MondoConfig.h \
  ControlStructures.o \
  OptionKeys.o

PrintParsed.o : \
  $(top_builddir)/MondoConfig.h \
  $(top_builddir)/MondoLocalVersion.h \
  ControlStructures.o

RQI.o : \
  Overlay.o \
  PunchHDF.o \
  ControlStructures.o \
  OptionKeys.o \
  Response.o

ControlStructures.o : \
  $(top_builddir)/MondoConfig.h

Numerics.o : \
  $(top_builddir)/MondoConfig.h

LennardJones.o :\
  ControlStructures.o

Harmonic.o : \
  ControlStructures.o
