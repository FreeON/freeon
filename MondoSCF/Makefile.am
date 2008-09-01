include $(top_srcdir)/make.inc

bin_PROGRAMS = MondoSCF

MondoSCF_SOURCES = \
  ControlStructures.f90 \
  ls_rmsd.f90 \
  NEB.f90 \
  Numerics.f90 \
  ParseCommands.f90 \
  OptionKeys.f90 \
  ConflictDetect.f90 \
  ParseOptions.f90 \
  DynamicsKeys.f90 \
  ParseDynamics.f90 \
  GeometryKeys.f90 \
  ParseProperties.f90 \
  InCoords.f90 \
  ParseGeometries.f90 \
  PeriodicKeys.f90 \
  ParsePeriodic.f90 \
  ParseBasis.f90 \
  ParallelKeys.f90 \
  ParseParallel.f90 \
  Massage.f90 \
  GeomOptKeys.f90 \
  ParseGeomOpt.f90 \
  ParseExtraCoords.f90 \
  ParseInput.f90 \
  PunchHDF.f90 \
  Overlay.f90 \
  RQI.f90 \
  PrintParsed.f90 \
  SCFKeys.f90 \
  SCFs.f90 \
  MLP.f90 \
  Response.f90 \
  Hessian.f90 \
  QUICCA.f90 \
  Optimizer.f90 \
  MDynamics.f90 \
  MonteCarlo.f90 \
  ZippyQuote.f90 \
  GetPID.c \
  Spawn.c \
  MondoSCF.f90

MondoSCF_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES)

SOURCES_DEPCHECK = $(MondoSCF_SOURCES)

EXTRA_INCLUDES = -I../MondoMods

MondoSCF_LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

CPPMISC = -DNEB2

# Link rule for MondoSCF.
MondoSCF_LINK = $(FCLINK)

# Dependencies.
ParseCommands.o : \
  ControlStructures.o \
  $(top_srcdir)/MondoConfig.h

ParseOptions.o : \
	Numerics.o \
  ConflictDetect.o \
	OptionKeys.o \
	ControlStructures.o

ParseGeometries.o : \
  ConflictDetect.o \
  OptionKeys.o \
  DynamicsKeys.o \
  GeometryKeys.o \
  ControlStructures.o \
  ls_rmsd.o \
  NEB.o

ParseExtraCoords.o : \
	ControlStructures.o \
	OptionKeys.o \
	GeometryKeys.o \
	ParseGeometries.o \
	InCoords.o \
	Massage.o

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
	ParseProperties.o

ParsePeriodic.o : \
	OptionKeys.o \
  PeriodicKeys.o \
	ControlStructures.o

ParseParallel.o : \
	ParallelKeys.o \
  ControlStructures.o

ParseBasis.o : \
  ControlStructures.o

ParseDynamics.o : \
  DynamicsKeys.o \
	ControlStructures.o \
	OptionKeys.o

ConflictDetect.o : \
  ControlStructures.o

Overlay.o : \
  PunchHDF.o

RQI.o : \
  PunchHDF.o \
  Overlay.o 

MDynamics.o : \
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
	ControlStructures.o \
	GeomOptKeys.o \
	PunchHDF.o \
	DynamicsKeys.o \
	MDynamics.o

Optimizer.o : \
  SCFs.o \
	QUICCA.o

MondoSCF.o : \
  SCFs.o \
  Response.o \
  PunchHDF.o \
  Optimizer.o \
  ParseInput.o \
  ZippyQuote.o \
  PrintParsed.o \
  MDynamics.o \
  MonteCarlo.o \
  RQI.o

QUICCA.o : \
  PunchHDF.o \
	Hessian.o

ZippyQuote.o : \
  ControlStructures.o

NEB.o : \
	ls_rmsd.o \
  ControlStructures.o

Massage.o : \
	OptionKeys.o \
        DynamicsKeys.o \
	GeometryKeys.o \
	ControlStructures.o

SCFs.o : \
	SCFKeys.o \
	Overlay.o \
	DynamicsKeys.o \
	PunchHDF.o \
	Numerics.o \
	OptionKeys.o \
	ControlStructures.o \
	NEB.o

InCoords.o : \
  ControlStructures.o

ParseGeomOpt.o : \
  GeomOptKeys.o

PunchHDF.o : \
  ControlStructures.o \
	OptionKeys.o

PrintParsed.o : \
  ControlStructures.o \
  $(top_srcdir)/MondoConfig.h
