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

FreeON_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES)

SOURCES_DEPCHECK = $(FreeON_SOURCES)

FreeON_LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

CPPMISC = -DNEB2

# Link rule for FreeON.
FreeON_LINK = $(FCLINK)

# Dependencies.
ConflictDetect.o: ControlStructures.o
ControlStructures.o: 
FreeON.o: MDynamics.o MonteCarlo.o Optimizer.o ParseInput.o PunchHDF.o RQI.o \
	Response.o SCFs.o ZippyQuote.o
GeomOptKeys.o: 
Hessian.o: ControlStructures.o GeomOptKeys.o InCoords.o PunchHDF.o
InCoords.o: ControlStructures.o
Massage.o: ControlStructures.o DynamicsKeys.o GeometryKeys.o OptionKeys.o
MDynamics.o: ControlStructures.o DynamicsKeys.o GeomOptKeys.o Numerics.o \
	PunchHDF.o SCFs.o
MLP.o: 
MonteCarlo.o: ControlStructures.o DynamicsKeys.o GeomOptKeys.o MDynamics.o \
	Numerics.o PunchHDF.o SCFs.o
NEB.o: ControlStructures.o ls_rmsd.o
Numerics.o: 
Optimizer.o: ControlStructures.o GeomOptKeys.o Hessian.o InCoords.o \
	Numerics.o PunchHDF.o QUICCA.o SCFs.o
Overlay.o: ControlStructures.o PunchHDF.o
ParseBasis.o: ControlStructures.o
ParseCommands.o: ControlStructures.o
ParseDynamics.o: ControlStructures.o DynamicsKeys.o OptionKeys.o
ParseExtraCoords.o: ControlStructures.o GeometryKeys.o InCoords.o Massage.o \
	OptionKeys.o ParseGeometries.o
ParseGeometries.o: ConflictDetect.o ControlStructures.o DynamicsKeys.o \
	GeometryKeys.o NEB.o OptionKeys.o ls_rmsd.o
ParseGeomOpt.o: ControlStructures.o GeomOptKeys.o OptionKeys.o
ParseInput-preprocessed.o: ControlStructures.o Massage.o ParseBasis.o \
	ParseCommands.o ParseDynamics.o ParseExtraCoords.o ParseGeomOpt.o \
	ParseGeometries.o ParseOptions.o ParseParallel.o ParsePeriodic.o \
	ParseProperties.o PrintParsed.o
ParseInput.o: ControlStructures.o Massage.o ParseBasis.o ParseCommands.o \
	ParseDynamics.o ParseExtraCoords.o ParseGeomOpt.o ParseGeometries.o \
	ParseOptions.o ParseParallel.o ParsePeriodic.o ParseProperties.o \
	PrintParsed.o
ParseOptions.o: ConflictDetect.o ControlStructures.o Numerics.o OptionKeys.o
ParseParallel.o: ControlStructures.o ParallelKeys.o
ParsePeriodic.o: ControlStructures.o OptionKeys.o PeriodicKeys.o
ParseProperties.o: ControlStructures.o
PrintParsed.o: ControlStructures.o
PunchHDF.o: ControlStructures.o OptionKeys.o
QUICCA.o: ControlStructures.o Hessian.o InCoords.o PunchHDF.o
Response.o: ControlStructures.o NEB.o Numerics.o OptionKeys.o Overlay.o \
	PunchHDF.o SCFKeys.o SCFs.o
RQI.o: ControlStructures.o OptionKeys.o Overlay.o PunchHDF.o
SCFs.o: ControlStructures.o DynamicsKeys.o NEB.o Numerics.o OptionKeys.o \
	Overlay.o PunchHDF.o SCFKeys.o
ZippyQuote.o: ControlStructures.o
