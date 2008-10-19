include $(top_srcdir)/make.inc

bin_PROGRAMS = FreeON

FreeON_SOURCES = \
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
  FreeON.f90

FreeON_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES)

SOURCES_DEPCHECK = $(MondoSCF_SOURCES)

EXTRA_INCLUDES = -I../MondoMods

FreeON_LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

CPPMISC = -DNEB2

# Link rule for FreeON.
FreeON_LINK = $(FCLINK)

# Dependencies.

ConflictDetect.o: ControlStructures.o
ControlStructures.o: 
GeomOptKeys.o: 
Hessian.o: ControlStructures.o GeomOptKeys.o InCoords.o PunchHDF.o
InCoords.o: ControlStructures.o
Massage.o: ControlStructures.o DynamicsKeys.o GeometryKeys.o OptionKeys.o
MDynamics.o: ControlStructures.o DynamicsKeys.o GeomOptKeys.o Numerics.o \
	PunchHDF.o SCFs.o
MLP.o: 
MondoSCF.o: MDynamics.o MonteCarlo.o Optimizer.o ParseInput.o PrintParsed.o \
	PunchHDF.o RQI.o Response.o SCFs.o ZippyQuote.o
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
ParseInput.o: ControlStructures.o Massage.o ParseBasis.o ParseCommands.o \
	ParseDynamics.o ParseExtraCoords.o ParseGeomOpt.o ParseGeometries.o \
	ParseOptions.o ParseParallel.o ParsePeriodic.o ParseProperties.o
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
