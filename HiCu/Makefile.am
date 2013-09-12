include $(top_srcdir)/make.inc

bin_PROGRAMS = HiCu XCForce

HiCu_SOURCES    = Grid.f90 Thresholds.f90 RhoTree.f90 CubeTree.f90 TreeWalk.f90 KxcGen.f90 ParaHiCu.f90 HiCu.f90
XCForce_SOURCES = Grid.f90 Thresholds.f90 RhoTree.f90 CubeTree.f90 TreeWalk.f90 KxcGen.f90 ParaHiCu.f90 dXCBlok.f90 XCForce.f90

HiCu_DEPENDENCIES    = $(LIBRARY_DEPENDENCY_FILES)
XCForce_DEPENDENCIES = $(LIBRARY_DEPENDENCY_FILES)

SOURCES_DEPCHECK = \
  $(HiCu_SOURCES) \
  $(XCForce_SOURCES)

EXTRA_DIST = \
  ExplicitBraElements.Inc \
  ExplicitLeafContribution.Inc \
  ExplicitLeafPopulation.Inc \
  GeneralLeafContribution.Inc \
  GeneralLeafPopulation.Inc \
  MMA/CubeRules/Rule5.Inc \
  MMA/CubeRules/Rule7.Inc \
  MMA/CubeRules/RuleB.Inc \
  MMA/CubeRules/RuleC.Inc \
  MMA/CubeRules/RuleD.Inc \
  MMA/CubeRules/RuleE.Inc

LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

EXTRA_INCLUDES = -I$(top_srcdir)/Modules/MMA/Functionals
AM_CPPFLAGS += -DRULE7 -DEXPLICIT_SOURCE

# Dependencies.
CubeTree.o : \
  Grid.o \
	RhoTree.o \
	Thresholds.o

TreeWalk.o : \
  CubeTree.o

KxcGen.o : \
  Thresholds.o \
	CubeTree.o \
	TreeWalk.o

XCForce.o : \
  CubeTree.o \
  dXCBlok.o

dXCBlok.o : \
  Thresholds.o \
	CubeTree.o \
	TreeWalk.o

RhoTree.o : \
  Thresholds.o \
  Grid.o

HiCu.o : \
  Thresholds.o \
	RhoTree.o \
	CubeTree.o \
  KxcGen.o \
  $(top_builddir)/MondoConfig.h

ParaHiCu.o : \
  $(top_builddir)/MondoConfig.h
