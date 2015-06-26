include $(top_srcdir)/make.inc

bin_PROGRAMS = Make1E MakeM MakeU MakeS MakeT MakeRho SForce TForce UForce MForce THess

Make1E_SOURCES  = Int1E.F90 Make1X.F90 Make1E.F90
MakeM_SOURCES   = MBlok.F90 MakeM.F90
MakeU_SOURCES   = UBlok.F90 MakeU.F90
MakeS_SOURCES   = SBlok.F90 MakeS.F90
MakeT_SOURCES   = TBlok.F90 MakeT.F90
MakeRho_SOURCES = RhoTools.F90 RhoBlok.F90 MakeRho.F90
SForce_SOURCES  = BlokTrWdS.F90 SForce.F90
TForce_SOURCES  = BlokTrPdT.F90 TForce.F90
MForce_SOURCES  = BlokTrPdM.F90 MForce.F90
THess_SOURCES   = BlokTrPd2T.F90 THess.F90
UForce_SOURCES  = UBlok.F90 UForce.F90

Make1E_DEPENDENCIES  = $(LIBRARY_DEPENDENCY_FILES)
MakeM_DEPENDENCIES   = $(LIBRARY_DEPENDENCY_FILES)
MakeU_DEPENDENCIES   = $(LIBRARY_DEPENDENCY_FILES)
MakeS_DEPENDENCIES   = $(LIBRARY_DEPENDENCY_FILES)
MakeT_DEPENDENCIES   = $(LIBRARY_DEPENDENCY_FILES)
MakeRho_DEPENDENCIES = $(LIBRARY_DEPENDENCY_FILES)
SForce_DEPENDENCIES  = $(LIBRARY_DEPENDENCY_FILES)
TForce_DEPENDENCIES  = $(LIBRARY_DEPENDENCY_FILES)
MForce_DEPENDENCIES  = $(LIBRARY_DEPENDENCY_FILES)
THess_DEPENDENCIES   = $(LIBRARY_DEPENDENCY_FILES)
UForce_DEPENDENCIES  = $(LIBRARY_DEPENDENCY_FILES)

SOURCES_DEPCHECK = \
  $(SMake1E_SOURCES) \
  $(MakeM_SOURCES) \
  $(MakeU_SOURCES) \
  $(MakeS_SOURCES) \
  $(MakeT_SOURCES) \
  $(MakeRho_SOURCES) \
  $(SForce_SOURCES) \
  $(TForce_SOURCES) \
  $(MForce_SOURCES) \
  $(THess_SOURCES) \
  $(UForce_SOURCES)

EXTRA_DIST = \
  MMA/Omega1.Inc \
  MMA/Omega2.Inc \
  MMA/QQuad128.Inc \
  MMA/QQuad32.Inc \
  MMA/QQuad64.Inc

LDADD = $(LDADD_MODULES) $(LAPACK_LIBS) $(BLAS_LIBS)

EXTRA_INCLUDES = -I$(srcdir)/MMA

TAGS_FILES = \
	$(top_srcdir)/Modules/*.c
	$(top_srcdir)/Modules/*.F
	$(top_srcdir)/Modules/*.F90

AM_CTAGSFLAGS = --Fortran-kinds=+i+L

# Dependencies.
MakeU.o   : UBlok.o
MakeS.o   : $(top_builddir)/MondoConfig.h SBlok.o
MakeT.o   : $(top_builddir)/MondoConfig.h TBlok.o
MakeM.o   : MBlok.o
Make1E.o  : Make1X.o
Make1X.o  : Int1E.o
MakeRho.o : $(top_builddir)/MondoConfig.h RhoBlok.o RhoTools.o
SForce.o  : $(top_builddir)/MondoConfig.h BlokTrWdS.o
TForce.o  : $(top_builddir)/MondoConfig.h BlokTrPdT.o
MForce.o  : BlokTrPdM.o
THess.o   : BlokTrPd2T.o
UForce.o  : UBlok.o
RhoBlok.o : RhoTools.o
