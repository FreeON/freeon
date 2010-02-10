include $(top_srcdir)/make.inc

bin_PROGRAMS = Make1E MakeM MakeU MakeS MakeT MakeRho SForce TForce UForce MForce THess

Make1E_SOURCES  = Int1E.f90 Make1X.f90 Make1E.f90
MakeM_SOURCES   = MBlok.f90 MakeM.f90
MakeU_SOURCES   = UBlok.f90 MakeU.f90
MakeS_SOURCES   = SBlok.f90 MakeS.f90
MakeT_SOURCES   = TBlok.f90 MakeT.f90
MakeRho_SOURCES = RhoTools.f90 RhoBlok.f90 MakeRho.f90
SForce_SOURCES  = BlokTrWdS.f90 SForce.f90
TForce_SOURCES  = BlokTrPdT.f90 TForce.f90
MForce_SOURCES  = BlokTrPdM.f90 MForce.f90
THess_SOURCES   = BlokTrPd2T.f90 THess.f90
UForce_SOURCES  = UBlok.f90 UForce.f90

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

LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

EXTRA_INCLUDES = -I$(srcdir)/MMA

# In case we have ctags installed, make the necessary tags files.
.PHONY : tags
tags :
	ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.f90 ../Modules/*.{c,f,f90}

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
