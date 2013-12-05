include $(top_srcdir)/make.inc

bin_PROGRAMS = IRInv ODA CPSCFStatus DDIIS FBuild TC2Response AInv LowdinO RHeqs SP2 TS4 SCFstats P2Use DIIS

IRInv_SOURCES       = IRInv.F90 DenMatMethods.F90
ODA_SOURCES         = ODA.F90 DenMatMethods.F90
CPSCFStatus_SOURCES = CPSCFStatus.F90
DDIIS_SOURCES       = DDIIS.F90
FBuild_SOURCES      = FBuild.F90
TC2Response_SOURCES = TC2Response.F90 DenMatMethods.F90 DenMatResponse.F90
AInv_SOURCES        = BlokAInv.F90 DenMatMethods.F90
LowdinO_SOURCES     = LowdinO.F90 DenMatMethods.F90
RHeqs_SOURCES       = RHeqs.F90 DenMatMethods.F90
SP2_SOURCES         = FillPBlok.F90 SP2.F90 DenMatMethods.F90
TS4_SOURCES         = FillPBlok.F90 TS4.F90 DenMatMethods.F90
SCFstats_SOURCES    = SCFStatus.F90
P2Use_SOURCES       = FillPBlok.F90 P2Use.F90 DenMatMethods.F90

if REALLY_OLD_DIIS
DIIS_SOURCES = DIIS_2004-01-16.F90
else
if OLD_DIIS
DIIS_SOURCES = DIIS_1.3.F90
else
DIIS_SOURCES = DIIS.F90
endif
endif

IRInv_DEPENDENCIES       = $(LIBRARY_DEPENDENCY_FILES)
ODA_DEPENDENCIES         = $(LIBRARY_DEPENDENCY_FILES)
CPSCFStatus_DEPENDENCIES = $(LIBRARY_DEPENDENCY_FILES)
DDIIS_DEPENDENCIES       = $(LIBRARY_DEPENDENCY_FILES)
FBuild_DEPENDENCIES      = $(LIBRARY_DEPENDENCY_FILES)
TC2Response_DEPENDENCIES = $(LIBRARY_DEPENDENCY_FILES)
AInv_DEPENDENCIES        = $(LIBRARY_DEPENDENCY_FILES)
LowdinO_DEPENDENCIES     = $(LIBRARY_DEPENDENCY_FILES)
RHeqs_DEPENDENCIES       = $(LIBRARY_DEPENDENCY_FILES)
SP2_DEPENDENCIES         = $(LIBRARY_DEPENDENCY_FILES)
TS4_DEPENDENCIES         = $(LIBRARY_DEPENDENCY_FILES)
SCFstats_DEPENDENCIES    = $(LIBRARY_DEPENDENCY_FILES)
P2Use_DEPENDENCIES       = $(LIBRARY_DEPENDENCY_FILES)
DIIS_DEPENDENCIES        = $(LIBRARY_DEPENDENCY_FILES)

SOURCES_DEPCHECK = \
  $(IRInv_SOURCES) \
  $(ODA_SOURCES) \
  $(CPSCFStatus_SOURCES) \
  $(DDIIS_SOURCES) \
  $(FBuild_SOURCES) \
  $(TC2Response_SOURCES) \
  $(AInv_SOURCES) \
  $(LowdinO_SOURCES) \
  $(RHeqs_SOURCES) \
  $(SP2_SOURCES) \
  $(TS4_SOURCES) \
  $(SCFstats_SOURCES) \
  $(P2Use_SOURCES) \
  $(DIIS_SOURCES)

LDADD = $(LDADD_MODULES)

TAGS_FILES = \
	$(top_srcdir)/Modules/*.c
	$(top_srcdir)/Modules/*.F
	$(top_srcdir)/Modules/*.F90

AM_CTAGSFLAGS = --Fortran-kinds=+i+L

# Dependencies.
IRInv.o          : DenMatMethods.o
ODA.o            : $(top_builddir)/MondoConfig.h DenMatMethods.o
TC2Response.o    : DenMatMethods.o DenMatResponse.o
SP2.o            : $(top_builddir)/MondoConfig.h DenMatMethods.o
TS4.o            : DenMatMethods.o
DenMatResponse.o : DenMatMethods.o
LowdinO.o        : DenMatMethods.o
P2Use.o          : $(top_builddir)/MondoConfig.h DenMatMethods.o FillPBlok.o
RHeqs.o          : DenMatMethods.o
BlokAInv.o       : $(top_builddir)/MondoConfig.h DenMatMethods.o
SCFstats.o       : $(top_builddir)/MondoConfig.h
DenMatMethods.o  : $(top_builddir)/MondoConfig.h
DIIS.o           : $(top_builddir)/MondoConfig.h
FBuild.o         : $(top_builddir)/MondoConfig.h
