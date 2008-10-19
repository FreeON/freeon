include $(top_srcdir)/make.inc

bin_PROGRAMS = IRInv ODA CPSCFStatus DDIIS FBuild TC2Response AInv LowdinO RHeqs SP2 TS4 SCFstats P2Use DIIS

IRInv_SOURCES       = IRInv.f90 DenMatMethods.f90
ODA_SOURCES         = ODA.f90 DenMatMethods.f90
CPSCFStatus_SOURCES = CPSCFStatus.f90
DDIIS_SOURCES       = DDIIS.f90
FBuild_SOURCES      = FBuild.f90
TC2Response_SOURCES = TC2Response.f90 DenMatMethods.f90 DenMatResponse.f90
AInv_SOURCES        = BlokAInv.f90 DenMatMethods.f90
LowdinO_SOURCES     = LowdinO.f90 DenMatMethods.f90
RHeqs_SOURCES       = RHeqs.f90 DenMatMethods.f90
SP2_SOURCES         = FillPBlok.f90 SP2.f90 DenMatMethods.f90
TS4_SOURCES         = FillPBlok.f90 TS4.f90 DenMatMethods.f90
SCFstats_SOURCES    = SCFStatus.f90
P2Use_SOURCES       = FillPBlok.f90 P2Use.f90 DenMatMethods.f90

if REALLY_OLD_DIIS
DIIS_SOURCES = DIIS_2004-01-16.f90
else
if OLD_DIIS
DIIS_SOURCES = DIIS_1.3.f90
else
DIIS_SOURCES = DIIS.f90
endif
endif

IRInv_DEPENDENCIES       = $(MondoMods_DEPENDENCY_FILES)
ODA_DEPENDENCIES         = $(MondoMods_DEPENDENCY_FILES)
CPSCFStatus_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES)
DDIIS_DEPENDENCIES       = $(MondoMods_DEPENDENCY_FILES)
FBuild_DEPENDENCIES      = $(MondoMods_DEPENDENCY_FILES)
TC2Response_DEPENDENCIES = $(MondoMods_DEPENDENCY_FILES)
AInv_DEPENDENCIES        = $(MondoMods_DEPENDENCY_FILES)
LowdinO_DEPENDENCIES     = $(MondoMods_DEPENDENCY_FILES)
RHeqs_DEPENDENCIES       = $(MondoMods_DEPENDENCY_FILES)
SP2_DEPENDENCIES         = $(MondoMods_DEPENDENCY_FILES)
TS4_DEPENDENCIES         = $(MondoMods_DEPENDENCY_FILES)
SCFstats_DEPENDENCIES    = $(MondoMods_DEPENDENCY_FILES)
P2Use_DEPENDENCIES       = $(MondoMods_DEPENDENCY_FILES)
DIIS_DEPENDENCIES        = $(MondoMods_DEPENDENCY_FILES)

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

LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

EXTRA_INCLUDES = -I../MondoMods -IMMA

# Dependencies.
IRInv.o          : DenMatMethods.o
ODA.o            : DenMatMethods.o
TC2Response.o    : DenMatMethods.o DenMatResponse.o
SP2.o            : DenMatMethods.o
TS4.o            : DenMatMethods.o
DenMatResponse.o : DenMatMethods.o
LowdinO.o        : DenMatMethods.o
P2Use.o          : $(top_srcdir)/MondoConfig.h DenMatMethods.o FillPBlok.o
RHeqs.o          : DenMatMethods.o
BlokAInv.o       : DenMatMethods.o
