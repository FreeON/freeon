#
#     Matt Challacombe
#     Los Alamos National Laboratory
#     Copyright 2000, The University of California
#
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

CPPMISC=
#-DNODIIS

Gues_Objs  = FillPBlok.o P2Use.o 
SDMM_Objs  = FillPBlok.o SDMM.o               
Stat_Objs  = SCFStatus.o          
RHEq_Objs  = RHeqs.o              
FBuild_Objs= FBuild.o       
DIIS_Objs  = DIIS.o               
Dev_Objs   = ToDev.o               

All:	Stat Gues SDMM FBuild  DIIS RHEq

Stat:	$(Stat_Objs)
	$(LD) $(LDFLAGS) SCFstats $(Stat_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv SCFstats $(MONDO_HOME)/Exec/SCFstats

Gues:	$(Gues_Objs)
	$(LD) $(LDFLAGS) P2Use $(Gues_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv P2Use $(MONDO_HOME)/Exec/P2Use

SDMM:	$(SDMM_Objs)
	$(LD) $(LDFLAGS) SDMM $(SDMM_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv SDMM $(MONDO_HOME)/Exec/SDMM

RHEq:	$(RHEq_Objs)
	$(LD) $(LDFLAGS) RHeqs $(RHEq_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv RHeqs $(MONDO_HOME)/Exec/RHeqs 

FBuild:	$(FBuild_Objs)
	$(LD) $(LDFLAGS) FBuild $(FBuild_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv FBuild $(MONDO_HOME)/Exec/FBuild

DIIS:	$(DIIS_Objs)
	$(LD) $(LDFLAGS) DIIS $(DIIS_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv DIIS $(MONDO_HOME)/Exec/DIIS

Dev:	$(Dev_Objs)
	$(LD) $(LDFLAGS) ToDev $(Dev_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv ToDev $(MONDO_HOME)/Exec/ToDev

clean:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/SDMM
	rm -f $(MONDO_HOME)/Exec/P2Use
	rm -f $(MONDO_HOME)/Exec/FBuild
	rm -f $(MONDO_HOME)/Exec/RHeqs
	rm -f $(MONDO_HOME)/Exec/SCFstats
	ln -s /dev/null core
	ls -l

DIIS.o:		            $(MONDO_HOME)/Libs/libMondoMods.a
FBuild.o:	            $(MONDO_HOME)/Libs/libMondoMods.a
P2Use.o:	FillPBlok.o $(MONDO_HOME)/Libs/libMondoMods.a
FillPBlok.o:	            $(MONDO_HOME)/Libs/libMondoMods.a
RHeqs.o:		    $(MONDO_HOME)/Libs/libMondoMods.a
SCFStatus.o:                $(MONDO_HOME)/Libs/libMondoMods.a
SDMM.o:		            $(MONDO_HOME)/Libs/libMondoMods.a
