############################################################################
#--  This code is part of the MondoSCF suite of programs for linear scaling 
#    electronic structure theory and ab initio molecular dynamics.
#
#--  Copyright (c) 2001, the Regents of the University of California.  
#    This SOFTWARE has been authored by an employee or employees of the 
#    University of California, operator of the Los Alamos National Laboratory 
#    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
#    The U.S. Government has rights to use, reproduce, and distribute this 
#    SOFTWARE.  The public may copy, distribute, prepare derivative works 
#    and publicly display this SOFTWARE without charge, provided that this 
#    Notice and any statement of authorship are reproduced on all copies.  
#    Neither the Government nor the University makes any warranty, express 
#    or implied, or assumes any liability or responsibility for the use of 
#    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
#    such modified SOFTWARE should be clearly marked, so as not to confuse 
#    it with the version available from LANL.  The return of derivative works
#    to the primary author for integration and general release is encouraged. 
#    The first publication realized with the use of MondoSCF shall be
#    considered a joint work.  Publication of the results will appear
#    under the joint authorship of the researchers nominated by their
#    respective institutions. In future publications of work performed
#    with MondoSCF, the use of the software shall be properly acknowledged,
#    e.g. in the form "These calculations have been performed using MondoSCF, 
#    a suite of programs for linear scaling electronic structure theory and
#    ab initio molecular dynamics", and given appropriate citation.  
############################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC=#-DNODIIS
#
Gues_Objs  = FillPBlok.o P2Use.o 
SDMM_Objs  = FillPBlok.o SDMM.o               
Stat_Objs  = SCFStatus.o          
RHEq_Objs  = RHeqs.o              
FBuild_Objs= FBuild.o       
DIIS_Objs  = DIIS.o               
Dev_Objs   = ToDev.o               
#
ReleaseFiles=FillPBlok.x P2Use.x FillPBlok.x \
             SDMM.x SCFStatus.x RHeqs.x      \
             FBuild.x DIIS.x ToDev.x
#
all:	Stat Gues SDMM FBuild  DIIS RHEq
#
release:clean $(ReleaseFiles)
#
purge:	clean 
	rm -r $(REMOVESRCE)
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/SDMM
	rm -f $(MONDO_HOME)/Exec/P2Use
	rm -f $(MONDO_HOME)/Exec/FBuild
	rm -f $(MONDO_HOME)/Exec/RHeqs
	rm -f $(MONDO_HOME)/Exec/SCFstats
	rm -f $(MONDO_HOME)/Exec/DIIS
	rm -f $(MONDO_HOME)/Exec/ToDev
	ln -s /dev/null core
	ls -l
#
Stat:	$(Stat_Objs)
	$(LD) $(LDFLAGS) SCFstats $(Stat_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv SCFstats $(MONDO_HOME)/Exec/SCFstats
#
Gues:	$(Gues_Objs)
	$(LD) $(LDFLAGS) P2Use $(Gues_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv P2Use $(MONDO_HOME)/Exec/P2Use
#
SDMM:	$(SDMM_Objs)
	$(LD) $(LDFLAGS) SDMM $(SDMM_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv SDMM $(MONDO_HOME)/Exec/SDMM
#
RHEq:	$(RHEq_Objs)
	$(LD) $(LDFLAGS) RHeqs $(RHEq_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv RHeqs $(MONDO_HOME)/Exec/RHeqs 
#
FBuild:	$(FBuild_Objs)
	$(LD) $(LDFLAGS) FBuild $(FBuild_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv FBuild $(MONDO_HOME)/Exec/FBuild
#
DIIS:	$(DIIS_Objs)
	$(LD) $(LDFLAGS) DIIS $(DIIS_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv DIIS $(MONDO_HOME)/Exec/DIIS
#
Dev:	$(Dev_Objs)
	$(LD) $(LDFLAGS) ToDev $(Dev_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv ToDev $(MONDO_HOME)/Exec/ToDev
#
DIIS.o:		            $(MONDO_HOME)/Libs/libMondoMods.a
FBuild.o:	            $(MONDO_HOME)/Libs/libMondoMods.a
P2Use.o:	FillPBlok.o $(MONDO_HOME)/Libs/libMondoMods.a
FillPBlok.o:	            $(MONDO_HOME)/Libs/libMondoMods.a
RHeqs.o:		    $(MONDO_HOME)/Libs/libMondoMods.a
SCFStatus.o:                $(MONDO_HOME)/Libs/libMondoMods.a
SDMM.o:		            $(MONDO_HOME)/Libs/libMondoMods.a
