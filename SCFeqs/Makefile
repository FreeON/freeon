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
CPPMISC=#-DWORK_ON_DEGENERACY #-DNODIIS #-DSPATIAL_THRESHOLDING #-DFIND_CONDA 
#
BlkAInvObjs = BlokAInv.o  
LowdinOObjs = LowdinO.o 
#
SDMM_Objs  = FillPBlok.o SDMM.o DenMatMethods.o
SP2_Objs   = FillPBlok.o SP2.o  DenMatMethods.o
SP4_Objs   = FillPBlok.o SP4.o  DenMatMethods.o
NT4_Objs   = FillPBlok.o NT4.o  DenMatMethods.o
PM_Objs    = FillPBlok.o PM.o   DenMatMethods.o
#
Gues_Objs  = FillPBlok.o P2Use.o 
Stat_Objs  = SCFStatus.o          
RHEq_Objs  = RHeqs.o  DenMatMethods.o
FBuild_Objs= FBuild.o                   
DIIS_Objs  = DIIS.o               
Dev_Objs   = ToDev.o  
#
ReleaseFiles=BlokAInv.x LowdinO.x RHeqs.x    \
             SDMM.x SP2.x SP4.x NT4.x PM.x   \
             FillPBlok.x P2Use.x FillPBlok.x \
             SCFStatus.x FBuild.x DIIS.x ToDev.x 
#
all:	Blk Lwd RHEq SDMM SP2 SP4 NT4 PM Stat Gues FBuild DIIS 
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
	rm -f $(MONDO_EXEC)/SDMM
	rm -f $(MONDO_EXEC)/P2Use
	rm -f $(MONDO_EXEC)/FBuild
	rm -f $(MONDO_EXEC)/RHeqs
	rm -f $(MONDO_EXEC)/SCFstats
	rm -f $(MONDO_EXEC)/DIIS
	rm -f $(MONDO_EXEC)/ToDev
	rm -f $(MONDO_EXEC)/SP2
	rm -f $(MONDO_EXEC)/SP4
	rm -f $(MONDO_EXEC)/NT4
	rm -f $(MONDO_EXEC)/PM
	rm -f $(MONDO_EXEC)/AInv
	rm -f $(MONDO_EXEC)/LowdinO
	ln -s /dev/null core
#
#	Orthogonalization transforms
#
Blk:	$(BlkAInvObjs)
	$(LD) $(LDFLAGS) AInv $(BlkAInvObjs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv AInv $(MONDO_EXEC)/AInv
#
Lwd:	$(LowdinOObjs)
	$(LD) $(LDFLAGS) LowdinO $(LowdinOObjs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv LowdinO $(MONDO_EXEC)/LowdinO
#
#       SCF Solvers
#
RHEq:	$(RHEq_Objs)
	$(LD) $(LDFLAGS) RHeqs $(RHEq_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv RHeqs $(MONDO_EXEC)/RHeqs 
#
SDMM:	$(SDMM_Objs)
	$(LD) $(LDFLAGS) SDMM $(SDMM_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv SDMM $(MONDO_EXEC)/SDMM
#
SP2:	$(SP2_Objs)
	$(LD) $(LDFLAGS) SP2  $(SP2_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv SP2 $(MONDO_EXEC)/SP2
#
SP4:	$(SP4_Objs)
	$(LD) $(LDFLAGS) SP4  $(SP4_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv SP4 $(MONDO_EXEC)/SP4
#
NT4:	$(NT4_Objs)
	$(LD) $(LDFLAGS) NT4  $(NT4_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv NT4 $(MONDO_EXEC)/NT4
#
PM:	$(PM_Objs)
	$(LD) $(LDFLAGS) PM  $(PM_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv PM $(MONDO_EXEC)/PM
#
#	Misc
#
Stat:	$(Stat_Objs)
	$(LD) $(LDFLAGS) SCFstats $(Stat_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv SCFstats $(MONDO_EXEC)/SCFstats
#
Gues:	$(Gues_Objs)
	$(LD) $(LDFLAGS) P2Use $(Gues_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv P2Use $(MONDO_EXEC)/P2Use
#
FBuild:	$(FBuild_Objs)
	$(LD) $(LDFLAGS) FBuild $(FBuild_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv FBuild $(MONDO_EXEC)/FBuild
#
DIIS:	$(DIIS_Objs)
	$(LD) $(LDFLAGS) DIIS $(DIIS_Objs) $(LibOpts) $(F90_LIBRARIES) $(LAPACK)
	mv DIIS $(MONDO_EXEC)/DIIS
#
Dev:	$(Dev_Objs)
	$(LD) $(LDFLAGS) ToDev $(Dev_Objs) $(LibOpts) $(F90_LIBRARIES) 
	mv ToDev $(MONDO_EXEC)/ToDev
#
DIIS.o:		            $(MONDO_HOME)/Libs/libMondoMods.a
FBuild.o:	            $(MONDO_HOME)/Libs/libMondoMods.a
P2Use.o:	FillPBlok.o $(MONDO_HOME)/Libs/libMondoMods.a
FillPBlok.o:	            $(MONDO_HOME)/Libs/libMondoMods.a
RHeqs.o:		    $(MONDO_HOME)/Libs/libMondoMods.a
SCFStatus.o:                $(MONDO_HOME)/Libs/libMondoMods.a
SDMM.o:		            $(MONDO_HOME)/Libs/libMondoMods.a
SP2.o:                      $(MONDO_HOME)/Libs/libMondoMods.a
SP4.o:                      $(MONDO_HOME)/Libs/libMondoMods.a
NT4.o:                      $(MONDO_HOME)/Libs/libMondoMods.a
PM.o:                       $(MONDO_HOME)/Libs/libMondoMods.a
LowdinO.o:		    $(MONDO_HOME)/Libs/libMondoMods.a
BlokAInv.o:		    $(MONDO_HOME)/Libs/libMondoMods.a
