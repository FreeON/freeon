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
ODA_Objs   = ODA.o DenMatMethods.o
MIX_Objs   = MIX.o DenMatMethods.o
SDMM_Objs  = FillPBlok.o SDMM.o DenMatMethods.o
SP2_Objs   = FillPBlok.o SP2.o  DenMatMethods.o
SP4_Objs   = FillPBlok.o SP4.o  DenMatMethods.o
TS4_Objs   = FillPBlok.o TS4.o  DenMatMethods.o
PM_Objs    = FillPBlok.o PM.o   DenMatMethods.o
#
Gues_Objs  = FillPBlok.o P2Use.o DenMatMethods.o
Stat_Objs  = SCFStatus.o 
RStat_Objs = CPSCFStatus.o
RHEq_Objs  = RHeqs.o  DenMatMethods.o
FBuild_Objs= FBuild.o                   
DIIS_Objs  = DIIS.o               
DDIIS_Objs = DDIIS.o                #vw
TC2Response_Objs = TC2Response.o DenMatMethods.o   #vw
Dev_Objs   = ToDev.o  
#
ReleaseFiles=BlokAInv.x LowdinO.x RHeqs.x    \
             SDMM.x SP2.x SP4.x TS4.x PM.x   \
             FillPBlok.x P2Use.x FillPBlok.x \
             SCFStatus.x FBuild.x DIIS.x ToDev.x \
	     DDIIS.x TC2Response.x CPSCFStatus.x \
             DenMatMethods.x 

#
all:	ODA RStat DDIIS FBuild TC2Response Blk Lwd RHEq SP2 TS4 Stat Gues DIIS  
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
	rm -f $(MONDO_EXEC)/CPSCFStatus
	rm -f $(MONDO_EXEC)/DIIS
	rm -f $(MONDO_EXEC)/DDIIS     #vw
	rm -f $(MONDO_EXEC)/ToDev
	rm -f $(MONDO_EXEC)/SP2
	rm -f $(MONDO_EXEC)/SP4
	rm -f $(MONDO_EXEC)/TS4
	rm -f $(MONDO_EXEC)/PM
	rm -f $(MONDO_EXEC)/AInv
	rm -f $(MONDO_EXEC)/LowdinO
	ln -s /dev/null core
#
#       EXTRAPOLATIONATORS
# 
DIIS:	$(DIIS_Objs)
	$(LD) $(LDFLAGS) DIIS $(DIIS_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv DIIS $(MONDO_EXEC)/DIIS
#
DDIIS:	$(DDIIS_Objs)
	$(LD) $(LDFLAGS) DDIIS $(DDIIS_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv DDIIS $(MONDO_EXEC)/DDIIS
#
ODA:	$(ODA_Objs)
	$(LD) $(LDFLAGS) ODA $(ODA_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv ODA $(MONDO_EXEC)/ODA

MIX:	$(MIX_Objs)
	$(LD) $(LDFLAGS) MIX $(MIX_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv MIX $(MONDO_EXEC)/MIX
#
#	Orthogonalization transforms
#
Blk:	$(BlkAInvObjs)
	$(LD) $(LDFLAGS) AInv $(BlkAInvObjs)  $(F90_LIBRARIES) $(LAPACK)
	mv AInv $(MONDO_EXEC)/AInv
#
Lwd:	$(LowdinOObjs)
	$(LD) $(LDFLAGS) LowdinO $(LowdinOObjs)  $(F90_LIBRARIES) $(LAPACK)
	mv LowdinO $(MONDO_EXEC)/LowdinO
#
#       SCF Solvers
#
RHEq:	$(RHEq_Objs)
	$(LD) $(LDFLAGS) RHeqs $(RHEq_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv RHeqs $(MONDO_EXEC)/RHeqs 
#
SDMM:	$(SDMM_Objs)
	$(LD) $(LDFLAGS) SDMM $(SDMM_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv SDMM $(MONDO_EXEC)/SDMM
#
SP2:	$(SP2_Objs)
	$(LD) $(LDFLAGS) SP2  $(SP2_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv SP2 $(MONDO_EXEC)/SP2
#
SP4:	$(SP4_Objs)
	$(LD) $(LDFLAGS) SP4  $(SP4_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv SP4 $(MONDO_EXEC)/SP4
#
TS4:	$(TS4_Objs)
	$(LD) $(LDFLAGS) TS4  $(TS4_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv TS4 $(MONDO_EXEC)/TS4
#
PM:	$(PM_Objs)
	$(LD) $(LDFLAGS) PM  $(PM_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv PM $(MONDO_EXEC)/PM
#
#	Misc
#
Stat:	$(Stat_Objs)
	$(LD) $(LDFLAGS) SCFstats $(Stat_Objs)  $(F90_LIBRARIES) 
	mv SCFstats $(MONDO_EXEC)/SCFstats
#
RStat:	$(RStat_Objs)
	$(LD) $(LDFLAGS) CPSCFStatus $(RStat_Objs)  $(F90_LIBRARIES) 
	mv CPSCFStatus $(MONDO_EXEC)/CPSCFStatus
#
Gues:	$(Gues_Objs)
	$(LD) $(LDFLAGS) P2Use $(Gues_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv P2Use $(MONDO_EXEC)/P2Use
#
FBuild:	$(FBuild_Objs)
	$(LD) $(LDFLAGS) FBuild $(FBuild_Objs)  $(F90_LIBRARIES) 
	mv FBuild $(MONDO_EXEC)/FBuild
#
TC2Response:	$(TC2Response_Objs)
	$(LD) $(LDFLAGS) TC2Response  $(TC2Response_Objs)  $(F90_LIBRARIES) $(LAPACK)
	mv TC2Response $(MONDO_EXEC)/TC2Response
#
Dev:	$(Dev_Objs)
	$(LD) $(LDFLAGS) ToDev $(Dev_Objs)  $(F90_LIBRARIES) 
	mv ToDev $(MONDO_EXEC)/ToDev
#
DIIS.o:		            $(MONDO_LIB)/libMondoMods.a
DDIIS.o:		    $(MONDO_LIB)/libMondoMods.a
TC2Response.o: DenMatMethods.o	$(MONDO_LIB)/libMondoMods.a
FBuild.o:	            $(MONDO_LIB)/libMondoMods.a
P2Use.o:	DenMatMethods.o FillPBlok.o $(MONDO_LIB)/libMondoMods.a 
FillPBlok.o:	            $(MONDO_LIB)/libMondoMods.a
RHeqs.o:		    $(MONDO_LIB)/libMondoMods.a
SCFStatus.o:                $(MONDO_LIB)/libMondoMods.a
ODA.o:	        DenMatMethods.o $(MONDO_LIB)/libMondoMods.a
MIX.o:	        DenMatMethods.o $(MONDO_LIB)/libMondoMods.a
SP2.o:          DenMatMethods.o $(MONDO_LIB)/libMondoMods.a
SP4.o:          DenMatMethods.o $(MONDO_LIB)/libMondoMods.a
TS4.o:          DenMatMethods.o $(MONDO_LIB)/libMondoMods.a
PM.o:           DenMatMethods.o $(MONDO_LIB)/libMondoMods.a
LowdinO.o:		    $(MONDO_LIB)/libMondoMods.a
BlokAInv.o:		    $(MONDO_LIB)/libMondoMods.a
