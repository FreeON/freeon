#####################################################################################
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
#####################################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC = -DRULE7 -DEXPLICIT_SOURCE 
#         -DOLDSPLIT
#
EXTRA_INCLUDES=-I$(MONDO_HOME)/MondoMods/MMA/Functionals
#
CommonObjs = Grid.o Thresholds.o RhoTree.o CubeTree.o TreeWalk.o 
HObjs      = $(CommonObjs) KxcGen.o ParaHiCu.o HiCu.o 
XFObjs     = $(CommonObjs) KxcGen.o ParaHiCu.o dXCBlok.o XCForce.o
#
VObjs  = VRhoMD.o 
#
ReleaseFiles=Grid.x Thresholds.x RhoTree.x    \
             ParaHiCu.x CubeTree.x TreeWalk.x KxcGen.x   \
             HiCu.x dXCBlok.x XCForce.x VRhoMD.x 
#
all:	HC XF
#
release:clean $(ReleaseFiles)
#
purge:	clean
	rm -f $(REMOVESRCE)
#
clean:	CHC CXF
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ls -lS
#-----------------------------------------------------------------------
HC:	$(HObjs)
	$(LD) $(LDFLAGS) HiCu $(HObjs) $(F90_LIBRARIES) $(LAPACK)
	mv HiCu $(MONDO_EXEC)
#
CHC:	
	rm -f $(MONDO_EXEC)/HiCu
#-----------------------------------------------------------------------
XF:     $(XFObjs)
	$(LD) $(LDFLAGS) XCForce $(XFObjs) $(F90_LIBRARIES) $(LAPACK)
	mv XCForce $(MONDO_EXEC)
#
CXF:
	rm -f $(MONDO_EXEC)/XCForce
#-----------------------------------------------------------------------
VR:	$(VObjs)
	$(LD) $(LDFLAGS) VRho $(VObjs) $(F90_LIBRARIES)
	./VRho
	rm -rf VRho VRho.o VRho.f90
#
CVR:	
	rm -rf Explicit*.Inc
#
ParaHiCu.o:	CubeTree.o TreeWalk.o KxcGen.o $(MONDO_LIB)/libMondoMods.a
CubeTree.o:	Grid.o GeneralLeafPopulation.Inc   \
		GeneralLeafContribution.Inc   $(MONDO_LIB)/libMondoMods.a
Grid.o:		MMA/CubeRules/Rule5.Inc  MMA/CubeRules/Rule7.Inc  $(MONDO_LIB)/libMondoMods.a
HiCu.o:		KxcGen.o TreeWalk.o \
                ParaHiCu.o CubeTree.o RhoTree.o  \
		$(MONDO_LIB)/libMondoMods.a
XCForce.o:	dXCBlok.o TreeWalk.o CubeTree.o KxcGen.o ParaHiCu.o RhoTree.o \
		$(MONDO_LIB)/libMondoMods.a
RhoTree.o:	$(MONDO_LIB)/libMondoMods.a
VRhoMD.o:	$(MONDO_LIB)/libMondoMods.a
