#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC = -DRULE7 -DEXPLICIT_SOURCE
#
EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/LookUpTables \
                -I$(MONDO_HOME)/MondoMods/MMA/Functionals
#
CommonObjs = Grid.o Thresholds.o RhoTree.o CubeTree.o TreeWalk.o 
HObjs      = $(CommonObjs) KxcGen.o HiCu.o 
XFObjs     = $(CommonObjs) dXCBlok.o XCForce.o
#
VObjs  = VRhoMD.o 
#
all:	VR HC XF
#
clean:	CHC CXF
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	ls -lS
#
purge:	clean CVR
#-----------------------------------------------------------------------
HC:	$(HObjs)
	$(LD) $(LDFLAGS) HiCu $(HObjs) $(LibOpts) $(F90_LIBRARIES)
	mv HiCu $(MONDO_HOME)/Exec
#
CHC:	
	rm -f $(MONDO_HOME)/Exec/HiCu
#-----------------------------------------------------------------------
XF:     $(XFObjs)
	$(LD) $(LDFLAGS) XCForce $(XFObjs) $(LibOpts) $(F90_LIBRARIES)
	mv XCForce $(MONDO_HOME)/Exec
#
CXF:
	rm -f $(MONDO_HOME)/Exec/XCForce
#-----------------------------------------------------------------------
VR:	$(VObjs)
	$(LD) $(LDFLAGS) VRho $(VObjs) $(LibOpts) $(F90_LIBRARIES)
	./VRho
	rm -rf VRho VRho.o VRho.f90
#
CVR:	
	rm -rf Explicit*.Inc
#
CubeTree.o:	Grid.o                      \
		GeneralLeafPopulation.Inc   \
		GeneralLeafContribution.Inc   $(MONDO_HOME)/Libs/libMondoMods.a
Grid.o:		MMA/CubeRules/Rule5.Inc \
                MMA/CubeRules/Rule7.Inc        $(MONDO_HOME)/Libs/libMondoMods.a
HiCu.o:		KxcGen.o TreeWalk.o \
                CubeTree.o RhoTree.o  \
		$(MONDO_HOME)/Libs/libMondoMods.a
XCForce.o:	dXCBlok.o TreeWalk.o \
		CubeTree.o RhoTree.o \
		$(MONDO_HOME)/Libs/libMondoMods.a
RhoTree.o:	$(MONDO_HOME)/Libs/libMondoMods.a
VRhoMD.o:	$(MONDO_HOME)/Libs/libMondoMods.a
