####################################################################################
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
####################################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
CPPMISC = 
#
EXTRA_INCLUDES= -I$(MONDO_HOME)/QCTC/ -I$(MONDO_HOME)/HiCu/
#
VisDXObjs=PotMapRho.o 
MisObjs=$(MONDO_HOME)/HiCu/CubeTree.o    \
	$(MONDO_HOME)/HiCu/RhoTree.o     \
	$(MONDO_HOME)/HiCu/Grid.o        \
	$(MONDO_HOME)/HiCu/Thresholds.o  \
	$(MONDO_HOME)/QCTC/TreeWalk.o    \
	$(MONDO_HOME)/QCTC/MondoPoles.o  \
	$(MONDO_HOME)/QCTC/PoleTree.o    \
	$(MONDO_HOME)/QCTC/PBCFarField.o \
	$(MONDO_HOME)/QCTC/PFFTensor.o   \
	$(MONDO_HOME)/QCTC/XLate.o       \
	$(MONDO_HOME)/QCTC/CTraX.o       \
	$(MONDO_HOME)/QCTC/Thresholds.o  \
	$(MONDO_HOME)/QCTC/Globals.o   
#
ReleaseFiles=PotMapRho.x 
#
all:	vis
#
release:clean $(ReleaseFiles)
#
vis:	$(VisDXObjs) 
	$(LD) $(LDFLAGS) PotMapRho $(VisDXObjs) $(MisObjs) $(F90_LIBRARIES) 
	mv PotMapRho $(MONDO_EXEC)/PotMapRho
#
purge:	clean
	rm -f $(REMOVESRCE)
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/PotMapRho
	ls 
