#############################################################################
#--  This source code is part of the MondoSCF suite of programs for 
#    linear scaling electronic structure theory and ab initio molecular 
#    dynamics.
#
#--  Matt Challacombe
#    Los Alamos National Laboratory
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
#############################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
Libs =  
#
HessObjs = BFGSMisc.o BFGSHess.o  
NewObjs  = NewStep.o
#
ReleaseFiles=BFGSMisc.x BFGSHess.x NewStep.x
#
all:	Hess New
#
release:clean $(ReleaseFiles)
	rm -rf OLD
#
purge:	clean 
	rm -f $(REMOVESRCE)
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/BFGSHess
	rm -f $(MONDO_EXEC)/NewStep
	ln -s /dev/null core
	ls -l 
#
Hess:	$(HessObjs)
	$(LD) $(LDFLAGS) BFGSHess $(HessObjs) $(F90_LIBRARIES) $(LAPACK)
	mv BFGSHess $(MONDO_EXEC)/BFGSHess
#
New:	$(NewObjs)
	$(LD) $(LDFLAGS) NewStep $(NewObjs) $(F90_LIBRARIES) $(LAPACK)
	mv NewStep $(MONDO_EXEC)/NewStep
#
BFGSHess.o:	$(MONDO_LIB)/libMondoMods.a

