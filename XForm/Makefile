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
CPPOPTS=#-DEXTREME_DEBUG -DEXTREME_PRINT_DEBUG
#
Libs =  $(F90_LIBRARIES) $(LAPACK)
#
BlkAInvObjs = BlokAInv.o  
LowdinOObjs = LowdinO.o 
#
ReleaseFiles=BlokAInv.x Lowdin.x
#
all:	Blk Lwd
#
release:clean $(RELEASEFILES)
#
purge:	clean
	rm -f $(REMOVESRCE)
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/AInv
	rm -f $(MONDO_HOME)/Exec/LOrthog
	ln -s /dev/null core
	ls -l 
#
Blk:	$(BlkAInvObjs)
	$(LD) $(LDFLAGS) AInv $(BlkAInvObjs) $(LibOpts) $(Libs)
	mv AInv $(MONDO_HOME)/Exec/AInv
#
Lwd:	$(LowdinOObjs)
	$(LD) $(LDFLAGS) LowdinO $(LowdinOObjs) $(LibOpts) $(Libs)
	mv LowdinO $(MONDO_HOME)/Exec/LowdinO
#
LowdinO.o:	$(MONDO_HOME)/Libs/libMondoMods.a
BlokAInv.o:	$(MONDO_HOME)/Libs/libMondoMods.a
