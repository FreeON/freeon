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
EXTRA_INCLUDES= 
#
CTraXObjs=CTraX107.o CTraX306.o CTraX406.o CTraX505.o CTraX603.o \
	CTraX607.o CTraX704.o CTraX206.o CTraX307.o CTraX407.o   \
	CTraX506.o CTraX604.o CTraX701.o CTraX705.o CTraX207.o   \
	CTraX404.o CTraX503.o CTraX507.o CTraX605.o CTraX702.o   \
        CTraX706.o CTraX305.o CTraX405.o CTraX504.o CTraX602.o   \
        CTraX606.o CTraX703.o CTraX707.o
#
HGTraXObjs=AuxInt10.o AuxInt8.o HGTraX304.o HGTraX405.o  \
	HGTraX505.o HGTraX604.o HGTraX702.o AuxInt11.o   \
	AuxInt9.o HGTraX305.o HGTraX406.o HGTraX506.o    \
	HGTraX605.o HGTraX703.o AuxInt12.o HGTraX106.o   \
	HGTraX306.o HGTraX502.o HGTraX601.o HGTraX606.o  \
	HGTraX704.o AuxInt13.o HGTraX205.o HGTraX403.o   \
	HGTraX503.o HGTraX602.o HGTraX700.o HGTraX705.o  \
	AuxInt7.o  HGTraX206.o HGTraX404.o HGTraX504.o   \
	HGTraX603.o HGTraX701.o HGTraX706.o
#
ModObjs=PoleNodeType.o ERIGlobals.o PoleGlobals.o 
#
ReleaseFiles=CTraX107.x CTraX306.x CTraX406.x CTraX505.x CTraX603.x \
	CTraX607.x CTraX704.x CTraX206.x CTraX307.x CTraX407.x   \
	CTraX506.x CTraX604.x CTraX701.x CTraX705.x CTraX207.x   \
	CTraX404.x CTraX503.x CTraX507.x CTraX605.x CTraX702.x   \
        CTraX706.x CTraX305.x CTraX405.x CTraX504.x CTraX602.x   \
        CTraX606.x CTraX703.x CTraX707.x                         \
	AuxInt10.x AuxInt8.x HGTraX304.x HGTraX405.x             \
	HGTraX505.x HGTraX604.x HGTraX702.x AuxInt11.x           \
	AuxInt9.x HGTraX305.x HGTraX406.x HGTraX506.x            \
	HGTraX605.x HGTraX703.x AuxInt12.x HGTraX106.x           \
	HGTraX306.x HGTraX502.x HGTraX601.x HGTraX606.x          \
	HGTraX704.x AuxInt13.x HGTraX205.x HGTraX403.x           \
	HGTraX503.x HGTraX602.x HGTraX700.x HGTraX705.x          \
	AuxInt7.x  HGTraX206.x HGTraX404.x HGTraX504.x           \
	HGTraX603.x HGTraX701.x HGTraX706.x                      \
        PoleNodeType.x ERIGlobals.x PoleGlobals.x 
#-----------------------------------------------------------------------------------------------
#
all:	$(ModObjs) $(CTraXObjs) $(HGTraXObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libQCTC.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libQCTC.a	
#
purge:	
	rm -f $(MONDO_LIB)/libQCTC.a
	rm -f $(REMOVESRCE)
#
release:purge $(ReleaseFiles)
#
clean:	
	rm -rf $(REMOVEMISC) 
	rm -rf $(REMOVEMODS)
	rm -rf *.o
#---------------------------------------------------------------------------------------------
#       DEPENDENCIES
#
#------------------------------------------------------------------------------
