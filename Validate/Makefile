##############################################################################
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
##############################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
Libs =  $(F90_LIBRARIES)
#
TstObjs=Test.o
#
ReleaseFiles=Test.x
#
all:	test small 
#
test:	$(TstObjs)
	$(LD) $(LDFLAGS) Test $(TstObjs) $(LibOpts) $(Libs) 
	rm -rf Validation.log
#------------------------------------------------------------------
#
small:
	$(MONDO_HOME)/Exec/MondoSCF SiF4.inp    SiF4.out;  \
                             ./Test SiF4.out    SiF4.ref
	$(MONDO_HOME)/Exec/MondoSCF Urea.inp    Urea.out;  \
                             ./Test Urea.out    Urea.ref
	$(MONDO_HOME)/Exec/MondoSCF Water.inp   Water.out;  \
                             ./Test Water.out   Water.ref 
	$(MONDO_HOME)/Exec/MondoSCF Biotin.inp  Biotin.out; \
                             ./Test Biotin.out  Biotin.ref
#
big:	small	
#----------------------------------------------
#       Submit validation suite to pbs queing 
#
small_pbs:
	chmod u+x sub_pbs
	./sub_pbs SiF4
	./sub_pbs Urea 
	./sub_pbs Water
	./sub_pbs Biotin
#
big_pbs:small_pbs
	./sub_pbs Caffiene
	./sub_pbs BetaCarotene
	./sub_pbs Chlorophyll
	./sub_pbs Vancomycin
#
#--------------------------------------------
#
release:clean $(ReleaseFiles) 
	rm -f *.out
	rm -f *.log
#
purge:	clean 
	rm -f $(REMOVESRCE)
	rm -f *.out
	rm -f *.log
	rm -f *.geo
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f *.o* *.pbs*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core
