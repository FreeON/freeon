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
full:	all  big
#
pbs:	test small_pbs periodic_pbs big_pbs
#------------------------------------------------------------------
#       Code to parse *.ref and *.out, checking numbers for correctness
test:	$(TstObjs)
	$(LD) $(LDFLAGS) Test $(TstObjs) $(LibOpts) $(Libs) 
	rm -rf Validation.log
#------------------------------------------------------------------
#       Run validation suite sequentialy
serial:	small periodic big
#
small:
	m SiF4.inp    SiF4.out;              ./Test SiF4.out    SiF4.ref
	m Urea.inp    Urea.out;              ./Test Urea.out    Urea.ref
	m Water.inp   Water.out;             ./Test Water.out   Water.ref 
	m Biotin.inp  Biotin.out;            ./Test Biotin.out  Biotin.ref
#
periodic:
	m Mg0.inp          MgO.out;          ./Test MgO.out           MgO.ref
#
big:	
	m Caffiene.inp     Caffiene.out;     ./Test Caffiene.out      Caffiene.ref
	m BetaCarotene.inp BetaCarotene.out; ./Test BetaCarotene.out  BetaCarotene.ref 
	m Chlorophyll.inp  Chlorophyll.out;  ./Test Chlorophyll.out   Chlorophyll.ref
	m Vancomycin.inp   Vancomycin.out;   ./Test Vancomycin.out    Vancomycin.ref
#-----------------------------------------------------------------------------------
#       Submit validation suite to pbs queing system 
#
small_pbs:
	chmod u+x sub_pbs
	./sub_pbs SiF4
	./sub_pbs Urea 
	./sub_pbs Water
	./sub_pbs Biotin
#
periodic_pbs:
	chmod u+x sub_pbs
	./sub_pbs MgO
#
big_pbs:
	./sub_pbs Caffiene
	./sub_pbs BetaCarotene
	./sub_pbs Chlorophyll
	./sub_pbs Vancomycin
#--------------------------------------------
release:purge $(ReleaseFiles) 
#
purge:	clean 
	rm -f $(REMOVESRCE)
	rm -f Test
	rm -f \#*
	rm -f *.out
	rm -f *.log
	rm -f *.car
	rm -f *.xyz
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f *.o* *.pbs* *.pdb*
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core
