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
all:	small 
#
full:	all big
#
#       Run validation suite sequentialy
#
serial:	small big
#
small:
	m Water.inp   Water.out      Water.log
	m NH3.inp     NH3.out        NH3.log
	m Biotin.inp  Biotin.out     Biotin.log    
	m Diamond.inp Diamond.out    Diamond.log
	m Na2Cl2.inp  Na2Cl2.out     Na2Cl2.log
#
big:	
	m Caffiene.inp     Caffiene.out      Caffine.log
	m MgO.inp          MgO.out           MgO.log
	m BetaCarotene.inp BetaCarotene.out  BetaCarotene.log
	m Chlorophyll.inp  Chlorophyll.out;  Chlorophyll.log
	m Vancomycin.inp   Vancomycin.out;   Vancomycin.log
#
#       Submit validation suite to pbs queing system 
#
pbs:
	chmod u+x sub_pbs
	./sub_pbs Water
	./sub_pbs NH3
	./sub_pbs Biotin
	./sub_pbs Diamond
	./sub_pbs Na2Cl2
	./sub_pbs MgO
	./sub_pbs Caffiene
	./sub_pbs BetaCarotene
	./sub_pbs Chlorophyll
	./sub_pbs Vancomycin
#
release:  purge 
	rm -rf BALA-PERIODIC
	rm -rf MgO-PERIODIC
#
purge:	clean 
	rm -f $(REMOVESRCE)
	rm -f Test
	rm -f \#*
	rm -f *.out
	rm -f *.log
	rm -f *.car
	rm -f *.xsf *.xyz *.pdb*
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f *.o* *.pbs* 
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core
