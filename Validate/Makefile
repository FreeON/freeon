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
	MondoSCF Ice.inp           Ice.out           Ice.log   
	MondoSCF Steam.inp         Steam.out         Steam.log   
	MondoSCF Water.inp         Water.out         Water.log
	MondoSCF Urea.inp          Urea.out          Urea.log
	MondoSCF Biotin.inp        Biotin.out        Biotin.log    
	MondoSCF Graphite.inp      Graphite.out      Graphite.log
	MondoSCF Peroxide.inp      Peroxide.out      Peroxide.log
	MondoSCF Histidine.inp     Histidine.out     Histidine.log 
big:	
	MondoSCF Caffiene.inp     Caffiene.out      Caffine.log
	MondoSCF Diamond.inp      Diamond.out       Diamond.log
	MondoSCF BetaCarotene.inp BetaCarotene.out  BetaCarotene.log
	MondoSCF Chlorophyll.inp  Chlorophyll.out;  Chlorophyll.log
	MondoSCF Vancomycin.inp   Vancomycin.out;   Vancomycin.log
#
#       Submit serial validation suite to pbs queing system 
#
pbs:
	chmod u+x sub_pbs
	./sub_pbs Ice            1
	./sub_pbs Steam          1
	./sub_pbs Water          1
	./sub_pbs Urea           1
	./sub_pbs Biotin         1
	./sub_pbs Graphite       1
	./sub_pbs Peroxide       1
	./sub_pbs Histidine      1
	./sub_pbs Caffiene       1
	./sub_pbs Vancomycin     1
	./sub_pbs BetaCarotene   1
	./sub_pbs Diamond        1
	./sub_pbs Chlorophyll    1
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
	rm -f *.xsf *.xyz *.pdb* *.MDO
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f *.o* *.pbs* 
	rm -f *.inf.ascii
	rm -rf  core core.* core*
	rm -f PI*
	rm -f *.usertime.* *.pcsamp.* *.*_Info64 *_Info64 *.aggregate.* 
	ln -s /dev/null core
