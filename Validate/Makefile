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
	MondoSCF NO+.inp           NO+out            NO+.log
	MondoSCF H2O.inp           H2O.out           H2O.log
	MondoSCF AcO-.inp          AcO-.out          AcO-.log
	MondoSCF PIce.inp          PIce.out          PIce.log
	MondoSCF IceIh.inp         IceIh.out         IceIh.log
	MondoSCF Steam.inp         Steam.out         Steam.log   
	MondoSCF Biotin.inp        Biotin.out        Biotin.log    
	MondoSCF Graphite.inp      Graphite.out      Graphite.log
	MondoSCF Peroxide.inp      Peroxide.out      Peroxide.log
	MondoSCF Histidine.inp     Histidine.out     Histidine.log 
big:	
	MondoSCF Quartz.inp       Quartz.out        Quartz.log
	MondoSCF WaterHexamer.inp WaterHexamer.out  WaterHexamer.log
	MondoSCF BoronNitride.inp BoronNitride.out  BoronNitride.log
	MondoSCF Caffiene.inp     Caffiene.out      Caffine.log
	MondoSCF Diamond.inp      Diamond.out       Diamond.log
	MondoSCF BetaCarotene.inp BetaCarotene.out  BetaCarotene.log
	MondoSCF Chlorophyll.inp  Chlorophyll.out   Chlorophyll.log
	MondoSCF Vancomycin.inp   Vancomycin.out    Vancomycin.log
#
#       Submit serial validation suite to pbs queing system 
#
pbs:
	chmod u+x sub_pbs
	./sub_pbs AcO-
	./sub_pbs BetaCarotene
	./sub_pbs Biotin
	./sub_pbs BoronNitride
	./sub_pbs Caffiene
	./sub_pbs CH2x-6
	./sub_pbs Chlorophyll
	./sub_pbs Diamond
	./sub_pbs Graphite
	./sub_pbs H2O
	./sub_pbs Histidine
	./sub_pbs IceIh
	./sub_pbs NH3
	./sub_pbs NO+
	./sub_pbs Peroxide
	./sub_pbs PIce
	./sub_pbs Quartz
	./sub_pbs Steam
	./sub_pbs Vancomycin
	./sub_pbs WaterHexamer

lfs_pink:
	chmod u+x sub_lsf_pink
	./sub_lfs_pink BetaCarotene   8
	./sub_lfs_pink Biotin         4
	./sub_lfs_pink BoronNitride   8
	./sub_lfs_pink Caffiene       4
	./sub_lfs_pink CH2x-6         8
	./sub_lfs_pink Chlorophyll    8
	./sub_lfs_pink Diamond        4
	./sub_lfs_pink Graphite       4
	./sub_lfs_pink Histidine      4
	./sub_lfs_pink IceIh          4
	./sub_lfs_pink Peroxide       2
	./sub_lfs_pink PIce           4
	./sub_lfs_pink Quartz         4
	./sub_lfs_pink Steam         32
	./sub_lfs_pink Vancomycin    16
	./sub_lfs_pink WaterHexamer   4
#
legacy:	purge 
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
