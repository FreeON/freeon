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
#    MAIN MAKEFILE FOR MondoSCF
#    Author: Matt Challacombe
#------------------------------------------------------------------------------
include $(MONDO_HOME)/Includes/RemoveAll
#
all:	 Env mm s x 1 2 e o n
#
clean:	CExec cmm cs cx c1 c2 ce cn 
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
#
purge:	clean PScr PWrk pmm p2 PLib
#
Env:	
	cat $(MONDO_HOME)/Includes/CopyrightNotice.txt
	sleep 2
#
2:	q h # o
#
c2:	cq ch # co
#
p2:	ph pq
#
#       LIBRARIES
#
mm:	
	$(MAKE) -C MondoMods all
#
cmm:	
	$(MAKE) -i -C MondoMods clean
#
pmm:	
	$(MAKE) -i -C MondoMods purge
#
PLib:	
	rm -rf $(MONDO_HOME)/Libs/*
#
#       FRONT END
#
s:	
	$(MAKE) -C MondoSCF
#
cs:	
	$(MAKE) -i -C MondoSCF clean
#
#       GENERATION OF ORTHOGONALIZATION MATRICES
#
x:	
	$(MAKE) -C XForm
#
cx:	
	$(MAKE) -i -C XForm clean
#
#       SOLVING THE SELF-CONSISTENT-FIELD EQUATIONS
#
e:	
	$(MAKE) -C SCFeqs
#
ce:	
	$(MAKE) -i -C SCFeqs clean
#
#       GEOMETRY OPTIMIZATION USING NEWTONS METHOD
#
n:	
	$(MAKE) -C QuNew
#
cn:	
	$(MAKE) -i -C QuNew clean
#
#       ONE ELECTRON ROUTINES
#
1:
	$(MAKE) -C OneE
#
c1:	
	$(MAKE) -i -C OneE clean
#
#       ORDER N EXCHANGE
#
o:	
	$(MAKE) -C ONX
#
co:	
	$(MAKE) -i -C ONX clean 
#
#       DIRECT J BUILD
#
j:	
	$(MAKE) -C DirectJ
#
cj:	
	$(MAKE) -C DirectJ clean
#
pj:	
	$(MAKE) -C DirectJ purge
#
#       QUANTUM CHEMICAL TREE CODE
#
q:	
	$(MAKE) -C QCTC
#
cq:	
	$(MAKE) -C QCTC clean
#
pq:	
	$(MAKE) -C QCTC purge
#
#       HIERARCHICAL CUBATURE
#
h:	
	$(MAKE) -C HiCu
#
ch:	
	$(MAKE) -i -C HiCu clean
#
ph:	
	$(MAKE) -i -C HiCu purge
#
#       CLEAN EXECUTABLES
#
CExec:	
	rm -rf $(MONDO_HOME)/Exec/*
#
#       PURGE CURRENT WORK DIRECTORY
#
PWrk:	
	$(MAKE) -i -C $(MONDO_WORK) clean
#
#       PURGE SCRATCH DIRECTORY
#
PScr:	
	rm  -rf $(MONDO_SCRATCH)/*.S
	rm  -rf $(MONDO_SCRATCH)/*.T
	rm  -rf $(MONDO_SCRATCH)/*.Rho
	rm  -rf $(MONDO_SCRATCH)/*.V
	rm  -rf $(MONDO_SCRATCH)/*.J
	rm  -rf $(MONDO_SCRATCH)/*.E
	rm  -rf $(MONDO_SCRATCH)/*.F
	rm  -rf $(MONDO_SCRATCH)/*.OrthoF
	rm  -rf $(MONDO_SCRATCH)/*.F_DIIS
	rm  -rf $(MONDO_SCRATCH)/*.D	
	rm  -rf $(MONDO_SCRATCH)/*.OrhtoD
	rm  -rf $(MONDO_SCRATCH)/*.Kxc
	rm  -rf $(MONDO_SCRATCH)/*
#
#       A RECURSIVELY GZIPED, DATE-TAGGED TARBALL
#   
backup:	purge
	cd $(MONDO_HOME)/.. ;\
        mv $(MONDO_HOME) MONDO_`date '+%B'`_`date '+%d'`_`date +%y` ;\
        tar -cvf MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar           \
                    MONDO_`date '+%B'`_`date '+%d'`_`date +%y`              ;\
        gzip     MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar          ;\
        mv  MONDO_`date '+%B'`_`date '+%d'`_`date +%y` $(MONDO_HOME) 
#
#---------------------------------------------------------------------------------