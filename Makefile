#========================================================================
#
#  This makefile is part of the MondoSCF suite of 
#  linear scaling electronic structure codes.  
#
#  Matt Challacombe
#  Los Alamos National Laboratory
#  Copywrite 2000, The University of California
#
#========================================================================
include $(MONDO_HOME)/Includes/RemoveAll
#
all:	 Env        mm  s  x  1  2  e
#
clean:	CExec cmm cs cx c1 c2 ce 
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
#
purge:	clean PScr PWrk pmm p2 
#
Env:	
	env | grep -i "mondo_"
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
	$(MAKE) -C MondoMods
#
cmm:	
	$(MAKE) -i -C MondoMods clean
#
pmm:	
	$(MAKE) -i -C MondoMods purge
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
        mv $(MONDO_HOME) $(MONDO_HOME)_`date '+%B'`_`date '+%d'`_`date +%y` ;\
        tar -cvf MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar           \
                    MONDO_`date '+%B'`_`date '+%d'`_`date +%y`              ;\
        gzip     MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar          ;\
        mv  $(MONDO_HOME)_`date '+%B'`_`date '+%d'`_`date +%y` $(MONDO_HOME) 
#
