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
	make -C MondoMods
#
cmm:	
	make -i -C MondoMods clean
#
pmm:	
	make -i -C MondoMods purge
#
#       FRONT END
#
s:	
	make -C MondoSCF
#
cs:	
	make -i -C MondoSCF clean
#
#       GENERATION OF ORTHOGONALIZATION MATRICES
#
x:	
	make -C XForm
#
cx:	
	make -i -C XForm clean
#
#       SOLVING THE SELF-CONSISTENT-FIELD EQUATIONS
#
e:	
	make -C SCFeqs
#
ce:	
	make -i -C SCFeqs clean
#
#       ONE ELECTRON ROUTINES
#
1:
	make -C OneE
#
c1:	
	make -i -C OneE clean
#
#       ORDER N EXCHANGE
#
o:	
	make -C ONX
#
co:	
	make -i -C ONX clean 
#
#       QUANTUM CHEMICAL TREE CODE
#
q:	
	make -C QCTC
#
cq:	
	make -C QCTC clean
#
pq:	
	make -C QCTC purge
#
#       HIERARCHICAL CUBATURE
#
h:	
	make -C HiCu
#
ch:	
	make -i -C HiCu clean
#
ph:	
	make -i -C HiCu purge
#
#       CLEAN EXECUTABLES
#
CExec:	
	rm -rf $(MONDO_HOME)/Exec/*
#
#       PURGE CURRENT WORK DIRECTORY
#
PWrk:	
	make -i -C $(MONDO_WORK) clean
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
#       MAKE A RECURSIVELY GZIPED, DATE-TAGGED TARBALL 
#   
backup:	purge
	cd $(MONDO_HOME)/.. ;\
        mv $(MONDO_HOME) $(MONDO_HOME)_`date '+%B'`_`date '+%d'`_`date +%y` ;\
        tar -cvf MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar           \
                    MONDO_`date '+%B'`_`date '+%d'`_`date +%y`              ;\
        gzip     MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar          ;\
        mv  $(MONDO_HOME)_`date '+%B'`_`date '+%d'`_`date +%y` $(MONDO_HOME) 
#
