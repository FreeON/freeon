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
#    MAIN MAKEFILE FOR MONDOSCF
##############################################################################
#
include $(MONDO_HOME)/Includes/RemoveAll
#
all:	CatCpy mm s e x n 1 2 
#
clean:	cmm cs ce cx cn c1 c2
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
#
purge:	pmm ps pe px pn p1 p2 pMisc
#
release: rmLegacy rmm rs re rx rn r1 r2 tarball
#
backup:	purge rmLegacy tarball
#
rmLegacy:
	rm -rf DirectJ
	rm -rf ONX/GMEONX
	rm -rf ONX/GSONX
	rm -rf ONX/MASONX
	rm -rf ONX/PONX*
	rm -rf ONX/SONX
	rm -rf Inpts/DEVELOP
	rm -rf Inpts/PROTEINS
	rm -rf Inpts/HMX
#
CatCpy:	
	cat $(MONDO_HOME)/Includes/CopyrightNotice.txt
	sleep 1
	cat $(MONDO_HOME)/README
	sleep 30
#----------------------------------------------
#   MondoMods
 mm:	
	$(MAKE)    -C MondoMods all
rmm:
	$(MAKE) -i -C MondoMods release
pmm:	
	$(MAKE) -i -C MondoMods purge
cmm:	
	$(MAKE) -i -C MondoMods clean
#----------------------------------------------
#   MondoSCF
 s:	
	$(MAKE)    -C MondoSCF
rs:	
	$(MAKE) -i -C MondoSCF release
ps:	
	$(MAKE) -i -C MondoSCF purge
cs:	
	$(MAKE) -i -C MondoSCF clean
#----------------------------------------------
#    Solving SCF Equantions
 e:	
	$(MAKE)    -C SCFeqs
re:	
	$(MAKE) -i -C SCFeqs release
pe:	
	$(MAKE) -i -C SCFeqs purge
ce:	
	$(MAKE) -i -C SCFeqs clean
#----------------------------------------------
#    Orthogonalization transformations
 x:	
	$(MAKE)    -C XForm
rx:	
	$(MAKE) -i -C XForm release
px:	
	$(MAKE) -i -C XForm purge
cx:	
	$(MAKE) -i -C XForm clean
#----------------------------------------------
#     Geometry optimization via Quasi Newton
 n:	
	$(MAKE)    -C QuNew
rn:	
	$(MAKE) -i -C QuNew release
pn:	
	$(MAKE) -i -C QuNew purge
cn:	
	$(MAKE) -i -C QuNew clean
#----------------------------------------------
#      One electron routines
 1:
	$(MAKE)    -C OneE
r1:	
	$(MAKE) -i -C OneE release
p1:	
	$(MAKE) -i -C OneE purge
c1:	
	$(MAKE) -i -C OneE clean
#----------------------------------------------
#     Two electron directories
#
 2:	q   h  o
r2:	rq rh  ro
c2:	cq ch  co
p2:	ph pq  po
#----------------------------------------------
#     ONX
 o:	
	$(MAKE)    -C ONX
ro:	
	$(MAKE) -i -C ONX release
po:	
	$(MAKE) -i -C ONX purge
co:	
	$(MAKE) -i -C ONX clean 
#----------------------------------------------
#     QCTC 
 q:	
	$(MAKE)    -C QCTC
rq:	
	$(MAKE) -i -C QCTC release
pq:	
	$(MAKE) -i -C QCTC purge
cq:	
	$(MAKE) -i -C QCTC clean
#----------------------------------------------
#      HiCu
 h:	
	$(MAKE) -C HiCu
rh:	
	$(MAKE) -i -C HiCu release
ph:	
	$(MAKE) -i -C HiCu purge
ch:	
	$(MAKE) -i -C HiCu clean
#----------------------------------------------
#   Cleaning of other directories 
#
pMisc:	pLib pWrk pInp pScr pPWD
#
pLib:	
	rm -rf $(MONDO_HOME)/Libs/*
#
pExec:	
	rm -rf $(MONDO_HOME)/Exec/*
#
pWrk:	
	$(MAKE) -i -C $(MONDO_WORK) clean
#
pInp:	
	$(MAKE) -i -C $(MONDO_HOME)/Inpts/WATER    purge
	$(MAKE) -i -C $(MONDO_HOME)/Inpts/Validate purge
#
pScr:	
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
	rm  -rf $(MONDO_SCRATCH)/*.PFFT
	rm  -rf $(MONDO_SCRATCH)/*
#
pPWD:
	rm -f $(REMOVEALL)
#---------------------------------------------------------------------------------
#     Date tagged tarball 
tarball:
	cd $(MONDO_HOME)/.. ;\
        mv $(MONDO_HOME) MONDO_`date '+%B'`_`date '+%d'`_`date +%y` ;\
        tar -zcvf MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tgz  \
                     MONDO_`date '+%B'`_`date '+%d'`_`date +%y`     ;\
        mv  MONDO_`date '+%B'`_`date '+%d'`_`date +%y` $(MONDO_HOME) 
#---------------------------------------------------------------------------------
