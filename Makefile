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
all:	CatCpy mm s e 1 2 t  #d #v 
#
clean:	cmm cs ce c1 c2 cd cv ct 
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
#
purge:	pmm ps pe p1 p2 pd pv pt pMisc  
#
release: clean rmLegacy swREADME rmm rs re r1 r2 rd rv rt tarball
#
backup:	purge rmLegacy tarball
#
rmLegacy:
	rm -rf DirectJ
	rm -rf ONX
	rm -rf XForm
	rm -rf Inpts
	rm -rf Scratch
	rm -rf Exec
	rm -rf Libs
	rm -rf QuNew
	rm -rf MondoSCF
	rm -rf Dynamo
	rm -rf Bader
	rm -rf Docs
	rm -rf LibTwoE
	rm -rf IntCoo
	rm -rf ChangeLog*
#
swREADME:
	rm -f  README; cp $(MONDO_HOME)/Includes/ALPHA_RELEASE README
#
CatCpy:	
	cat $(MONDO_HOME)/Includes/CopyrightNotice.txt
	sleep 1
	cat $(MONDO_HOME)/README
	sleep 2
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
	$(MAKE)    -C MondoSCF2
rs:	
	$(MAKE) -i -C MondoSCF2 release
ps:	
	$(MAKE) -i -C MondoSCF2 purge
cs:	
	$(MAKE) -i -C MondoSCF2 clean
#----------------------------------------------
#    Solving SCF Equations
 e:	
	$(MAKE)    -C SCFeqs
re:	
	$(MAKE) -i -C SCFeqs release
pe:	
	$(MAKE) -i -C SCFeqs purge
ce:	
	$(MAKE) -i -C SCFeqs clean
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
 2:	q   h  o2
r2:	rq rh  ro2
c2:	cq ch  co2
p2:	ph pq  po2
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
#     ONX2
 o2:	
	$(MAKE)    -C ONX2
ro2:	
	$(MAKE) -i -C ONX2 release
po2:	
	$(MAKE) -i -C ONX2 purge
co2:	
	$(MAKE) -i -C ONX2 clean 
#----------------------------------------------
#     TwoE
 t:	
	$(MAKE)    -C TwoE
rt:	
	$(MAKE) -i -C TwoE release
pt:	
	$(MAKE) -i -C TwoE purge
ct:	
	$(MAKE) -i -C TwoE clean 
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
#     DX Visualization support
 d:	
	$(MAKE)    -C VisDX
rd:	
	$(MAKE) -i -C VisDX release
pd:	
	$(MAKE) -i -C VisDX purge
cd:	
	$(MAKE) -i -C VisDX clean
#----------------------------------------------
#     Validation suite
 v:	
	$(MAKE)    -C Validate 
vp:	
	$(MAKE)    -C Validate pbs
rv:	
	$(MAKE) -i -C Validate release
pv:	
	$(MAKE) -i -C Validate purge
cv:	
	$(MAKE) -i -C Validate clean
#----------------------------------------------
#   Cleaning of other directories 
#
pMisc:	pInp #pScr pPWD
#
pInp:	
	$(MAKE) -i -C $(MONDO_HOME)/Benchmarks/PROTEINS purge
	$(MAKE) -i -C $(MONDO_HOME)/Benchmarks/WATER purge
	$(MAKE) -i -C $(MONDO_HOME)/Validate purge
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
        tar -cvf MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar  \
                    MONDO_`date '+%B'`_`date '+%d'`_`date +%y`     ;\
	gzip MondoSCF_`date '+%B'`_`date '+%d'`_`date +%y`.tar     ;\
        mv  MONDO_`date '+%B'`_`date '+%d'`_`date +%y` $(MONDO_HOME) 
#---------------------------------------------------------------------------------
