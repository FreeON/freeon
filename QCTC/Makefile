####################################################################################
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
####################################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =-DEXPLICIT_SOURCE
#
EXTRA_INCLUDES=-I$(MONDO_HOME)/OneE/ -I LibQCTC/
#
QCTCLibs=$(F90_LIBRARIES) -lQCTC
#
CommonObjs=XLate.o CTraX.o Globals.o MondoPoles.o Thresholds.o    \
           PoleTree.o PFFTensor.o PBCFarField.o TreeWalk.o 
#
QCTCObjs=$(CommonObjs) JGen.o NukE.o QCTC.o 
JFrcObjs=$(CommonObjs) BlokTrPdJ.o JForce.o
PFFTObjs=$(CommonObjs) MakePFFT.o
#
ReleaseFiles=XLate.x CTraX.x Globals.x MondoPoles.x \
             Thresholds.x PoleTree.x PFFTensor.x    \
             PBCFarField.x TreeWalk.x JGen.x NukE.x \
             QCTC.x BlokTrPdJ.x JForce.x
#
all:	lQCTC QCTC JFrc PFFT
#
release:clean $(ReleaseFiles)
	rm -rf MMA/Moments
	rm -rf MMA/PAC
	rm -rf MMA/MD3TRR
	rm -rf MMA/XLate
#
clean:	clQCTC cQCTC cJFrc cMMA cPFFT
#
purge:	plQCTC cQCTC cJFrc pMMA cPFFT
	rm -f $(REMOVESRCE)
#
lQCTC:	
	$(MAKE) -C LibQCTC
plQCTC:
	$(MAKE) -C LibQCTC purge
clQCTC:
	$(MAKE) -C LibQCTC clean
rlQCTC:
	$(MAKE) -C LibQCTC release
#
QCTC:	$(QCTCObjs)
	$(LD) $(LDFLAGS) QCTC $(QCTCObjs) $(QCTCLibs) 
	mv QCTC $(MONDO_EXEC)/QCTC
cQCTC:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/QCTC
	ls 
#
JFrc:	$(JFrcObjs)
	$(LD) $(LDFLAGS) JForce $(JFrcObjs) $(QCTCLibs) 
	mv JForce $(MONDO_EXEC)/JForce
cJFrc:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/JForce
	ls 
#
PFFT:	$(PFFTObjs)
	$(LD) $(LDFLAGS) MakePFFT $(PFFTObjs) $(QCTCLibs) 
	mv MakePFFT $(MONDO_EXEC)/MakePFFT
cPFFT:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/MakePFFT
	ls 
#
MMA:	ctrax hgtosp hgtrax 
#
cMMA:	cctrax chgtosp chgtrax
#
pMMA:	cMMA
#
ctrax:
	$(MAKE) -i -C MMA/CTraX
cctrax:
	$(MAKE) -i -C MMA/CTraX  clean
#	
hgtosp:
	$(MAKE) -i -C MMA/HGToSP
chgtosp:
	$(MAKE) -i -C MMA/HGToSP clean
#
hgtrax:
	$(MAKE) -i -C MMA/HGTraX
chgtrax:
	$(MAKE) -i -C MMA/HGTraX clean
#
# dependencies are laking here, needs work!
#
