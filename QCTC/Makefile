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
EXTRA_INCLUDES= -I$(MONDO_HOME)/QCTC/MMA/CTraX   \
                -I$(MONDO_HOME)/QCTC/MMA/HGToSP  \
                -I$(MONDO_HOME)/QCTC/MMA/HGTraX  \
                -I$(MONDO_HOME)/OneE/
#
CommonObjs=XLate.o CTraX.o Globals.o MondoPoles.o Thresholds.o  \
           PoleTree.o PFFTensor.o PBCFarField.o TreeWalk.o 
QCTCObjs=$(CommonObjs) JGen.o NukE.o QCTC.o 
JFrcObjs=$(CommonObjs) BlokTrPdJ.o JForce.o
PFFTObjs=$(CommonObjs) MakePFFT.o
#
ReleaseFiles=XLate.x CTraX.x Globals.x MondoPoles.x \
             Thresholds.x PoleTree.x PFFTensor.x    \
             PBCFarField.x TreeWalk.x JGen.x NukE.x \
             QCTC.x BlokTrPdJ.x JForce.x
#
all:	QCTC JFrc PFFT
#
release:clean $(ReleaseFiles)
	rm -rf MMA/Moments
	rm -rf MMA/PAC
	rm -rf MMA/MD3TRR
	rm -rf MMA/XLate
#
clean:	CQCTC CJFrc cmma 
#
purge:	CQCTC CJFrc cmma 
	rm -f $(REMOVESRCE)
#
QCTC:	$(QCTCObjs)
	$(LD) $(LDFLAGS) QCTC $(QCTCObjs) $(F90_LIBRARIES) 
	mv QCTC $(MONDO_EXEC)/QCTC
#
JFrc:	$(JFrcObjs)
	$(LD) $(LDFLAGS) JForce $(JFrcObjs) $(F90_LIBRARIES) 
	mv JForce $(MONDO_EXEC)/JForce
#
PFFT:	$(PFFTObjs)
	$(LD) $(LDFLAGS) MakePFFT $(PFFTObjs) $(F90_LIBRARIES) 
	mv MakePFFT $(MONDO_EXEC)/MakePFFT
#
CQCTC:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/QCTC
	ls 
#
CJFrc:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/JForce
	ls 
#
CPFFT:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/MakePFFT
	ls 
#
mma:	ctrax hgtosp hgtrax 
#
cmma:	cctrax chgtosp chgtrax
#
pmma:	cmma
	rm -f MMA/HGToSP/HGToSP.Inc
	rm -f MMA/CTraX/CTraX.Inc MMA/CTraX/IrRegular.Inc
	rm -f MMA/HGTraX/HGTraX.Inc
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
# dependencies are laking here.  This needs work!
#
QCTC.o:		MondoPoles.o TreeWalk.o CTraX.o Globals.o Thresholds.o  $(MONDO_LIB)/libMondoMods.a
MPoles.o:	MMA/HGToSP/HGToSP.Inc $(MONDO_LIB)/libMondoMods.a
TreeWalk.o:	MMA/HGTraX/HGTraX.Inc \
		MMA/CTraX/CTraX.Inc MMA/CTraX/IrRegulars.Inc $(MONDO_LIB)/libMondoMods.a
MakePFFT.o:	MondoPoles.o PFFTensor.o $(MONDO_LIB)/libMondoMods.a