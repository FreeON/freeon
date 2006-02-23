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
CPPMISC = -DOVERLAP_PROJECTION # -DPAPIEX  -DPAC_DEBUG  #-DMAC_DEBUG # -DEXPLICIT_SOURCE 
#

EXTRA_INCLUDES=-I$(MONDO_HOME)/OneE/ -I$(MONDO_HOME)/MondoMods/MMA/Functions #-I./LibQCTC/
#
QCTCLibs=$(F90_LIBRARIES) $(LAPACK)  #-lQCTC 
#
CommonObjs=PoleNodeType.o Indexing.o MMML.o  CTraX.o \
CTraX4.o CTraX6.o CTraX8.o CTraX10.o  CTraX12.o  CTraX14.o  \
HGTraX0.o  HGTraX1.o  HGTraX2.o  HGTraX3.o HGTraX4.o \
Globals.o Gamma.o KTrax.o MD3TRR.o AA95.o AA95c.o MondoPoles.o Thresholds.o    \
PoleTree.o TreeWalk.o PFFTensor.o DivPFFTensor.o PBCFarField.o 
#
QCTCObjs=$(CommonObjs) JGen.o NukE.o QCTC.o 
JFrcObjs=$(CommonObjs) BlokTrPdJ.o JForce.o
PFFTObjs=$(CommonObjs) MakePFFT.o
#
ReleaseFiles=
#
#all:	lQCTC QCTC JFrc PFFT
all:	QCTC PFFT
#

tarball:
	tar -cvf JWALKING_`date '+%B'`_`date '+%d'`_`date +%y`.tar  \
		*.F* Makefile CTraX HGToSP MACError figs* *.Inc *.m  	         
	gzip JWALKING_`date '+%B'`_`date '+%d'`_`date +%y`.tar     ;\

release: clean $(ReleaseFiles) rlQCTC 

#
clean:	cQCTC cJFrc cPFFT
#
purge:	cQCTC cJFrc cPFFT
	rm -f $(REMOVESRCE)
#
QCTC: $(QCTCObjs)
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
#--------------------------------------------------------------------
#	RULES FOR FORTRAN77 GENERATION
#--------------------------------------------------------------------
#	
77:	MkG MkK	MkR
c77:	cG  cK  cR
#
MkG:	cG VAuxWriter.o VOvrWriter.o
	$(LD) $(LDFLAGS) VAG VAuxWriter.o $(F90_LIBRARIES)
	./VAG
	rm ./VAG
	$(LD) $(LDFLAGS) OVR VOvrWriter.o $(F90_LIBRARIES)
	./OVR
	rm ./OVR
cG:	
	rm -f VAuxWriter.f90 VAuxWriter.o Gamma.F

MkH:	cH HGWriter.o QSortI.o
	$(LD) $(LDFLAGS) VH HGWriter.o QSortI.o $(F90_LIBRARIES)
	./VH
	rm ./VH
cH:	
	rm -f HGWriter.f90 HGWriter.o HGTraX*.F 



MkK:	cK VKWriter.o QSortI.o
	$(LD) $(LDFLAGS) VK VKWriter.o QSortI.o $(F90_LIBRARIES)
	./VK
	rm ./VK
cK:	
	rm -f VKWriter.f90 VKWriter.o KTrax.F


MkR:	cR VRWriter.o
	$(LD) $(LDFLAGS) VR VRWriter.o $(F90_LIBRARIES)
	./VR
	rm ./VR
cR:	
	rm -f VRWriter.f90 VRWriter.o MD3TRR.F

#
QCTC.o:		JGen.o PoleTree.o NukE.o PoleNodeType.o 
JGen.o:		TreeWalk.o PoleNodeType.o 
NukE.o:		TreeWalk.o PoleNodeType.o 	
TreeWalk.o:	MondoPoles.o CTraX.Inc	
PoleTree.o:	MondoPoles.o PoleNodeType.o 
MondoPoles:	Indexing.o PoleNodeType.o 

