###################################################################################
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

all: NONX XForce

Objects1 = ONXParameters.o ONXInit.o ONXGet.o ONXFillOut.o ONXMemory.o        \
           ONXCtrSclg.o ONXVrr.o ONXHrr.o ONXPut.o ONXRGen.o Sizes.o          \
           Stats.o ONXInLoop.o ONXComptK.o ONXRng.o ONXDOrder.o               \
           QuickSortDis.o PartDrv.o ONXMemInit.o CCDriver.o VRRDriver.o ONX.o

Objects2 = ONXParameters.o ONXMemory.o ContractionScaling.o MatFilter.o         \
           Sizes.o Stats.o MemInit.o InitExchangeMatrix.o GradONX.o Contract.o  \
           HrrBra.o HrrKet.o VRRl.o RGen1C.o CCDriver.o VRRDriver.o             \
           DisOrderGrad.o RangeOfDensity.o VRRs.o QuickSortDis.o PutDis.o       \
           GDLoader.o GDrivers.o ComputeXForce.o ContractG.o RGen.o GetSubBlk.o \
           GetGradient.o InitSubInd.o DigestGradient.o GetAdr.o GetIntSpace.o

ReleaseFiles= \
           ONXParameters.x ONXInit.x ONXGet.x ONXFillOut.x ONXMemory.x          \
           ONXCtrSclg.x ONXVrr.x ONXHrr.x ONXPut.x ONXRGen.x Sizes.x            \
           Stats.x ONXInLoop.x ONXComptK.x ONXRng.x ONXDOrder.x                 \
           QuickSortDis.x PartDrv.x ONXMemInit.x CCDriver.x VRRDriver.x ONX.x   \
           ONXParameters.x ONXMemory.x ContractionScaling.x MatFilter.x         \
           Sizes.x Stats.x MemInit.x InitExchangeMatrix.x GradONX.x Contract.x  \
           HrrBra.x HrrKet.x VRRl.x RGen1C.x CCDriver.x VRRDriver.x             \
           DisOrderGrad.x RangeOfDensity.x VRRs.x QuickSortDis.x PutDis.x       \
           GDLoader.x GDrivers.x ComputeXForce.x ContractG.x RGen.x GetSubBlk.x \
           GetGradient.x InitSubInd.x DigestGradient.x GetAdr.x GetIntSpace.x
#
Libs =  $(F90_LIBRARIES) $(LAPACK) -lVRR
#
NONX:	$(Objects1)
	$(LD) $(LDFLAGS) NONX $(Objects1) $(Libs)
	mv NONX $(MONDO_EXEC)/ONX
#
XForce: $(Objects2)
	$(LD) $(LDFLAGS) XForce $(Objects2) $(Libs)
	mv XForce $(MONDO_EXEC)/XForce
#
release:purge $(ReleaseFiles)
	rm -f HRR.tar
	rm -f Notes.txt
	rm -f README
#
purge:	clean
	rm -f $(REMOVESRCE)
#
clean:	
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/ONX
	rm -f $(MONDO_EXEC)/XForce
	ls
#
doc.txt:
	cat *.F90 | grep "^!H" | sed "s/^!H/  /" > doc.txt;


