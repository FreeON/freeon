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

all: NONX XForce

Objects1 = ONXParameters.o ONXMemory.o ContractionScaling.o MatFilter.o \
           InitExchangeMatrix.o Sizes.o Stats.o MemInit.o \
           ONX.o Contract.o HRRBra.o HRRKet.o \
           VRRl.o RGen1C.o CCDriver.o VRRDriver.o \
           DisOrder.o RangeOfDensity.o VRRs.o QuickSortDis.o PutDis.o \
           RangeOfExchange.o Tables.o GetTables.o \
           GetSubBlk.o XPose.o InitSubInd.o \
           RGen.o Digest.o Scatter.o GetAdr.o PutSubBlk.o Fillout.o \
           ComputeKe.o ComputeKg.o GetIntCode.o GetIntSpace.o

Objects2 = ONXParameters.o ONXMemory.o ContractionScaling.o        \
           Sizes.o Stats.o MemInit.o GradONX.o Contract.o          \
           HRRBra.o HRRKet.o VRRl.o RGen1C.o CCDriver.o            \
           VRRDriver.o DisOrderGrad.o RangeOfDensity.o VRRs.o      \
           QuickSortDis.o PutDis.o Tables.o GetTables.o GDLoader.o \
           GDrivers.o ComputeXForce.o ContractG.o RGen.o           \
           HRRBraGrad.o HRRKetGrad.o GetSubBlk.o GetGradient.o     \
           InitSubInd.o DigestGradient.o GetAdr.o

ReleaseFiles=ONXParameters.x ONXMemory.x ContractionScaling.x MatFilter.x \
             InitExchangeMatrix.x Sizes.x Stats.x MemInit.x \
             ONX.x GradONX.x Contract.x HRRBra.x HRRKet.x \
             VRRl.x RGen1C.x CCDriver.x VRRDriver.x \
             DisOrder.x RangeOfDensity.x VRRs.x QuickSortDis.x PutDis.x \
             RangeOfExchange.x Tables.x GetTables.x \
             GetSubBlk.x XPose.x InitSubInd.x \
             RGen.x Digest.x Scatter.x GetAdr.x PutSubBlk.x Fillout.x \
             ComputeKe.x ComputeKg.x GetGradient.x     \
             InitSubInd.x DigestGradient.x GetAdr.x GDLoader.x \
             BSearch.x  HRRBraGrad.x  HRRKetGrad.x  LSearch.x \
             ComputeXForce.x  ContractG.x DisOrderGrad.x GDrivers.x 
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
