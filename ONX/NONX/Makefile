############################################################################
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
############################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
Objects = ONXParameters.o ONXMemory.o ContractionScaling.o MatFilter.o \
          InitExchangeMatrix.o Sizes.o Stats.o MemInit.o \
          ONX.o Contract.o HRRBra.o HRRKet.o \
          vrrl.o RGen1C.o CCDriver.o VRRDriver.o \
          DisOrder.o RangeOfDensity.o VRRs.o QuickSortDis.o PutDis.o \
          RangeOfExchange.o Tables.o GetTables.o \
          ComputeK.o GetSubBlk.o XPose.o \
          RGen.o Digest.o Scatter.o GetAdr.o PutSubBlk.o Fillout.o
#
ReleaseFiles=ONXParameters.x ONXMemory.x ContractionScaling.x MatFilter.x \
             InitExchangeMatrix.x Sizes.x Stats.x MemInit.x \
             ONX.x Contract.x HRRBra.x HRRKet.x \
             vrrl.x RGen1C.x CCDriver.x VRRDriver.x \
             DisOrder.x RangeOfDensity.x VRRs.x QuickSortDis.x PutDis.x \
             RangeOfExchange.x Tables.x GetTables.x \
             ComputeK.x GetSubBlk.x XPose.x \
             RGen.x Digest.x Scatter.x GetAdr.x PutSubBlk.x Fillout.x
#
Libs =  $(F90_LIBRARIES) $(LAPACK) 
#
NONX:	$(Objects)
	$(LD) $(LDFLAGS) NONX $(Objects) $(LibOpts) $(Libs)
	mv NONX $(MONDO_HOME)/Exec/NONX
#
release:clean $(RELEASEFILES)
#
purge:	clean
	rm -f $(REMOVESRCE)
#
clean:	
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/NONX
	ls


