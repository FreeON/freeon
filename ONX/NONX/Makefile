include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

all: NONX XForce

Objects1 = ONXParameters.o ONXMemory.o ContractionScaling.o MatFilter.o \
           InitExchangeMatrix.o Sizes.o Stats.o MemInit.o \
           ONX.o Contract.o HRRBra.o HRRKet.o \
           vrrl.o RGen1C.o CCDriver.o VRRDriver.o \
           DisOrder.o RangeOfDensity.o VRRs.o QuickSortDis.o PutDis.o \
           RangeOfExchange.o Tables.o GetTables.o \
           GetSubBlk.o XPose.o InitSubInd.o \
           RGen.o Digest.o Scatter.o GetAdr.o PutSubBlk.o Fillout.o \
           ComputeKe.o ComputeKg.o

Objects2 = ONXParameters.o ONXMemory.o ContractionScaling.o        \
           Sizes.o Stats.o MemInit.o GradONX.o Contract.o          \
           HRRBra.o HRRKet.o vrrl.o RGen1C.o CCDriver.o            \
           VRRDriver.o DisOrderGrad.o RangeOfDensity.o VRRs.o      \
           QuickSortDis.o PutDis.o Tables.o GetTables.o GDLoader.o \
           GDrivers.o ComputeXForce.o ContractG.o RGen.o           \
           HRRBraGrad.o HRRKetGrad.o GetSubBlk.o GetGradient.o     \
           InitSubInd.o DigestGradient.o GetAdr.o

Libs =  $(F90_LIBRARIES) $(LAPACK) -lVRR

NONX:	$(Objects1)
	$(LD) $(LDFLAGS) NONX $(Objects1) $(LibOpts) $(Libs)
	mv NONX $(MONDO_HOME)/Exec/ONX

XForce: $(Objects2)
	$(LD) $(LDFLAGS) XForce $(Objects2) $(LibOpts) $(Libs)
	mv XForce $(MONDO_HOME)/Exec/XForce

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/ONX
	rm -f $(MONDO_HOME)/Exec/XForce
	ls
