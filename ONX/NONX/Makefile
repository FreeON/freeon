include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

Objects = ONXParameters.o ONXMemory.o ContractionScaling.o MatFilter.o \
          InitExchangeMatrix.o Sizes.o Stats.o MemInit.o \
          ONX.o Contract.o hrrbra_old.o hrrket_old.o \
          vrrl.o RGen1C.o CCDriver.o VRRDriver.o \
          DisOrder.o RangeOfDensity.o VRRs.o


Libs =  $(F90_LIBRARIES) $(LAPACK) 

NONX:	$(Objects)
	$(LD) $(LDFLAGS) NONX $(Objects) $(LibOpts) $(Libs)
	mv NONX $(MONDO_HOME)/Exec/NONX

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/NONX
	ls
	
