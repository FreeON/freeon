include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

Objects = ONXParameters.o ContractionScaling.o MatFilter.o \
          InitExchangeMatrix.o IDmn.o NFinal.o iT.o \
          ONX.o contract.o hrrbra_old.o hrrket_old.o \
          vrrl.o GetBigInt.o RGen1C.o CCDriver.o VRRDriver.o


Libs =  $(F90_LIBRARIES) $(LAPACK) 

ONX:	$(Objects)
	$(LD) $(LDFLAGS) ONX $(Objects) $(LibOpts) $(Libs)
	mv ONX $(MONDO_HOME)/Exec/ONX

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/ONX
	ls
	
