include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

Objects = ONXParameters.o ContractionScaling.o InitExchangeMatrix.o \
          MinMax.o MatFilter.o ONX.o ONXBuf.o ONXSch.o RLoader.o   \
          MMML.o Misc.o hgtftosp.o onxsort2.o AddDis.o GetID.o     \
          onxio.o RangeOfExchange.o ONXLoop.o CCDriver.o           \
          VRRDriver.o ErrorF.o BuffTbl.o GetSubBlk.o CalcK.o       \
          skipout.o BSearch.o LSearch.o $(MAIN__)                  \
          onx1111.o onx1121.o onx1122.o onx1161.o onx1162.o        \
          onx1166.o onx2111.o onx2121.o onx2122.o onx2161.o        \
          onx2162.o onx2211.o onx2221.o onx2222.o onx2261.o        \
          onx6111.o onx6121.o onx6122.o onx6161.o onx6211.o        \
          onx6221.o onx6611.o vrrl.o                               \
          hrrket.o hrrbra.o digest.o tgen0.o tgen1t2.o vrrs.o      \
          tgen5t6.o tgen7t8.o contract.o tgen3t4.o getadr.o \
          putsubblk.o hrrket_old.o hrrbra_old.o           \
          MondoHalt.o OpenAppend.o LineEnd.o            \
          RangeOfDensity.o GatherK.o XPose.o Fillout.o 


Libs =  $(F90_LIBRARIES) $(LAPACK) -lPONXMults

ONX:	$(Objects)
	$(LD) $(LDFLAGS) ONX $(Objects) $(LibOpts) $(Libs)
	mv ONX $(MONDO_HOME)/Exec/ONX

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/ONX
	ls
	
