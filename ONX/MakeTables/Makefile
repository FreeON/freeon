include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions 

LD=$(F90)
Objects = Convert.o GDump.o EDump.o HDump.o GetGammaTable.o \
          ConvertG0.o  ConvertG1.o  ConvertG2.o  ConvertG3.o ConvertG4.o \
          ConvertG5.o  ConvertG6.o  ConvertG7.o  ConvertG8.o ConvertG9.o \
          ConvertG10.o ConvertG11.o ConvertG12.o \
          ConvertE.o 

makes:	$(Objects)
	$(LD) $(LDFLAGS) ConvertTables $(Objects) $(LibOpts)
	./ConvertTables

clean:
	rm -f $(REMOVEALL)
	rm -f ConvertTables
	rm -f $(MONDO_HOME)/MondoMods/MMA/Functions/*.bin
	ls

LinK.o:	$(MONDO_HOME)/Includes/Wrk.Inc

