include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions 

all: ConvertTables 

LD=$(F90)
Objects1 = Convert.o GDump.o EDump.o HDump.o ADump.o GetGammaTable.o \
           ConvertG0.o  ConvertG1.o  ConvertG2.o  ConvertG3.o ConvertG4.o \
           ConvertG5.o  ConvertG6.o  ConvertG7.o  ConvertG8.o ConvertG9.o \
           ConvertG10.o ConvertG11.o ConvertG12.o \
           ConvertE.o 

Objects2 = MakeGD.o GCMaker.o

ConvertTables:	$(Objects1)
	$(LD) $(LDFLAGS) ConvertTables $(Objects1) $(LibOpts)
	./ConvertTables

MakeGDrivers:	$(Objects2)
	$(LD) $(LDFLAGS) MakeGDrivers $(Objects2) $(LibOpts)
	./MakeGDrivers

clean:
	rm -f $(REMOVEALL)
	rm -f ConvertTables MakeGDrivers
	rm -f $(MONDO_HOME)/MondoMods/MMA/Functions/*.bin
	ls

