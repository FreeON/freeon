include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

all: MakeGD

LD=$(F90)
Objects = MakeGD.o GCMaker.o 

MakeGD:	$(Objects)
	$(LD) $(LDFLAGS) MakeGD $(Objects) $(LibOpts) 
	./MakeGD

clean:
	rm -f $(REMOVEALL)
	rm -f $(MONDO_HOME)/ONX/GDrivers/MakeGD
	rm -f $(MONDO_HOME)/ONX/GDrivers/*.binary

