#================================================================================
#
#     Matt Challacombe
#     Los Alamos National Laboratory
#     Copyright 2000, The University of California
#
#================================================================================
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
include VMDObjs 

EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions 

all:	$(VObjectFiles)
	$(AR) $(ARFLAGS) MD.a $(?:.f=.o)
	$(RANLIB) MD.a 
	mv MD.a $(MONDO_HOME)/Libs/libMcMurchieDavidson.a

vmd:	expVMD.o
	$(LD) $(LDFLAGS) vmd expVMD.o
	./vmd

smd:	expSMD.o
	$(LD) $(LDFLAGS) smd expSMD.o
	./smd

clean: 
	rm -f $(REMOVEMISC) $(REMOVEMODS) 
	rm -f vmd smd
	rm -f $(MONDO_HOME)/Libs/libMcMurchieDavidson.a
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
