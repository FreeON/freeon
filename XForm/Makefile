#
#     Matt Challacombe
#     Los Alamos National Laboratory
#     Copyright 2000, The University of California
#
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

CPPOPTS=-DEXTREME_DEBUG -DEXTREME_PRINT_DEBUG

Libs =  $(F90_LIBRARIES) $(LAPACK)

BlkAInvObjs = BlokAInv.o  
LowdinOObjs = LowdinO.o 

All:	Blk Lwd

Blk:	$(BlkAInvObjs)
	$(LD) $(LDFLAGS) AInv $(BlkAInvObjs) $(LibOpts) $(Libs)
	mv AInv $(MONDO_HOME)/Exec/AInv

Lwd:	$(LowdinOObjs)
	$(LD) $(LDFLAGS) LowdinO $(LowdinOObjs) $(LibOpts) $(Libs)
	mv LowdinO $(MONDO_HOME)/Exec/LowdinO

clean:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/AInv
	rm -f $(MONDO_HOME)/Exec/LOrthog
	ln -s /dev/null core
	ls -l 

LowdinO.o:	$(MONDO_HOME)/Libs/libMondoMods.a
BlokAInv.o:	$(MONDO_HOME)/Libs/libMondoMods.a
