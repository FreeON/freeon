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

JObjs =  MultipoleLib.o MMoments.o JBlok.o MakeDJ.o 

all:	MD QCTC

clean:	CQCTC

purge:	CMD CQCTC 

MD:	
	make -i -C VMD

QCTC:	$(JObjs)
	$(LD) $(LDFLAGS) QCTC $(JObjs) $(LibOpts) $(F90_LIBRARIES) -lMcMurchieDavidson 
	mv QCTC $(MONDO_HOME)/Exec/QCTC

CQCTC:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/MakeDJ
	ls 

CMD:	
	make -i -C VMD clean


MakeDJ.o: 	JBlok.o MMoments.o MultipoleLib.o $(MONDO_HOME)/Libs/libMondoMods.a
JBlok.o:	MMoments.o VMD/VMDBlok.Inc        $(MONDO_HOME)/Libs/libMondoMods.a
MMoments.o:	                                  $(MONDO_HOME)/Libs/libMondoMods.a
