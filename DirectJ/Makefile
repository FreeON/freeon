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

EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions  \
                -I$(MONDO_HOME)/DirectJ/MMA/Moments

JObjs =  Moments.o MMML_OLD.o MMML.o Multipoles.o JBlok.o MakeDJ.o 

all:	MD DJ

clean:	CDJ

purge:	CMD CDJ

MD:	
	$(MAKE) -i -C VMD

DJ:	$(JObjs)
	$(LD) $(LDFLAGS) DirectJ $(JObjs) $(LibOpts) $(F90_LIBRARIES) -lMcMurchieDavidson 
	mv DirectJ $(MONDO_HOME)/Exec/DirectJ

CDJ:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/MakeDJ
	ls 

CMD:	
	$(MAKE) -i -C VMD clean


MakeDJ.o: 	  JBlok.o Multipoles.o                $(MONDO_HOME)/Libs/libMondoMods.a
JBlok.o:	  Multipoles.o VMD/VMDBlok.Inc        $(MONDO_HOME)/Libs/libMondoMods.a
MultipoleLib.o:   Moments.o  MMML_OLD.o MMML.o        $(MONDO_HOME)/Libs/libMondoMods.a
Moments.o:	                                      $(MONDO_HOME)/Libs/libMondoMods.a

