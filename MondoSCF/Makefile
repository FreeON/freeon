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
#
ObjectFiles = SCFLocals.o            \
              ParsingKeys.o          \
              Spawn.o                \
              GetPID.o               \
              Overlay.o              \
              ParseInput.o           \
              SetSCFs.o              \
              DrvSCFs.o              \
              DrvFrcs.o              \
              ChkSCFs.o              \
              MondoSCF.o
#              MondoSCF_Grads.o

all:	$(ObjectFiles) 
	$(LD) $(LDFLAGS) MondoSCF $(ObjectFiles) $(LibOpts) $(F90_LIBRARIES)
	mv MondoSCF $(MONDO_HOME)/Exec/MondoSCF


Test:	$(TestObjs) 
	$(LD) $(LDFLAGS) Test $(TestObjs) $(LibOpts) $(Libs)
	mv Test $(MONDO_HOME)/Exec/Test


MondoSCF.o:	DrvSCFs.o ChkSCFs.o SetSCFs.o ParseInput.o Overlay.o $(MONDO_HOME)/Libs/libMondoMods.a
ChkSCFs.o:	SCFLocals.o                                          $(MONDO_HOME)/Libs/libMondoMods.a
DrvSCFs.o:	SCFLocals.o Overlay.o                                $(MONDO_HOME)/Libs/libMondoMods.a
ParseInput.o:	SCFLocals.o                                          $(MONDO_HOME)/Libs/libMondoMods.a
SetSCFs.o:	SCFLocals.o Overlay.o                                $(MONDO_HOME)/Libs/libMondoMods.a

clean:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/MondoSCF
	ln -s /dev/null core
	ls -l
