include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=
#
MiscObjs=ShellPairStruct.o ShellPrint.o Scratch.o GetAtomPair.o InitDbl.o 
#
ReleaseFiles=ShellPairStruct.x ShellPrint.x Scratch.x GetAtomPair.x InitDbl.x 
#
all:    Misc
#
clean:  CMisc
#
purge:clean
	rm -f $(REMOVESRCE)
#
release: clean $(ReleaseFiles)
#
Misc:	$(MiscObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libTwoEMisc.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libTwoEMisc.a
#
CMisc:	
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
#
