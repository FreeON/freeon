include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
ObjectFiles = \
 mm_double_NN_1.o  mm_double_NT_1.o  mm_double_TN_1.o  mm_double_TT_1.o   \
 mm_double_NN_c.o  mm_double_NT_c.o  mm_double_TN_c.o  mm_double_TT_c.o  
#
all:	$(ObjectFiles)
	$(AR) $(ARFLAGS) $(MONDO_HOME)/Libs/libMondoMods.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_HOME)/Libs/libMondoMods.a	
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
#
