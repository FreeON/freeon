include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=
#
SPObObjs=ShellPairStruct.o
#
TwoEObjs=  \
	Int1111.o \
	Int3111.o \
	Int3311.o \
	Int3131.o \
	Int3331.o \
	Int3333.o 
#
TwoERels=  \
	Int1111.x \
	Int3111.x \
	Int3311.x \
	Int3131.x \
	Int3331.x \
	Int3333.x 
#
all:    TwoE
#
clean:  CTwoE
#
purge:clean 
	rm -f $(MONDO_LIB)/libTwoE.a
	rm -f $(REMOVESRCE)
#
TwoE:$(SPObObjs) $(TwoEObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libTwoE.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libTwoE.a
#
CTwoE:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
