include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=
#
TwoELibs=$(F90_LIBRARIES) -lTwoE $(LAPACK) 
#
SPObObjs=ShellPairStruct.o
#
TwoEObjs=VScratch.o  \
 dInt1111.o \
 dInt3111.o \
 dInt3311.o \
 dInt3131.o \
 dInt3331.o \
 dInt3333.o 
# Int1111.o \
# Int3111.o \
# Int3311.o \
# Int3131.o \
# Int3331.o \
# Int3333.o \


#
all:    TwoE #N4K dN4K
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

Source: Ints Grads
#
Ints:
	math<Integrals.m
#
Grads:
	math<Gradients.m
#
CTwoE:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
