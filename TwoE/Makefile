include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=
#

TwoELibs=$(F90_LIBRARIES) -lTwoE $(LAPACK) 

SPObObjs=ShellPairStruct.o

N4KObjs=GetAtomPair.o N4K.o

dN4KObjs=GetAtomPair.o dN4K.o
#
TwoEObjs=VScratch.o  \
 Int1111.o \
 Int3111.o \
 Int1311.o \
 Int3311.o \
 Int1131.o \
 Int3131.o \
 Int1331.o \
 Int3331.o \
 Int1113.o \
 Int3113.o \
 Int1313.o \
 Int3313.o \
 Int1133.o \
 Int3133.o \
 Int1333.o \
 Int3333.o \
 dInt1111.o \
 dInt3111.o \
 dInt1311.o \
 dInt3311.o \
 dInt1131.o \
 dInt3131.o \
 dInt1331.o \
 dInt3331.o \
 dInt1113.o \
 dInt3113.o \
 dInt1313.o \
 dInt3313.o \
 dInt1133.o \
 dInt3133.o \
 dInt1333.o \
 dInt3333.o 

#
all:    TwoE N4K dN4K
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

N4K: TwoE $(N4KObjs) ERIInclude.Inc
	$(LD) $(LDFLAGS) N4K $(N4KObjs)  $(TwoELibs)
	mv N4K $(MONDO_EXEC)/N4K

dN4K: TwoE $(dN4KObjs) dERIInclude.Inc
	$(LD) $(LDFLAGS) dN4K $(dN4KObjs)  $(TwoELibs)
	mv dN4K $(MONDO_EXEC)/dN4K

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
