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
 Int6111.o \
 Int3311.o \
 Int6311.o \
 Int6611.o \
 Int3131.o \
 Int6131.o \
 Int3331.o \
 Int6331.o \
 Int6631.o \
 Int6161.o \
 Int6361.o \
 Int6661.o \
 Int6133.o \
 Int3333.o \
 Int6333.o \
 Int6633.o \
 Int6363.o \
 Int6663.o \
 Int6666.o \
 dInt1111.o \
 dInt3111.o \
 dInt6111.o \
 dInt3311.o \
 dInt6311.o \
 dInt6611.o \
 dInt3131.o \
 dInt6131.o \
 dInt3331.o \
 dInt6331.o \
 dInt6631.o \
 dInt6161.o \
 dInt6361.o \
 dInt6661.o \
 dInt6133.o \
 dInt3333.o \
 dInt6333.o \
 dInt6633.o \
 dInt6363.o \
 dInt6663.o \
 dInt6666.o 


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
