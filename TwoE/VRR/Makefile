include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=-I$(MONDO_HOME)/TwoE/MSC/
#
VRRLibs=$(F90_LIBRARIES) -lVRR $(LAPACK)
#
VRRObjs=VRRd0d0.o\
VRRd0p0.o\
VRRd0s0.o\
VRRp0d0.o\
VRRp0p0.o\
VRRp0s0.o\
VRRs0d0.o\
VRRs0p0.o
#
VRRRels=VRRd0d0.x\
VRRd0p0.x\
VRRd0s0.x\
VRRp0d0.x\
VRRp0p0.x\
VRRp0s0.x\
VRRs0d0.x\
VRRs0p0.x
#
all:    VRR
#
clean:  CVRR
#
purge: clean
	rm -f $(MONDO_LIB)/libVRR.a
	rm -f $(REMOVESRCE)
ppurge: purge
	rm -f VRR*.F90 
#
source:
	math<VRR.m>VRR.out
	rm VRR.out
#
VRR: $(VRRObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libVRR.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libVRR.a
#
CVRR:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
