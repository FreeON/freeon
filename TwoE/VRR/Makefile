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
VRRd0f0.o\
VRRd0g0.o\
VRRd0h0.o\
VRRd0i0.o\
VRRd0j0.o\
VRRd0p0.o\
VRRd0s0.o\
VRRf0d0.o\
VRRf0f0.o\
VRRf0g0.o\
VRRf0h0.o\
VRRf0i0.o\
VRRf0j0.o\
VRRf0p0.o\
VRRf0s0.o\
VRRg0d0.o\
VRRg0f0.o\
VRRg0g0.o\
VRRg0h0.o\
VRRg0i0.o\
VRRg0j0.o\
VRRg0p0.o\
VRRg0s0.o\
VRRh0d0.o\
VRRh0f0.o\
VRRh0g0.o\
VRRh0h0.o\
VRRh0i0.o\
VRRh0j0.o\
VRRh0p0.o\
VRRh0s0.o\
VRRi0d0.o\
VRRi0f0.o\
VRRi0g0.o\
VRRi0h0.o\
VRRi0i0.o\
VRRi0j0.o\
VRRi0p0.o\
VRRi0s0.o\
VRRj0d0.o\
VRRj0f0.o\
VRRj0g0.o\
VRRj0h0.o\
VRRj0i0.o\
VRRj0j0.o\
VRRj0p0.o\
VRRj0s0.o\
VRRp0d0.o\
VRRp0f0.o\
VRRp0g0.o\
VRRp0h0.o\
VRRp0i0.o\
VRRp0j0.o\
VRRp0p0.o\
VRRp0s0.o\
VRRs0d0.o\
VRRs0f0.o\
VRRs0g0.o\
VRRs0h0.o\
VRRs0i0.o\
VRRs0j0.o\
VRRs0p0.o
#
VRRRels=VRRd0d0.x\
VRRd0f0.x\
VRRd0g0.x\
VRRd0h0.x\
VRRd0i0.x\
VRRd0j0.x\
VRRd0k0.x\
VRRd0p0.x\
VRRd0s0.x\
VRRf0d0.x\
VRRf0f0.x\
VRRf0g0.x\
VRRf0h0.x\
VRRf0i0.x\
VRRf0j0.x\
VRRf0k0.x\
VRRf0p0.x\
VRRf0s0.x\
VRRg0d0.x\
VRRg0f0.x\
VRRg0g0.x\
VRRg0h0.x\
VRRg0i0.x\
VRRg0j0.x\
VRRg0k0.x\
VRRg0p0.x\
VRRg0s0.x\
VRRh0d0.x\
VRRh0f0.x\
VRRh0g0.x\
VRRh0h0.x\
VRRh0i0.x\
VRRh0j0.x\
VRRh0k0.x\
VRRh0p0.x\
VRRh0s0.x\
VRRi0d0.x\
VRRi0f0.x\
VRRi0g0.x\
VRRi0h0.x\
VRRi0i0.x\
VRRi0j0.x\
VRRi0k0.x\
VRRi0p0.x\
VRRi0s0.x\
VRRj0d0.x\
VRRj0f0.x\
VRRj0g0.x\
VRRj0h0.x\
VRRj0i0.x\
VRRj0j0.x\
VRRj0k0.x\
VRRj0p0.x\
VRRj0s0.x\
VRRk0d0.x\
VRRk0f0.x\
VRRk0g0.x\
VRRk0h0.x\
VRRk0i0.x\
VRRk0j0.x\
VRRk0k0.x\
VRRk0p0.x\
VRRk0s0.x\
VRRp0d0.x\
VRRp0f0.x\
VRRp0g0.x\
VRRp0h0.x\
VRRp0i0.x\
VRRp0j0.x\
VRRp0k0.x\
VRRp0p0.x\
VRRp0s0.x\
VRRs0d0.x\
VRRs0f0.x\
VRRs0g0.x\
VRRs0h0.x\
VRRs0i0.x\
VRRs0j0.x\
VRRs0k0.x\
VRRs0p0.x
#
all:    VRR
#
clean:  CVRR
#
release: clean $(ReleaseFiles)
	rm -f *.F90
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
