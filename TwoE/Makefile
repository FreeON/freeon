include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=-I$(MONDO_HOME)/TwoE/MSC/ -I$(MONDO_HOME)/TwoE/INT/ -I$(MONDO_HOME)/TwoE/GRD/ 
#
TwoELibs=$(F90_LIBRARIES) -lINT -lVRR -lHRR -lTwoEMisc
#
all:    Misc HRRs VRRs Ints Grad
#
clean:  cMisc cHRRs cVRRs cInts cGrad
	rm -f N4K dN4K
#
source: sHRRs sVRRs sInts sGrad
#
purge:  clean pMisc pHRRs pVRRs pInts pGrad
#
ppurge: purge ppMain ppHRRs ppVRRs ppInts ppGrad 
#
Misc:	
	$(MAKE) -C MSC
cMisc:	
	$(MAKE) -C MSC clean
pMisc:
	$(MAKE) -C MSC purge
ppVRRs:
	$(MAKE) -C MSC ppurge
#
VRRs:	
	$(MAKE) -C VRR
cVRRs:	
	$(MAKE) -C VRR clean
pVRRs:
	$(MAKE) -C VRR purge
ppVRRs:
	$(MAKE) -C VRR ppurge
sVRRs:
	$(MAKE) -C VRR source
#
HRRs:	
	$(MAKE) -C HRR
cHRRs:	
	$(MAKE) -C HRR clean
pHRRs:
	$(MAKE) -C HRR purge
ppHRRs:
	$(MAKE) -C HRR ppurge
sHRRs:
	$(MAKE) -C HRR source
#
Ints:	
	$(MAKE) -C INT
cInts:	
	$(MAKE) -C INT clean
pInts:	
	$(MAKE) -C INT purge
ppInts:	
	$(MAKE) -C INT ppurge
sInts:	
	$(MAKE) -C INT source
#
Grad:	
	$(MAKE) -C GRD
cGrad:	
	$(MAKE) -C GRD clean
pGrad:	
	$(MAKE) -C GRD purge
ppGrad:	
	$(MAKE) -C GRD ppurge
sGrad:	
	$(MAKE) -C GRD source
#
N4K: Misc HRRs VRRs Ints N4K.o
	$(LD) $(LDFLAGS) N4K N4K.o $(TwoELibs)
#
dN4K: Misc HRRs VRRs Ints Grad dN4K.o
	$(LD) $(LDFLAGS) dN4K dN4K.o -lGRD $(TwoELibs) 
#
ppMain:
	rm -f dInt*.F90 Int*.F90

N4K.o:	INT/ERIInterfaceB.Inc

dN4K.o:	INT/ERIInterfaceB.Inc GRD/dERIInterfaceB.Inc
