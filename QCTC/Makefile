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
CPPMISC = -DEXPLICIT_SOURCE
#
EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions  \
                -I$(MONDO_HOME)/QCTC/MMA/Moments
#
QCTCObjs=XLate.o CTraX.o Globals.o MondoPoles.o PoleTree.o JGen.o NukE.o QCTC.o
GObjs=XLate.o CTraX.o Globals.o MondoPoles.o PoleTree.o GradJGen.o GradNukE.o GradQCTC.o
#
all:	QCTC GradQCTC
#
clean:	CQCTC cmma CGQCTC
#
purge:	CQCTC cmma CGQCTC
#
QCTC:	$(QCTCObjs)
	$(LD) $(LDFLAGS) QCTC $(QCTCObjs) $(LibOpts) $(F90_LIBRARIES) 
	mv QCTC $(MONDO_HOME)/Exec/QCTC
#
GradQCTC:	$(GObjs)
	$(LD) $(LDFLAGS) GradQCTC $(GObjs) $(LibOpts) $(F90_LIBRARIES) 
	mv GradQCTC $(MONDO_HOME)/Exec/GradQCTC
#
CQCTC:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/QCTC
	ls 
#
#
CGQCTC:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/GradQCTC
	ls 
#
mma:	ctrax hgtosp pac md3trr
#
cmma:	cctrax chgtosp cpac cmd3trr 
#
pmma:	cmma
	rm -f HGToSP.Inc
	rm -f CTraX.Inc 
	rm -f PFunk.Inc
	rm -f MD3TRR.Inc
#
ctrax:
	$(MAKE) -i -C MMA/CTraX
cctrax:
	$(MAKE) -i -C MMA/CTraX  clean
#	
hgtosp:
	$(MAKE) -i -C MMA/HGToSP
chgtosp:
	$(MAKE) -i -C MMA/HGToSP clean
#	
pac:
	$(MAKE) -i -C MMA/PAC
cpac:
	$(MAKE) -i -C MMA/PAC    clean
#
md3trr:
	$(MAKE) -i -C MMA/MD3TRR
cmd3trr:
	$(MAKE) -i -C MMA/MD3TRR clean
#
MPoles.o:	HGToSP.Inc

