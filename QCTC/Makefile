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

CPPMISC = -DEXPLICIT_SOURCE

EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions  \
                -I$(MONDO_HOME)/QCTC/MMA/Moments

QCTCObjs=XLate.o CTraX.o Globals.o MondoPoles.o PoleTree.o JGen.o NukE.o QCTC.o  

all:	QCTC

clean:	CQCTC

purge:	CQCTC

mma:	ctrax hgtosp md3trr

ctrax:
	$(MAKE) -i -C MMA/CTraX
hgtosp:
	$(MAKE) -i -C MMA/HGToSP
md3trr:
	$(MAKE) -i -C MMA/MD3TRR

cmma:	
	rm -f HGToSP.Inc
	rm -f CTraX.Inc 
	rm -f MD3TRR.Inc

QCTC:	$(QCTCObjs)
	$(LD) $(LDFLAGS) QCTC $(QCTCObjs) $(LibOpts) $(F90_LIBRARIES) 
	mv QCTC $(MONDO_HOME)/Exec/QCTC

CQCTC:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/QCTC
	ls 

CMD:	
	$(MAKE) -i -C VMD clean

MPoles.o:	HGToSP.Inc

