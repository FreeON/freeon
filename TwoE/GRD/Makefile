include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=-I$(MONDO_HOME)/TwoE/MSC/
#
GRDLibs=$(F90_LIBRARIES) -lGRD $(LAPACK)
#
IntObjs=BraHRR11.o\
BraHRR12.o\
BraHRR13.o\
BraHRR21.o\
BraHRR22.o\
BraHRR23.o\
BraHRR31.o\
BraHRR32.o\
BraHRR33.o\
dIntB1010101.o\
dIntB2010101.o\
dIntB2010201.o\
dIntB2020101.o\
dIntB2020201.o\
dIntB2020202.o\
dIntB3010101.o\
dIntB3010201.o\
dIntB3010202.o\
dIntB3010301.o\
dIntB3020101.o\
dIntB3020201.o\
dIntB3020202.o\
dIntB3020301.o\
dIntB3020302.o\
dIntB3030101.o\
dIntB3030201.o\
dIntB3030202.o\
dIntB3030301.o\
dIntB3030302.o\
dIntB3030303.o\
KetHRR21.o\
KetHRR22.o\
KetHRR32.o\
KetHRR33.o
#
IntRels=BraHRR11.x\
BraHRR12.x\
BraHRR13.x\
BraHRR21.x\
BraHRR22.x\
BraHRR23.x\
BraHRR31.x\
BraHRR32.x\
BraHRR33.x\
dIntB1010101.x\
dIntB2010101.x\
dIntB2010201.x\
dIntB2020101.x\
dIntB2020201.x\
dIntB2020202.x\
dIntB3010101.x\
dIntB3010201.x\
dIntB3010202.x\
dIntB3010301.x\
dIntB3020101.x\
dIntB3020201.x\
dIntB3020202.x\
dIntB3020301.x\
dIntB3020302.x\
dIntB3030101.x\
dIntB3030201.x\
dIntB3030202.x\
dIntB3030301.x\
dIntB3030302.x\
dIntB3030303.x\
KetHRR21.x\
KetHRR22.x\
KetHRR32.x\
KetHRR33.x
#
all:    GRD
#
clean:  CGRD
#
purge:clean
	rm -f $(MONDO_LIB)/libGRD.a
	rm -f $(REMOVESRCE)
ppurge:purge
	rm -f *.Inc dIntB*.F90 Ket*.F90 Bra*.F90 VRR*.F90 
#
source:
	math<GRD.m>GRD.out
	rm GRD.out
#
GRD:$(MiscObjs) $(IntObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libGRD.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libGRD.a
#
CGRD:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
