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
IntB1010101.o\
IntB2010101.o\
IntB2010201.o\
IntB2020101.o\
IntB2020201.o\
IntB2020202.o\
IntB3010101.o\
IntB3010201.o\
IntB3010202.o\
IntB3010301.o\
IntB3020101.o\
IntB3020201.o\
IntB3020202.o\
IntB3020301.o\
IntB3020302.o\
IntB3030101.o\
IntB3030201.o\
IntB3030202.o\
IntB3030301.o\
IntB3030302.o\
IntB3030303.o\
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
IntB1010101.x\
IntB2010101.x\
IntB2010201.x\
IntB2020101.x\
IntB2020201.x\
IntB2020202.x\
IntB3010101.x\
IntB3010201.x\
IntB3010202.x\
IntB3010301.x\
IntB3020101.x\
IntB3020201.x\
IntB3020202.x\
IntB3020301.x\
IntB3020302.x\
IntB3030101.x\
IntB3030201.x\
IntB3030202.x\
IntB3030301.x\
IntB3030302.x\
IntB3030303.x\
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
