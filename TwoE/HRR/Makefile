include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC =
#
EXTRA_INCLUDES=-I$(MONDO_HOME)/TwoE/MSC/
#
HRRLibs=$(F90_LIBRARIES) -lHRR $(LAPACK)
#
HRRObjs=KetHRR1010.o\
KetHRR101.o\
KetHRR102.o\
KetHRR103.o\
KetHRR104.o\
KetHRR106.o\
KetHRR1510.o\
KetHRR1515.o\
KetHRR151.o\
KetHRR152.o\
KetHRR153.o\
KetHRR154.o\
KetHRR156.o\
KetHRR21.o\
KetHRR22.o\
KetHRR32.o\
KetHRR33.o\
KetHRR41.o\
KetHRR42.o\
KetHRR43.o\
KetHRR44.o\
KetHRR61.o\
KetHRR62.o\
KetHRR63.o\
KetHRR64.o\
KetHRR66.o
#
HRRRels=KetHRR1010.x\
KetHRR101.x\
KetHRR102.x\
KetHRR103.x\
KetHRR104.x\
KetHRR106.x\
KetHRR1510.x\
KetHRR1515.x\
KetHRR151.x\
KetHRR152.x\
KetHRR153.x\
KetHRR154.x\
KetHRR156.x\
KetHRR21.x\
KetHRR22.x\
KetHRR32.x\
KetHRR33.x\
KetHRR41.x\
KetHRR42.x\
KetHRR43.x\
KetHRR44.x\
KetHRR61.x\
KetHRR62.x\
KetHRR63.x\
KetHRR64.x\
KetHRR66.x
#
all:    HRR
#
clean:  CHRR
#
purge: clean
	rm -f $(MONDO_LIB)/libHRR.a
	rm -f $(REMOVESRCE)
ppurge: purge
	rm -f HRR*.F90 
#
source:
	math<HRR.m>HRR.out
	rm HRR.out
#
HRR: $(HRRObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libHRR.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libHRR.a
#
CHRR:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
