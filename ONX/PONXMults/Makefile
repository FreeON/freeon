include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

FFLAGS = $(FMISC) $(FOPTIMIZE) $(FDEBUG) $(FINCLUDE) -I$(MONDO_HOME)/PONX

objects = md0x0.o md0x1.o md0x2.o md0x3.o md0x4.o md1x0.o md1x1.o md1x2.o \
          md1x3.o md1x4.o md2x0.o md2x1.o md2x2.o md2x3.o md2x4.o md3x0.o \
          md3x1.o md3x2.o md3x3.o md3x4.o md4x0.o md4x1.o md4x2.o md4x3.o \
          md4x4.o ms0x0.o ms0x1.o ms0x2.o ms0x3.o ms0x4.o ms1x0.o ms1x1.o \
          ms1x2.o ms1x3.o ms1x4.o ms2x0.o ms2x1.o ms2x2.o ms2x3.o ms2x4.o \
          ms3x0.o ms3x1.o ms3x2.o ms3x3.o ms3x4.o ms4x0.o ms4x1.o ms4x2.o \
          ms4x3.o ms4x4.o

all:    $(objects)
	$(AR) $(ARFLAGS) Mults.a $(?:.f=o)
	$(RANLIB) Mults.a
	mv Mults.a $(MONDO_HOME)/Libs/libPONXMults.a

clean:
	rm -f $(REMOVEALL)
	rm -f $(MONDO_HOME)/Libs/libPONXMults.a


