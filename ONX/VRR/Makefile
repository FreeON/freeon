include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

EXTRA_INCLUDES = -I$(MONDO_HOME)/ONX/NONX

Objects = Int1111.o Int2111.o Int1121.o Int1122.o Int2211.o Int2121.o \
          Int6111.o Int1161.o Int6611.o Int1166.o Int2122.o Int2221.o \
          Int2222.o Int6161.o Int6211.o Int1162.o Int6121.o Int2161.o \
          Int6221.o Int2162.o Int6122.o Int2261.o

all:	$(Objects)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libVRR.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libVRR.a

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Libs/libVRR.a
